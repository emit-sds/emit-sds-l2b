using Logging: @info
using LinearAlgebra: transpose
using Statistics: mean
using ArgParse2
using Logging
using ArchGDAL
using CSV
using DataFrames

function write_envi_header(output_file, samples, lines, bands, band_names, map_info=nothing, css=nothing)
    output_header =
"ENVI
samples = $(samples)
lines = $(lines)
bands = $(bands)
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bil
data ignore value = -9999
byte order = 0
band names = { $(join(band_names," , ")) }
"
    if !isnothing(map_info)
        output_header *= "\nmap_info = $(map_info)"
    end
    if !isnothing(css)
        output_header *= "\ncoordinate_system_string = $(css)"
    end
    open(output_file, "w") do fout
        write(fout, output_header)
    end
end

function filter_by_names(in_mat, mineral_names)
    out_mat = Matrix(coalesce.(in_mat[!,mineral_names], nothing))
    out_mat[isnothing.(out_mat)] .= 0
    out_mat[out_mat .== -1] .= 1
    return out_mat
end

function write_output(outdat, outfile, outuncdat, outuncfile)
    # write as BIL interleave
    @info "Writing output file $outfile of shape $(size(outdat))"
    outdat = permutedims(outdat, (2, 3, 1))
    open(outfile, "w") do fout
        write(fout, reinterpret(Float32, outdat))
    end

    if !isnothing(outuncdat)
        # write as BIL interleave
        outuncdat = permutedims(outuncdat, (2, 3, 1))
        open(outuncfile, "w") do fout
            write(fout, reinterpret(Float32, outuncdat))
        end
    end
end

function mean_optical_path_length(absorption_coefficient::Float64, band_depth::Float64)
    BD = 1 .- band_depth
    delta_k = absorption_coefficient
    MOPL = max.(log.(BD) ./ (-1.0 .* delta_k), 1e-12) # Mean Optical Path Length.  T_slab = MOPL/2
    return MOPL
end


function model_3(mineral_groupings::Matrix, abundance_metadata::DataFrame, band_depth::Vector)

    abs_abundance = zeros(length(abundance_metadata.mineral))
    m_return = zeros(length(abundance_metadata.mineral))
    T_return = zeros(length(abundance_metadata.mineral))
    if band_depth[2] > 0
        present_in = findall(x -> x > 0, mineral_groupings[Int(band_depth[2]), :])
        for mineral in present_in
            MOPL = mean_optical_path_length(abundance_metadata.delta_k_absorption_coefficient[mineral], band_depth[1])
            D = abundance_metadata.grain_diameter_microns[mineral] / 10000.0
            m = MOPL / (2 * D)
            T = pi/4 * D * sqrt(m/2)
            abs_abundance[mineral] = 3.0 * T * 100.0^2 * abundance_metadata.density_g_cc[mineral]
            m_return[mineral] = m
            T_return[mineral] = T
        end
    end

    if band_depth[4] > 0
        present_in = findall(x -> x > 0, mineral_groupings[Int(band_depth[4]), :])
        for mineral in present_in
            MOPL = mean_optical_path_length(abundance_metadata.delta_k_absorption_coefficient[mineral], band_depth[3])
            D = abundance_metadata.grain_diameter_microns[mineral] / 10000.0
            m = MOPL / (2 * D)
            T = pi/4 * D * sqrt(m/2)
            abs_abundance[mineral] = 3.0 * T * 100.0^2 * abundance_metadata.density_g_cc[mineral]
            m_return[mineral] = m
            T_return[mineral] = T
        end
    end

    mean_m = mean(m_return[findall(x -> x != 0, m_return)])
    if isnan(mean_m)
        mean_m = 0
    end

    return abs_abundance, mean_m, m_return, T_return

end

function model_4_external_qf(m3::Vector, mineral_groupings::Matrix, abundance_metadata::DataFrame, band_depth::Vector)

    rel_abundance = zeros(length(abundance_metadata.mineral))
    if band_depth[2] > 0
        present_in = findall(x -> x > 0, mineral_groupings[Int(band_depth[2]), :])
        for mineral in present_in
            rel_abundance[mineral] = m3[mineral]^0.5 * abundance_metadata.density_g_cc[mineral]
        end
    end

    if band_depth[4] > 0
        present_in = findall(x -> x > 0, mineral_groupings[Int(band_depth[4]), :])
        for mineral in present_in
            rel_abundance[mineral] = m3[mineral]^0.5 * abundance_metadata.density_g_cc[mineral]
        end
    end

    relsum = sum(rel_abundance)
    if sum(rel_abundance) > 0
        rel_abundance = rel_abundance ./ relsum
    end 

    return rel_abundance

end


function model_4_plus(m::Vector, mineral_groupings::Matrix, abundance_metadata::DataFrame, band_depth::Vector, q_dens::Float64, q_grainsize::Float64)

    rel_abundance = zeros(length(abundance_metadata.mineral))
    present_in = []
    if band_depth[2] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[2]), :]))
    end
    if band_depth[4] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[4]), :]))
    end

    if length(present_in) == 0
        return rel_abundance
    end
    
    den_sum = 0
    for mineral in present_in
        D_i = abundance_metadata.grain_diameter_microns[mineral]/10000 # cm
        row_i = abundance_metadata.density_g_cc[mineral]
        den_sum += D_i^1.5 * m[mineral]^0.25 * row_i
    end

    m3m = mean(m[findall(x -> x != 0, m)])
    D_q = q_grainsize / 10000
 
    for mineral in present_in
        row_i = abundance_metadata.density_g_cc[mineral]
        D_i = abundance_metadata.grain_diameter_microns[mineral]/10000 # cm
        rel_abundance[mineral] = D_i^1.5 * m[mineral]^0.25 * row_i / (den_sum + D_q^1.5 * m3m^0.25 * q_dens)
    end

    rel_abundance_qf = D_q^1.5 * m3m^0.25 * q_dens / (den_sum + D_q^1.5 * m3m^0.25 * q_dens)

    return rel_abundance, rel_abundance_qf

end


function model_4_bulk(m::Vector, mineral_groupings::Matrix, abundance_metadata::DataFrame, band_depth::Vector, bulk_dens_g_cc::Float64, bulk_grainsize_micron::Float64)

    rel_abundance = zeros(length(abundance_metadata.mineral))
    present_in = []
    if band_depth[2] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[2]), :]))
    end
    if band_depth[4] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[4]), :]))
    end

    if length(present_in) == 0
        return rel_abundance
    end
    
    D_p = bulk_grainsize_micron / 10000
    

    for mineral in present_in
        row_i = abundance_metadata.density_g_cc[mineral]
        D_i = abundance_metadata.grain_diameter_microns[mineral]/10000 # cm
        x_i = m[mineral] * 2 * D_i # back to MOPL

        #println(D_i ^ 0.5 * pi, " ", 8 * x_i ^ 0.5)
        rel_abundance[mineral] = row_i * D_i^2 * (D_i ^ 0.5 * pi - 4*sqrt(2) * x_i ^ 0.5) / 
                                 (bulk_dens_g_cc * D_p^2 * 4*sqrt(2) * x_i^0.5 + row_i * D_i ^ (5/2) * pi - 4*sqrt(2) * row_i * D_i^2 * x_i ^ 0.5)

    end

    # shouldn't ever get here, but just in case
    #relsum = sum(rel_abundance)
    #if sum(rel_abundance) > 0
    #    rel_abundance = rel_abundance ./ relsum
    #end 

    return rel_abundance, 1 - sum(rel_abundance)

end


function model_4_bulk_p(m::Vector, mineral_groupings::Matrix, abundance_metadata::DataFrame, band_depth::Vector, bulk_dens_g_cc::Float64, bulk_grainsize_micron::Float64)

    rel_abundance = zeros(length(abundance_metadata.mineral))
    present_in = []
    if band_depth[2] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[2]), :]))
    end
    if band_depth[4] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[4]), :]))
    end

    if length(present_in) == 0
        return rel_abundance
    end
    
    D_p = bulk_grainsize_micron / 10000
    

    for mineral in present_in
        row_i = abundance_metadata.density_g_cc[mineral]
        D_i = abundance_metadata.grain_diameter_microns[mineral]/10000 # cm
        x_i = m[mineral] * 2 * D_i # back to MOPL

        rel_abundance[mineral] = 4 * 2^(1/2) * x_i^(1/2) * row_i * D_i^(3/2) / (pi * bulk_dens_g_cc * D_p^2) 
    end

    # shouldn't ever get here, but just in case
    #relsum = sum(rel_abundance)
    #if sum(rel_abundance) > 0
    #    rel_abundance = rel_abundance ./ relsum
    #end 

    return rel_abundance, 1 - sum(rel_abundance)

end

function model_4_bulk_non(m::Vector, mineral_groupings::Matrix, abundance_metadata::DataFrame, band_depth::Vector, bulk_dens_g_cc::Float64, bulk_grainsize_micron::Float64)

    rel_abundance = zeros(length(abundance_metadata.mineral))
    present_in = []
    if band_depth[2] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[2]), :]))
    end
    if band_depth[4] > 0
        append!(present_in, findall(x -> x > 0, mineral_groupings[Int(band_depth[4]), :]))
    end

    if length(present_in) == 0
        return rel_abundance
    end
    
    D_p = bulk_grainsize_micron / 10000
    

    for mineral in present_in
        row_i = abundance_metadata.density_g_cc[mineral]
        D_i = abundance_metadata.grain_diameter_microns[mineral]/10000 # cm
        x_i = m[mineral] * 2 * D_i # back to MOPL

        rel_abundance[mineral] = 4 * 2^(1/2) * x_i^(1/2) * D_i^(1/2) * row_i / (pi * bulk_dens_g_cc * D_p) 
    end

    # shouldn't ever get here, but just in case
    #relsum = sum(rel_abundance)
    #if sum(rel_abundance) > 0
    #    rel_abundance = rel_abundance ./ relsum
    #end 

    return rel_abundance, 1 - sum(rel_abundance)

end




function main()
    parser = ArgumentParser(description="Translate to Rrs. and/or apply masks")
    add_argument!(parser, "output_base", type=String, metavar="OUTPUT")
    add_argument!(parser, "band_depth_file", type=String, metavar="Band Depth file.  4 bands (G1 BD, G1 Ref, G2 BD, G2 Ref)")
    add_argument!(parser, "quartz_size_microns", type=Float64)
    add_argument!(parser, "average_mingrains", type=Float64)
    add_argument!(parser, "--quartz_size_file", type=String, default=nothing)
    add_argument!(parser, "--quartz_file_scaling_factor", type=Float64, default=1.0)
    add_argument!(parser, "--quartz_uncertainty_delta", type=Float64, default=0.0)
    add_argument!(parser, "--quartz_gs_lb", type=Float64, default=0.0)
    add_argument!(parser, "--quartz_massfrac_file", type=String, default=nothing)
    add_argument!(parser, "--band_depth_unc_file", type=String, default=nothing)
    add_argument!(parser, "--model_style", type=String, default="m3")
    add_argument!(parser, "--q_dens", type=Float64, default=2.65)
    add_argument!(parser, "--mineral_groupings_matrix", type=String, default="data/mineral_grouping_matrix_20231113.csv")
    add_argument!(parser, "--abundance_metadata", type=String, default="data/abundance_metadata_20240123.csv")
    add_argument!(parser, "--sum_grainsize", action="store_true")
    add_argument!(parser, "--output_m", action="store_true")
    add_argument!(parser, "--log_file", type=String, default=nothing)
    add_argument!(parser, "--log_level", type=String, default="INFO")
    args = parse_args(parser)

    if isnothing(args.log_file)
        logger = Logging.SimpleLogger()
    else
        logger = Logging.SimpleLogger(open(args.log_file, "w+"))
    end

    band_depth_ds = ArchGDAL.read(args.band_depth_file)
    band_depth = permutedims(convert(Array{Float64}, ArchGDAL.readraster(args.band_depth_file)), (2,1,3))
    band_depth_header = readlines(args.band_depth_file * ".hdr")
    header_keys = [strip(split(x,"=")[1]) for x in band_depth_header]
    map_info_idx = findfirst(k->occursin("map info",k), [split(x,"=")[1] for x in band_depth_header])
    map_info = nothing
    if !isnothing(map_info_idx)
        map_info = split(band_depth_header[map_info_idx], "=")[2]
    end

    css_idx = findfirst(k->occursin("coordinate system string",k), [split(x,"=")[1] for x in band_depth_header])
    css = nothing
    if !isnothing(css_idx)
        css = split(band_depth_header[css_idx], "=")[2]
    end


    


    if !isnothing(args.band_depth_unc_file)
        band_depth_unc_ds = ArchGDAL.read(args.band_depth_unc_file)
        band_depth_unc = permutedims(convert(Array{Float64}, ArchGDAL.readraster(args.band_depth_unc_file)), (2,1,3))
    else
        band_depth_unc = nothing
    end

    if !isnothing(args.quartz_massfrac_file) && !isnothing(args.quartz_size_file)
        @error "Can only have one of quartz_massfrac_file or quartz_size_file"
        exit(1)
    end
    
    if args.model_style == "m4" && isnothing(args.quartz_massfrac_file) 
        @error "Current M4 implementation requires quartz_massfrac_file"
        exit(1)
    end

    quartz_massfrac = nothing
    if !isnothing(args.quartz_massfrac_file)
        quartz_massfrac_ds = ArchGDAL.read(args.quartz_massfrac_file)
        quartz_massfrac = permutedims(convert(Array{Float64}, ArchGDAL.readraster(args.quartz_massfrac_file))[:,:,1], (2,1)) / 100.
    end
    
    quartz_size = nothing
    if !isnothing(args.quartz_size_file)
        quartz_size_ds = ArchGDAL.read(args.quartz_size_file)
        quartz_size = (permutedims(convert(Array{Float64}, ArchGDAL.readraster(args.quartz_size_file))[:,:,end], (2,1))  .+ args.quartz_uncertainty_delta ) .* args.quartz_file_scaling_factor
    end

    abundance_metadata = CSV.read(args.abundance_metadata, DataFrame)
    mineral_groupings = CSV.read(args.mineral_groupings_matrix, DataFrame)

    # get header form mineral_groupings
    mi_header = names(mineral_groupings)
    mineral_names = [x for (_x, x) in enumerate(mi_header) if _x >= findall(mi_header .== "Calcite")[1] && _x <= findall(mi_header .== "Vermiculite")[1]]
    mineral_groupings = filter_by_names(mineral_groupings, mineral_names)

    num_minerals = length(mineral_names) + 1

    (rows, cols, __) = size(band_depth)

    @info "Loading complete, set up output file(s)"

    out_file = args.output_base * "_abs_abundance"
    out_file_rel = args.output_base * "_rel_abundance"
    out_unc_file = args.output_base * "_mfabundance_unc"
    out_scat_file = args.output_base * "_scatter"
    out_light_file = args.output_base * "_light"

    out_mineral_names = vcat(mineral_names, ["Quartz+Feldspar"])
    write_envi_header(out_file * ".hdr", cols, rows, num_minerals, out_mineral_names, map_info, css)
    write_envi_header(out_file_rel * ".hdr", cols, rows, num_minerals, out_mineral_names, map_info, css)
    if args.output_m
        write_envi_header(out_scat_file * ".hdr", cols, rows, num_minerals, out_mineral_names, map_info, css)
        write_envi_header(out_light_file * ".hdr", cols, rows, num_minerals, out_mineral_names, map_info, css)
    end
    if !isnothing(args.band_depth_unc_file)
        write_envi_header(out_unc_file * ".hdr", cols, rows, num_minerals, out_mineral_names, map_info, css)
    end
    out_data = zeros(Float32, rows, cols, num_minerals)
    if args.output_m
        out_m = zeros(Float32, rows, cols, num_minerals)
        out_T = zeros(Float32, rows, cols, num_minerals)
    end

    out_unc = nothing
    if !isnothing(args.band_depth_unc_file)
        out_unc = zeros(Float32, rows, cols, num_minerals)
    end


    Threads.@threads for _y in 1:rows
        println(_y, " ", size(band_depth)[1])
        for _x in 1:cols

            (m3, m3m, m_return, T_return) = model_3(mineral_groupings, abundance_metadata, band_depth[_y, _x, :])
            if args.output_m
                out_m[_y, _x, 1:end-1] .= m_return
                out_T[_y, _x, 1:end-1] .= T_return
                out_m[_y, _x, end] = m3m
            end

            if isnothing(args.quartz_size_file)
                qs = args.quartz_size_microns
            else
                qs = quartz_size[_y, _x]
            end

            if qs < args.quartz_gs_lb
                qs = args.quartz_gs_lb
            end

            if args.model_style == "m3"
                out_data[_y, _x, 1:end-1] = m3

                if !isnothing(args.quartz_massfrac_file)
                    non_quartz_mass = sum(m3)
                    quartz_mass = non_quartz_mass  * quartz_massfrac[_y, _x] / (1 - quartz_massfrac[_y, _x])
                    out_data[_y,_x,end] = quartz_mass
                else

                    T_q = pi/4 * (qs / 10000.0) / args.average_mingrains * sqrt(m3m / 2)
                    qmass = 3 * T_q * 100.0^2 * args.q_dens
                    out_data[_y,_x,end] = qmass

                end


            elseif args.model_style == "m4"
                if !isnothing(args.quartz_massfrac_file)
                    m4 = model_4_external_qf(m3, mineral_groupings, abundance_metadata, band_depth[_y, _x, :])
                    out_data[_y, _x, 1:end-1] .= m4 * (1 - quartz_massfrac[_y, _x])
                    out_data[_y, _x, end] = quartz_massfrac[_y, _x]
                else
                    m4 = model_4(m3, mineral_groupings, abundance_metadata, band_depth[_y, _x, :])
                    qm = sum(m4)
                    out_data[_y, _x, 1:end-1] .= m4 
                    out_data[_y, _x, end] = qm
                end
            elseif args.model_style == "m4+"
                m4, m4q = model_4_plus(m_return, mineral_groupings, abundance_metadata, band_depth[_y, _x, :], args.q_dens, qs)
                out_data[_y, _x, 1:end-1] .= m4 
                out_data[_y, _x, end] = m4q
            elseif args.model_style == "m4b"
                m4, m4q = model_4_bulk(m_return, mineral_groupings, abundance_metadata, band_depth[_y, _x, :], args.q_dens, qs)
                out_data[_y, _x, 1:end-1] .= m4 
                out_data[_y, _x, end] = m4q
            elseif args.model_style == "m4bp"
                m4, m4q = model_4_bulk_p(m_return, mineral_groupings, abundance_metadata, band_depth[_y, _x, :], args.q_dens, qs)
                out_data[_y, _x, 1:end-1] .= m4
                out_data[_y, _x, end] = m4q
            elseif args.model_style == "m4bn"
                m4, m4q = model_4_bulk_non(m_return, mineral_groupings, abundance_metadata, band_depth[_y, _x, :], args.q_dens, qs)
                out_data[_y, _x, 1:end-1] .= m4
                out_data[_y, _x, end] = m4q
            end

            

        end
    end

    if args.model_style == "m3"
        total_abundance = sum(out_data, dims=3)[:,:,1]
        out_data[total_abundance .<= 0, :] .= -9999.0
        write_output(out_data, out_file, out_unc, out_unc_file)
        for _y in 1:rows
            for _x in 1:cols
                if total_abundance[_y,_x] > 0
                    out_data[_y,_x,:] ./= total_abundance[_y,_x]
                end
            end
        end
        write_output(out_data, out_file_rel, nothing, nothing)
    elseif args.model_style in ["m4", "m4+", "m4b", "m4bp", "m4bn"]
        write_output(out_data, out_file_rel, nothing, nothing)
    end
    if args.output_m
        write_output(out_m, out_scat_file, nothing, nothing)
        write_output(out_T, out_light_file, nothing, nothing)
    end
 

end


#main()