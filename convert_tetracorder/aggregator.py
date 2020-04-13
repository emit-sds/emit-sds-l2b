# David Thompson


import argparse
import gzip
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import spectral.io.envi as envi
import emit_utils.file_checks
import gdal
import logging
import emit_utils.common_logs


def calculate_band_depth(wavelengths: np.array, reflectance: np.array, feature: tuple):
    """ Calculate the Clark, 2003 continuum normalized band depth of a particular feature.
    Args:
        wavelengths: an array of wavelengths corresponding to given reflectance values
        reflectance: an array of reflectance values to calculate the band depth from
        feature: definition of the feature to calculate band depth for, with the first two and last two values defining 
                 the averaging windows used identify the feature of interest. 
    :Returns
        band_depth: the band depth as defined in Clark, 2003.
    """

    left_inds = np.where(np.logical_and(wavelengths >= feature[0], wavelengths <= feature[1]))[0]
    left_x = wavelengths[int(left_inds.mean())]
    left_y = reflectance[left_inds].mean()

    right_inds = np.where(np.logical_and(wavelengths >= feature[2], wavelengths <= feature[3]))[0]
    right_x = wavelengths[int(right_inds.mean())]
    right_y = reflectance[right_inds].mean()

    continuum = interp1d([left_x, right_x], [left_y, right_y],
                         bounds_error=False, fill_value='extrapolate')(wavelengths)
    feature_inds = np.logical_and(wavelengths >= feature[0], wavelengths <= feature[3])

    # Band Depth definition from Clark, 2003 - max over this range will be taken as the 'band depth'
    depths = 1.0 - reflectance[feature_inds] / continuum[feature_inds]

    return max(depths)


def calculate_uncertainty(wavelengths: np.array, observed_reflectance: np.array, observed_reflectance_uncertainty: np.array,
                library_reflectance: np.array, feature: tuple):
    """ Calculate the uncertainty of the Clark, 2003 continuum normalized band depth of a particular feature.
    Args:
        wavelengths: an array of wavelengths corresponding to given reflectance values
        observed_reflectance: an array of observed reflectance values
        observed_reflectance_uncertainty: an array of uncertainties of the observed reflectance values
        library_reflectance: an array of library reference reflectance values
        feature: definition of the feature to calculate band depth for, with the first two and last two values defining
                 the averaging windows used identify the feature of interest.
    :Returns
        band_depth: the uncertainty for the band depth as defined in Clark, 2003.
    """
    inds = np.where(np.logical_and(wavelengths >= feature[1], wavelengths <= feature[2]))[0]

    num = (len(inds) * np.sum(observed_reflectance[:, inds, :]*library_reflectance[np.newaxis,inds,np.newaxis], axis=1) - \
            np.sum(observed_reflectance[:, inds, :], axis=1)*np.sum(library_reflectance[np.newaxis,inds,np.newaxis]))
    den = len(inds) * np.sum(np.power(library_reflectance[inds],2)) - np.power(np.sum(library_reflectance[inds]),2)

    a = num/den

    left_term = np.power(1 / (np.sum(np.power(a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis],2),axis=1) - np.power(np.sum(a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis],axis=1),2)),2)

    ## looped way, depcreated to below vectorization, preserved for context only
    #right_term = 0
    #for w in inds:
    #    right_term += np.power(len(inds)*a*library_reflectance[w] - np.sum(a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis],axis=1),2) * np.power(observed_reflectance_uncertainty[:,w,:],2)

    right_term2 = np.sum(np.power(len(inds)*a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis] - np.sum(a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis],axis=1)[:,np.newaxis,:],2) * np.power(observed_reflectance_uncertainty[:,inds,:],2),axis=1)

    psi = np.sqrt(left_term*right_term)

    return psi


# parse the command line (perform the correction on all command line arguments)
def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('tetra_expert_file', type=str, metavar='TETRA_EXPERT_SYSTEM')
    parser.add_argument('tetra_library_file', type=str, metavar='TETRA_LIBRARY_FILE')
    parser.add_argument('tetra_output_base', type=str, metavar='TETRA_OUTPUT_DIR')
    parser.add_argument('translation_file', type=str, metavar='MINERAL_FRACTIONS_FILE')
    parser.add_argument('output', type=str, metavar='OUTPUT')
    parser.add_argument('-calculate_uncertainty', type=int, choices=[0,1], metavar='CALCULATE_UNCERTAINTY')
    parser.add_argument('-reflectance_file', type=str, metavar='REFLECTANCE_FILE')
    parser.add_argument('-reflectance_uncertainty_file', type=str, metavar='REFLECTANCE_UNCERTAINTY_FILE')
    parser.add_argument('-log_file', type=str, default=None)
    parser.add_argument('-log_level', type=str, default='INFO')
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    emit_utils.common_logs.logtime()

    library = envi.open(args.tetra_library_file+'.hdr', args.tetra_library_file)
    library_reflectance = library.spectra.copy()
    library_records = [int(q) for q in library.metadata['record']]
    names = [q.strip() for q in library.metadata['spectra names']]
    wavelengths = np.array([float(q) for q in library.metadata['wavelength']])

    num_minerals = 10
    out_data = None

    # Read the EMIT mineral fraction csv
    translation_file_df = pd.read_csv(args.translation_file)
    emit_band_names = np.array(list(translation_file_df))[2:].tolist()
    tetra_record_numbers = np.array(translation_file_df['#Record']).tolist()
    emit_10_mixture_fractions = np.array(translation_file_df[emit_band_names])

    if args.calculate_uncertainty == 1:
        args.calculate_uncertainty = True
        emit_utils.file_checks.check_raster_files([args.reflectance_file, args.reflectance_uncertainty_file], map_space=False)
        refl_dataset = gdal.Open(args.reflectance_file, gdal.GA_ReadOnly)
        observed_reflectance = np.memmap(args.reflectance_file, mode='r',
                                         shape=(refl_dataset.RasterYSize, refl_dataset.RasterCount,
                                                refl_dataset.RasterXSize), dtype=np.float32)
        observed_reflectance_uncertainty = np.memmap(args.reflectance_uncertainty_file, mode='r',
                                         shape=(refl_dataset.RasterYSize, refl_dataset.RasterCount,
                                                refl_dataset.RasterXSize), dtype=np.float32)
    else:
        args.calculate_uncertainty = False


    # read expert system file and strip comments
    with open(args.tetra_expert_file, 'r') as fin:
        expert_file_commented = fin.readlines()

    expert_file_text, orig_lineno = [], []
    for line_index, line in enumerate(expert_file_commented):
        if not line.strip().startswith('\#'):
            orig_lineno.append(line_index)
            expert_file_text.append(line)
    del expert_file_commented

    # Go through expert system file one line at a time, after initializing key variables
    expert_line_index, group, spectrum, output_data, header, out_hdr, rows, cols = \
        0, None, None, None, True, None, 0, 0
    while expert_line_index < len(expert_file_text):

        # The Header flag excludes the definitions at the start
        if expert_file_text[expert_line_index].startswith('BEGIN SETUP'):
            header = False
        elif header:
            expert_line_index = expert_line_index + 1
            continue

        # if keyword 'group' appears, define the current group name
        if expert_file_text[expert_line_index].startswith('group'):
            group = int(expert_file_text[expert_line_index].strip().split()[1])

        # if we've gotten to the end of the record, time to pull everything together and write our output
        if expert_file_text[expert_line_index].startswith('endaction'):
            if group in [1, 2]:

                # get band depth of spectrum library file
                try:
                    library_band_depth = calculate_band_depth(
                        wavelengths, library_reflectance[library_records.index(record), :], features[0])

                except ValueError:
                    expert_line_index = expert_line_index + 1
                    continue

                # get band depths from tetracorder output map
                groupdir = args.tetra_output_base + '/group.' + str(group) + 'um/'
                hdrpath = groupdir + filename + '.depth.gz.hdr'
                datapath = groupdir + filename + '.depth.gz'
                logging.info('loading: {}'.format(hdrpath))
                try:
                    # TODO update IO
                    # read header, arrange output data files
                    hdr = envi.read_envi_header(hdrpath)
                    offs = int(hdr['header offset'])
                    if out_data is None:
                        out_hdr = hdr.copy()
                        if ('file compression' in out_hdr.keys()):
                            del out_hdr['file compression']
                        out_hdr['interleave'] = 'bil'
                        out_hdr['data type'] = 4
                        out_hdr['wavelengths'] = '{'+','.join([str(q) for q in wavelengths])+'}'
                        out_hdr['bands'] = num_minerals
                        out_hdr['header offset'] = 0
                        out_hdr['band names'] = emit_band_names
                        cols = int(hdr['samples'])
                        rows = int(hdr['lines'])
                        out_data = np.zeros((rows, cols, num_minerals), dtype=np.float32)
                        if args.calculate_uncertainty:
                            out_uncertainty = np.zeros((rows, cols, num_minerals), dtype=np.float32)

                    # read band depth
                    with open(datapath, 'rb') as fin:
                        compressed = fin.read()
                    decompressed = gzip.decompress(compressed)
                    band_depth = np.frombuffer(decompressed, dtype=np.uint8, count=(rows*cols)+offs)
                    band_depth = band_depth[offs:]  # one line offset by convention?
                    band_depth = band_depth.reshape((rows, cols))

                    # convert data type
                    band_depth = band_depth.astype(dtype=np.float32) / 255.0 * data_type_scaling

                    # normalize to the depth of the library spectrum, translating to aerial fractions
                    library_normalized_band_depth = band_depth / library_band_depth

                    # convert values < 0, > 1, or bad (nan/inf) to 0
                    library_normalized_band_depth[np.logical_not(
                        np.isfinite(library_normalized_band_depth))] = 0
                    library_normalized_band_depth[library_normalized_band_depth < 0] = 0
                    library_normalized_band_depth[library_normalized_band_depth > 1] = 1

                    # determine the mix of EMIT minerals
                    current_mixture_fractions = emit_10_mixture_fractions[tetra_record_numbers.index(
                        record)].copy()
                    current_mixture_fractions = current_mixture_fractions.reshape(
                        1, 1, num_minerals)
                    out_data = out_data + \
                        library_normalized_band_depth.reshape(
                            (rows, cols, 1)) @ current_mixture_fractions

                    # Calculate uncertainty
                    if args.calculate_uncertainty:
                        mixture_uncertainty = calculate_uncertainty(wavelengths, observed_reflectance,
                                                                    observed_reflectance_uncertainty,
                                                                    library_reflectance[library_records.index(record), :], features[0])
                        out_uncertainty = out_uncertainty + \
                                          mixture_uncertainty.reshape((rows, cols, 1)) @ current_mixture_fractions



                except FileNotFoundError:
                    logging.error(filename+'.depth.gz.hdr not found')

        # SMALL keyword tells us to find the library record number
        if 'SMALL' in expert_file_text[expert_line_index]:
            record = int(expert_file_text[expert_line_index].strip().split()[3])

        # 'define output' keyword tells us to get the 8 DN 255 scaling factor
        if 'define output' in expert_file_text[expert_line_index]:
            filename = expert_file_text[expert_line_index+2].strip().split()[0]
            data_type_scaling = float(expert_file_text[expert_line_index+3].strip().split()[4])

        # 'define features' means we've found the location to get the critical feature elements:
        #  the requisite wavelengths for now.  currently continuum removal threshold ct and lct/rct ignored
        if 'define features' in expert_file_text[expert_line_index]:
            feature_line_index = expert_line_index + 1
            features = []
            while ('endfeatures' not in expert_file_text[feature_line_index]):
                toks = expert_file_text[feature_line_index].strip().split()
                if len(toks) > 5 and toks[0].startswith('f') and toks[1] == 'DLw':
                    features.append([float(f) for f in toks[2:6]])
                feature_line_index = feature_line_index + 1
        expert_line_index = expert_line_index + 1

    # write as BIL interleave
    out_data = np.transpose(out_data, (0, 2, 1))
    with open(args.output, 'wb') as fout:
        fout.write(out_data.astype(dtype=np.float32).tobytes())
    envi.write_envi_header(args.output+'.hdr', out_hdr)

    out_uncertainty = np.transpose(out_uncertainty, (0, 2, 1))
    with open(args.output + '_uncert', 'wb') as fout:
        fout.write(out_uncertainty.astype(dtype=np.float32).tobytes())
    envi.write_envi_header(args.output+'_uncert.hdr', out_hdr)
    emit_utils.common_logs.logtime()

if __name__ == "__main__":
    main()
