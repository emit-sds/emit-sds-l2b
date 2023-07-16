# Philip G. Brodrick


import argparse
import gzip
import numpy as np
import spectral.io.envi as envi
import emit_utils.file_checks
from emit_utils.file_checks import envi_header
import logging
import emit_utils.common_logs
import pandas as pd



def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('output_base', type=str, metavar='OUTPUT')
    parser.add_argument('band_depth_file', type=str, metavar='Band Depth file.  4 bands (G1 BD, G1 Ref, G2 BD, G2 Ref)')
    parser.add_argument('--band_depth_unc_file', type=str, default=None)
    parser.add_argument('--mineral_groupings_matrix', type=str, default='data/mineral_grouping_matrix_20230503.csv')
    parser.add_argument('--style', type=str, default='EMIT10', choices=['EMIT10', 'Categories', 'All'])
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--log_level', type=str, default='INFO')
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    emit_utils.common_logs.logtime()
    band_depth_ds = envi.open(envi_header(args.band_depth_file))
    band_depth = band_depth_ds.open_memmap(interleave='bip').copy()
    if args.band_depth_unc_file is not None:
        band_depth_unc_ds = envi.open(envi_header(args.band_depth_unc_file))
        band_depth_unc = band_depth_unc_ds.open_memmap(interleave='bip').copy()

    logging.info('Loading complete, set up output file(s)')
    output_header = band_depth_ds.metadata.copy()

    output_header['interleave'] = 'bil'
    output_header['data type'] = 10
    output_header['header offset'] = 0
    output_header['data type'] = 4

    cols = int(output_header['samples'])
    rows = int(output_header['lines'])
    
    out_min_file = f'{args.output_base}_mineral'
    out_minunc_file = f'{args.output_base}_mineraluncert'
    
    out_cat_file = f'{args.output_base}_category'
    out_catunc_file = f'{args.output_base}_categoryuncert'
    
    
    # Read in the mineral groupings
    mineral_groupings = pd.read_csv(args.mineral_groupings_matrix)
    constituent_groups = mineral_groupings['Group']
    if args.style == 'EMIT10' or args.style == 'All':
        mineral_names = [x for _x, x in enumerate(list(mineral_groupings)) if _x > 6 and _x < 17]
        mineral_abundance = np.array(mineral_groupings[mineral_names])
        mineral_abundance[np.isnan(mineral_abundance)] = 0
        mineral_abundance[mineral_abundance == -1] = 1
        num_minerals = len(mineral_names)
        output_header['bands'] = num_minerals
        output_header['band names'] = mineral_names
        envi.write_envi_header(envi_header(out_min_file), output_header)
        if args.band_depth_unc_file is not None:
            envi.write_envi_header(envi_header(out_minunc_file), output_header)
        out_min_data = np.zeros((rows, cols, num_minerals), dtype=np.float32)

    if args.style == 'Categories' or args.style == 'All':
        category_names = [x for _x, x in enumerate(list(mineral_groupings)) if _x >= 17]
        category_abundance = np.array(mineral_groupings[category_names])
        category_abundance[np.isnan(category_abundance)] = 0
        category_abundance[category_abundance == -1] = 1
        category_names.append('Other')

        other_abun = np.zeros((category_abundance.shape[0],1))
        other_abun[np.sum(mineral_abundance,axis=1) + np.sum(category_abundance,axis=1) == 0,-1] = 1
        category_abundance = np.hstack([category_abundance, other_abun])
        num_categories = len(category_names)
    
   
        output_header['bands'] = num_categories
        output_header['band names'] = category_names
        envi.write_envi_header(envi_header(out_cat_file), output_header)
        if args.band_depth_unc_file is not None:
            envi.write_envi_header(envi_header(out_catunc_file), output_header)
        out_cat_data = np.zeros((rows, cols, num_categories), dtype=np.float32)
    
    
    if args.band_depth_unc_file is not None:
        out_min_unc = np.zeros((rows, cols, num_minerals), dtype=np.float32)
        out_cat_unc = np.zeros((rows, cols, num_categories), dtype=np.float32)

    print(band_depth.shape)
    for _c in range(mineral_abundance.shape[0]):
        print(_c)
        ## Now pull out the rest of the binary file and reshape
        #band_depth = np.frombuffer(decompressed, dtype=np.uint8, count=(rows * cols), offset=header_size)
        #band_depth = band_depth.reshape((rows, cols))

        if constituent_groups[_c] == 1:
            bd_band = 0
            ref_band = 1
        if constituent_groups[_c] == 2:
            bd_band = 2
            ref_band = 3
        ref_match = band_depth[:,:,ref_band:ref_band+1] == (_c+1)
        if np.sum(ref_match) == 0:
            continue

        # determine the mix of EMIT minerals
        if args.style == 'EMIT10' or args.style == 'All':
            current_min_fractions = mineral_abundance[_c, :].reshape(1, 1, -1)
            current_min_fractions[current_min_fractions != 0] = 1
            unmixed_outputs = (ref_match * band_depth[:,:,bd_band:bd_band+1]) @ current_min_fractions
            out_min_data = out_min_data + unmixed_outputs
            if args.band_depth_unc_file is not None:
                unmixed_outputs = (ref_match * band_depth_unc[:,:,bd_band:bd_band+1]) @ current_min_fractions
                out_min_unc = out_min_unc + unmixed_outputs
        
        if args.style == 'Categories' or args.style == 'All':
            current_cat_fractions = category_abundance[_c, :].reshape(1, 1, -1)
            current_cat_fractions[current_cat_fractions != 0] = 1
            unmixed_outputs = (ref_match * band_depth[:,:,bd_band:bd_band+1]) @ current_cat_fractions
            out_cat_data = out_cat_data + unmixed_outputs
            if args.band_depth_unc_file is not None:
                unmixed_outputs = (ref_match * band_depth_unc[:,:,bd_band:bd_band+1]) @ current_cat_fractions
                out_cat_unc = out_cat_unc + unmixed_outputs


    if args.style == 'EMIT10' or args.style == 'All':
        logging.info('Writing output')
        # write as BIL interleave
        out_min_data = np.transpose(out_min_data, (0, 2, 1))
        with open(out_min_file, 'wb') as fout:
            fout.write(out_min_data.astype(dtype=np.float32).tobytes())

        if args.band_depth_unc_file is not None:
            # write as BIL interleave
            out_min_unc = np.transpose(out_min_unc, (0, 2, 1))
            with open(out_minunc_file, 'wb') as fout:
                fout.write(out_min_unc.astype(dtype=np.float32).tobytes())

    if args.style == 'Categories' or args.style == 'All':
        
        logging.info('Writing output')
        # write as BIL interleave
        out_cat_data = np.transpose(out_cat_data, (0, 2, 1))
        with open(out_cat_file, 'wb') as fout:
            fout.write(out_cat_data.astype(dtype=np.float32).tobytes())
    
        if args.band_depth_unc_file is not None:
            # write as BIL interleave
            out_cat_unc = np.transpose(out_cat_unc, (0, 2, 1))
            with open(out_catunc_file, 'wb') as fout:
                fout.write(out_cat_unc.astype(dtype=np.float32).tobytes())

    emit_utils.common_logs.logtime()


if __name__ == "__main__":
    main()
