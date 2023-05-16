# David R. Thompson and Philip G. Brodrick


import argparse
import gzip
import numpy as np
from scipy.interpolate import interp1d
import spectral.io.envi as envi
import emit_utils.file_checks
from emit_utils.file_checks import envi_header
import logging
import emit_utils.common_logs
import os
import tetracorder
from collections import OrderedDict
import json
import matplotlib.pyplot as plt
import pandas as pd

# TODO: Get these from....direct input?  Configuration file?
MINERAL_FRACTION_FILES = [\
    'calcite.group2.txt',
    'chlorite.group2.txt',
    'dolomite.group2.txt',
    'goethite-all-for-reference.group1.txt',
    'gypsum.group2.txt',
    'hematite-all-for-reference.group1.txt',
    'illite+muscovite.group2.txt',
    'kaolinite.group2.txt',
    'montmorillonite.group2.txt',
    'vermiculite.group2.txt',
    ]


def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('tetracorder_output_base', type=str, metavar='TETRA_OUTPUT_DIR')
    parser.add_argument('mineral_groupings_matrix', type=str)
    parser.add_argument('output_base', type=str, metavar='OUTPUT')
    parser.add_argument('--expert_system_file', type=str, default='cmd.lib.setup.t5.27c1', metavar='EXPERT_SYS_FILE')
    parser.add_argument('--reflectance_file', type=str, metavar='REFLECTANCE_FILE')
    parser.add_argument('--reflectance_uncertainty_file', type=str, metavar='REFLECTANCE_UNCERTAINTY_FILE')
    parser.add_argument('--calculate_uncertainty', action='store_true')
    parser.add_argument('--reference_library', '-slib', type=str, default='Spectral-Library-Reader-master/s06av18a_envi', metavar='REFERENCE_LIBRARY')
    parser.add_argument('--research_library', '-rlib', type=str, default='Spectral-Library-Reader-master/r06av18a_envi', metavar='RESEARCH_LIBRARY')
    parser.add_argument('--diagnostic_files', action='store_true')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--log_level', type=str, default='INFO')
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    emit_utils.common_logs.logtime()

    if args.calculate_uncertainty:
        args.calculate_uncertainty = True
        emit_utils.file_checks.check_raster_files([args.reflectance_file, args.reflectance_uncertainty_file], map_space=False)

        refl_dataset = envi.open(envi_header(args.reflectance_file))
        # bulk-copy in the observed reflectance dataset.  This pre-supposes sufficient memory, but reduces IO costs on non-SSDs.
        observed_reflectance = refl_dataset.open_memmap(interleave='bip', writable=False).copy()
        wavelengths = np.array([float(x) for x in refl_dataset.metadata['wavelength']])

        # bulk-copy in the observed reflectance uncertainty dataset.  This pre-supposes sufficient memory, but reduces IO costs on non-SSDs.
        observed_reflectance_uncertainty_dataset = envi.open(envi_header(args.reflectance_uncertainty_file))
        observed_reflectance_uncertainty = observed_reflectance_uncertainty_dataset.open_memmap(interleave='bip', writable=False).copy()
    else:
        args.calculate_uncertainty = False

    expert_system_file = os.path.join(args.tetracorder_output_base, args.expert_system_file)
    if os.path.isfile(expert_system_file) is False:
        logging.error(f'No expert system file found, expected at: {expert_system_file}. Look for candidates in'
                      f'{args.tetracorder_output_base} that start with cmd.lib.setup')
        raise AttributeError('Could not find expert system file, see log for details.')

    decoded_expert = tetracorder.decode_expert_system(expert_system_file, log_file=args.log_file,
                                                      log_level=args.log_level)

    mineral_groupings = pd.read_csv(args.mineral_groupings_matrix)
    mineral_names = [x for _x, x in enumerate(list(mineral_groupings)) if _x >= list(mineral_groupings).index('Calcite') and _x <= list(mineral_groupings).index('Vermiculite')]
    category_names = [x for _x, x in enumerate(list(mineral_groupings)) if _x > list(mineral_groupings).index('Vermiculite')]

    library_names = mineral_groupings['Library']
    library_records = mineral_groupings['Record']
    el = envi.open(envi_header(args.reference_library))
    libraries = {'splib06': {'library_reflectance': el.spectra.copy(), 'library_records': [int(q) for q in el.metadata['record']]}}
    el = envi.open(envi_header(args.research_library))
    libraries['sprlb06'] = {'library_reflectance': el.spectra.copy(), 'library_records': [int(q) for q in el.metadata['record']]}
    del el

    mineral_abundance = np.array(mineral_groupings[mineral_names])
    category_abundance = np.array(mineral_groupings[category_names])
    mineral_abundance[np.isnan(mineral_abundance)] = 0
    category_abundance[np.isnan(category_abundance)] = 0
    mineral_abundance[mineral_abundance == -1] = 1
    category_abundance[category_abundance == -1] = 1
    group = np.array(mineral_groupings['Group'])
    
    category_names.append('Other')
    other_abun = np.zeros((category_abundance.shape[0],1))
    other_abun[np.sum(mineral_abundance,axis=1) + np.sum(category_abundance,axis=1) == 0,-1] = 1
    category_abundance = np.hstack([category_abundance, other_abun])
    
    filenames = np.array(mineral_groupings['Filename']).tolist()
    depth_filenames = [x + '.depth.gz' for x in filenames]
    fit_filenames = [x + '.fit.gz' for x in filenames]
    
    # Should be brought in from tetracorder, but not conveniently unpacked from expert system, and currrently univerally 0.5
    scaling = [0.5 for x in range(len(depth_filenames))]

    logging.info('Loading complete, set up output file(s)')
    # Set up output files
    input_header = None
    for uf in depth_filenames:
        header_name = os.path.join(args.tetracorder_output_base, envi_header(uf))
        print(header_name)
        if os.path.isfile(header_name):
            input_header = envi.read_envi_header(header_name)
            logging.info(f'using {header_name} as header base')
            break
        else:
            logging.debug(f'{header_name} not found - skipping in base header search')
    output_header = input_header.copy()
    if 'file compression' in output_header.keys():
        del output_header['file compression']
    if 'wavelengths' in output_header.keys():
        del output_header['wavelengths']

    output_header['interleave'] = 'bil'
    output_header['data type'] = 4
    output_header['bands'] = 4
    output_header['header offset'] = 0
    output_header['band names'] = ['Group 1 Band Depth', 'Group 1 Index', 'Group 2 Band Depth', 'Group 2 Index']
    cols = int(input_header['samples'])
    rows = int(input_header['lines'])
    
    out_complete_file = f'{args.output_base}_complete'
    envi.write_envi_header(envi_header(out_complete_file), output_header)

    output_header['band names'] = ['Group 1 Band Depth Uncertainty', 'Group 1 Fit', 'Group 2 Band Depth Uncertainty', 'Group 2 Fit']
    unc_complete_file = f'{args.output_base}_complete_uncertainty'
    envi.write_envi_header(envi_header(unc_complete_file), output_header)

    out_complete = np.zeros((rows, cols, 4), dtype=np.float32)
    unc_complete = np.zeros((rows, cols, 4), dtype=np.float32)
    
    for _c, constituent_file in enumerate(depth_filenames):
        if _c %10 == 0:
            print(f'{_c}/{len(depth_filenames)}')

        band_depth_file = os.path.join(args.tetracorder_output_base, constituent_file)
        fit_file = os.path.join(args.tetracorder_output_base, fit_filenames[_c])
        if os.path.isfile(band_depth_file) is False:
            logging.info(f'{band_depth_file} not found: skipping')
            continue

        # read band depth
        band_depth = read_tetra_file(band_depth_file, rows, cols, scaling[_c])
        fit = read_tetra_file(fit_file, rows, cols, scaling[_c])

        valid_pixels = band_depth > 0


        if library_records[_c] in libraries[library_names[_c]]['library_records']:
            local_lib_refl = (libraries[library_names[_c]]['library_reflectance'])[libraries[library_names[_c]]['library_records'].index(library_records[_c]),:]
            local_features = decoded_expert[filepath_to_key(constituent_file)]['features']  
            if args.calculate_uncertainty and np.sum(valid_pixels) > 0:
                loc_unc = calculate_uncertainty(wavelengths, observed_reflectance[valid_pixels,:], observed_reflectance_uncertainty[valid_pixels,:], local_lib_refl, local_features)

            if group[_c] == 1:
                out_complete[valid_pixels,0] = band_depth[band_depth > 0]
                out_complete[valid_pixels,1] = _c + 1
                if args.calculate_uncertainty and np.sum(valid_pixels) > 0:
                    unc_complete[valid_pixels,0] = loc_unc
                unc_complete[valid_pixels,1] = fit[band_depth > 0]

            if group[_c] == 2:
                out_complete[valid_pixels,2] = band_depth[band_depth > 0]
                out_complete[valid_pixels,3] = _c + 1
                if args.calculate_uncertainty and np.sum(valid_pixels) > 0:
                    unc_complete[valid_pixels,2] = loc_unc
                unc_complete[valid_pixels,3] = fit[band_depth > 0]
        else:
            logging.warning(f'Library record {library_records[_c]} for {constituent_file} not found in {library_names[_c]}')
        

    logging.info('Writing output')
    # write as BIL interleave
    out_complete = np.transpose(out_complete, (0, 2, 1))
    with open(out_complete_file, 'wb') as fout:
        fout.write(out_complete.astype(dtype=np.float32).tobytes())

    unc_complete = np.transpose(unc_complete, (0, 2, 1))
    with open(unc_complete_file, 'wb') as fout:
        fout.write(unc_complete.astype(dtype=np.float32).tobytes())
        
    emit_utils.common_logs.logtime()


def read_tetra_file(filename: str, rows: int, cols: int, scaling: float) -> np.ndarray:
    """
    Reads a compressed tetracorder file and returns a numpy array

    Args:
        filename (str): The name of the file to read
        rows (int): The number of rows in the file
        cols (int): The number of columns in the file
        scaling (float): A scaling factor to apply to the data

    Raises:
        AttributeError: Incorrect file format, no LBLSIZE found in VICAR header

    Returns:
        np.ndarray: A numpy array containing the band depth data
    """
    # read band depth
    with open(filename, 'rb') as fin:
        compressed = fin.read()
    decompressed = gzip.decompress(compressed)

    band_depth_header = envi.read_envi_header(filename + '.hdr')
    offs = max(50,int(band_depth_header['header offset']))
    vicar = decompressed[:offs].decode('ascii').split(' ')[0]
    if vicar[:7] != 'LBLSIZE':
        raise AttributeError(f'Incorrect file format {filename},'
                             'no LBLSIZE found in VICAR header')
    # Read the header size from the VICAR header
    header_size = int(vicar.split('=')[-1])

    # Now pull out the rest of the binary file and reshape
    read_data = np.frombuffer(decompressed, dtype=np.uint8, count=(rows * cols), offset=header_size).copy()
    read_data = read_data.reshape((rows, cols))

    # convert data type
    read_data = read_data.astype(dtype=np.float32) / 255.0 * scaling

    return read_data


def calculate_band_depth(wavelengths: np.array, reflectance: np.array, feature: tuple, record: int = None, name: str = None):
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

    if record is not None and name is not None:
        fig = plt.figure()
        plt.plot(wavelengths, reflectance)
        plt.plot(wavelengths[feature_inds], reflectance[feature_inds])
        plt.plot(wavelengths[feature_inds], continuum[feature_inds])
        bd_ind = np.argmax(depths)
        plt.plot([wavelengths[feature_inds][bd_ind], wavelengths[feature_inds][bd_ind]], [reflectance[feature_inds][bd_ind], continuum[feature_inds][bd_ind]])
        plt.title(f'{record}   |||   {name}\nBD = {max(depths)}')
        plt.savefig(f'/beegfs/scratch/brodrick/emit/reflectance_analyses/figs/{record}_{name.replace("/","-")}.png',dpi=200,bbox_inches='tight')

    return max(depths)


def cont_rem(wavelengths: np.array, reflectance: np.array, feature: tuple):
    """This function removes the continuum from a reflectance spectrum
    using the continuum removal method of Kaufman and Tanre (1992)

    Args:
        wavelengths (array): An array of wavelengths
        reflectance (array): An array of reflectance values
        feature (float): The feature value to use for removing the continuum

    Returns:
        array: An array of continuum removed reflectance values
    """
    left_inds = np.where(np.logical_and(wavelengths >= feature[0], wavelengths <= feature[1]))[0]
    left_x = wavelengths[int(left_inds.mean())]
    left_y = reflectance[:,left_inds].mean()

    right_inds = np.where(np.logical_and(wavelengths >= feature[2], wavelengths <= feature[3]))[0]
    right_x = wavelengths[int(right_inds.mean())]
    right_y = reflectance[:,right_inds].mean()
    
    feature_inds = np.where(np.logical_and(wavelengths >= feature[0], wavelengths <= feature[3]))[0]

    continuum_fun = interp1d([left_x, right_x], [left_y, right_y], bounds_error=False, fill_value='extrapolate')
    continuum = continuum_fun( np.ones((reflectance.shape[0], len(feature_inds))) * wavelengths[feature_inds][np.newaxis,:])
    depths = 1 - reflectance[:,feature_inds] / continuum, feature_inds
    return depths


def calculate_uncertainty(wavelengths: np.array, observed_reflectance: np.array, observed_reflectance_uncertainty: np.array, library_reflectance: np.array, feature_set: dict):
    """ Calculate the uncertainty of the Clark, 2003 continuum normalized band depth of a particular feature.

        Assumed that band depth = a * L(w_star)
        where a = argmin [ sum_w (a L(w) + b - D(w))^2 ]
                = n sum_w (L(w) D(w)) - sum_w (L(w)) sum_w (D(w)) / n sum_w (L(w)^2) - (sum_w (L(w)))^2

        Assuming that the uncertainty of the continuum removed reflectances equals the uncertainty of the reflectance at that wavelength (for simplicity),
        and given that reflectance uncertainties are assumed to be independent, and that the uncertainty of sum_i (a_i A_i) = sqrt(sum_i (a_i^2 sigma_Ai^2 )):

        bd_unc = L(w_star) / (n sum_w (L(w)^2) - (sum_w (L(w)))^2) * sqrt( n^2 sum_w (L(w)^2 sigma_r(w)^2) + sum_w ( sum_w ( L(w) )^2 * sigma_r(w) ^ 2) ) 

        band depth uncertainties are summed for each feature, proportionately to the number of wavelengths included in that feature

    Args:
        wavelengths: an array of wavelengths corresponding to given reflectance values
        observed_reflectance: an array of observed reflectance values
        observed_reflectance_uncertainty: an array of uncertainties of the observed reflectance values
        library_reflectance: an array of library reference reflectance values
        feature: definition of the feature to calculate band depth for.
    :Returns
        uncertainty: the uncertainty for the band depth as defined in Clark, 2003.
    """

    if np.any(wavelengths > 100):
        wl_micron = wavelengths / 1000
    else:
        wl_micron = wavelengths

    unc_output = []
    for feature in feature_set:

        if feature['feature_type'] in ['MLw', 'DLw']:

            # check to make sure that there are at least some valid reflectances insdide the given continuum
            if np.any(observed_reflectance[0, np.where(np.logical_and(wl_micron >= feature['continuum'][0], wl_micron <= feature['continuum'][1]))] > 0) and \
               np.any(observed_reflectance[0, np.where(np.logical_and(wl_micron >= feature['continuum'][2], wl_micron <= feature['continuum'][3]))] > 0):


                # Continuum removal
                lib_cont, wl_inds = cont_rem(wl_micron, library_reflectance.reshape((1,-1)), feature['continuum'])
                obs_cont, wl_inds = cont_rem(wl_micron, observed_reflectance, feature['continuum'])

                # remove extra library dimension
                lib_cont = np.squeeze(lib_cont)

                # subset uncertainty for convenience
                obs_rfl_unc = observed_reflectance_uncertainty[:, wl_inds]

                # some frequenty summation terms
                lib_cont_sum = np.sum(lib_cont)
                lib_cont_squared_sum = np.sum(lib_cont**2)

                # find the band-depth index from the library
                w_star = np.argmax(lib_cont) 
        
                # calculate a 
                const =  lib_cont[w_star] / (len(lib_cont) * lib_cont_squared_sum - lib_cont_sum**2)

                inner_term = np.sum(np.power(obs_rfl_unc,2) * np.power(lib_cont * const * len(lib_cont),2)[np.newaxis,:],axis=1) +\
                             np.sum(np.power(const * lib_cont_sum,2) * np.power(obs_rfl_unc,2), axis=1)

                loc_unc = np.sqrt(inner_term)

                # for analysis only
                #slope = (np.sum(obs_cont * lib_cont[np.newaxis,:],axis=1)*len(lib_cont) - np.sum(obs_cont, axis=1)*np.sum(lib_cont)) /( np.sum(lib_cont**2) * len(lib_cont) - np.sum(lib_cont)**2 )
                #offset = (np.sum(obs_cont, axis=1) - slope * lib_cont_sum) / len(lib_cont)

                unc_output.append([loc_unc, len(lib_cont)])

    total_bands = np.sum([x[1] for x in unc_output])
    uncertainty = np.sum(np.vstack([x[0] * x[1] for x in unc_output]),axis=0) / total_bands
    return uncertainty


def filepath_to_key(value: str):
    """ Convert a filepath to a dictionary key in a systematic way.
    Args:
        value: filepath to convert

    Returns:
        string key to dictionary

    """
    return value.split('.depth.gz')[0]


if __name__ == "__main__":
    main()

