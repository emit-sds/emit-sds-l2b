# David R. Thompson and Philip G. Brodrick


import argparse
import gzip
import numpy as np
from scipy.interpolate import interp1d
import spectral.io.envi as envi
import emit_utils.file_checks
import logging
import emit_utils.common_logs
import os
import tetracorder
from collections import OrderedDict

# TODO: Get these from....direct input?  Configuration file?
MINERAL_FRACTION_FILES = [\
    'kaolinite.group2.txt',
    'calcite.group2.txt',
    'dolomite.group2.txt',
    'vermiculite.group2.txt',
    'hematite-all-for-reference.group1.txt',
    'goethite-all-for-reference.group1.txt',
    'gypsum.group2.txt',
    'chlorite.group2.txt',
    'illite.group2.txt',
    'montmorillonite.group2.txt',
    ]

SPECTRAL_REFERENCE_LIBRARY = {\
    'splib06': 'Spectral-Library-Reader-master/s06av18a_envi',
    'sprlb06': 'Spectral-Library-Reader-master/r06av18a_envi',
    }


def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('tetracorder_output_base', type=str, metavar='TETRA_OUTPUT_DIR')
    parser.add_argument('output_base', type=str, metavar='OUTPUT')
    parser.add_argument('-spectral_reference_library_config', type=str, metavar='TETRA_LIBRARY_CONFIG_FILE')
    parser.add_argument('-expert_system_file', type=str, default='cmd.lib.setup.t5.2d4', metavar='EXPERT_SYS_FILE')
    parser.add_argument('-calculate_uncertainty', type=int, choices=[0,1], metavar='CALCULATE_UNCERTAINTY')
    parser.add_argument('-reflectance_file', type=str, metavar='REFLECTANCE_FILE')
    parser.add_argument('-reflectance_uncertainty_file', type=str, metavar='REFLECTANCE_UNCERTAINTY_FILE')
    parser.add_argument('-detailed_outputs', type=int, choices=[0, 1], default=0, metavar='DETAILED_OUTPUTS')
    parser.add_argument('-log_file', type=str, default=None)
    parser.add_argument('-log_level', type=str, default='INFO')
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    emit_utils.common_logs.logtime()

    args.detailed_outputs = args.detailed_outputs == 1

    if args.calculate_uncertainty == 1:
        args.calculate_uncertainty = True
        emit_utils.file_checks.check_raster_files([args.reflectance_file, args.reflectance_uncertainty_file], map_space=False)

        refl_dataset = envi.open(args.reflectance_file + '.hdr')
        observed_reflectance = refl_dataset.open_memmap(interleave='bil', writable=False)

        observed_reflectance_uncertainty_dataset = envi.open(args.reflectance_uncertainty_file + '.hdr')
        observed_reflectance_uncertainty = observed_reflectance_uncertainty_dataset.open_memmap(interleave='bil',
                                                                                                writable=False)
    else:
        args.calculate_uncertainty = False

    expert_system_file = os.path.join(args.tetracorder_output_base, args.expert_system_file)
    if os.path.isfile(expert_system_file) is False:
        logging.error(f'No expert system file found, expected at: {expert_system_file}. Look for candidates in'
                      f'{args.tetracorder_output_base} that start with cmd.lib.setup')
        raise AttributeError('Could not find expert system file, see log for details.')

    decoded_expert = tetracorder.decode_expert_system(expert_system_file, log_file=args.log_file,
                                                      log_level=args.log_level)

    mff = [os.path.join(args.tetracorder_output_base, 'cmds.abundances', 'lists.of.files.by.mineral', x) for x
           in MINERAL_FRACTION_FILES]
    mineral_fractions = tetracorder.read_mineral_fractions(mff)
    num_minerals = len(mineral_fractions.keys())

    logging.info('Organizing files to aggregate')
    unique_file_names, fractions, scaling, library_names, records = unique_file_fractions(mineral_fractions, decoded_expert)

    logging.info('Loading complete, set up output file(s)')
    # Set up output files
    input_header = envi.read_envi_header(os.path.join(args.tetracorder_output_base, unique_file_names[0] + '.hdr'))
    output_header = input_header.copy()
    if 'file compression' in output_header.keys():
        del output_header['file compression']
    if 'wavelengths' in output_header.keys():
        del output_header['wavelengths']

    output_header['interleave'] = 'bil'
    output_header['data type'] = 4
    output_header['bands'] = num_minerals
    output_header['header offset'] = 0
    output_header['band names'] = mineral_fractions.keys()
    cols = int(input_header['samples'])
    rows = int(input_header['lines'])
    envi.write_envi_header(f'{args.output_base}.hdr', output_header)
    envi.write_envi_header(f'{args.output_base}_uncert.hdr', output_header)

    out_data = np.zeros((rows, cols, num_minerals), dtype=np.float32)
    if args.calculate_uncertainty:
        out_uncertainty = np.zeros((rows, cols, num_minerals), dtype=np.float32)

    #TODO: include option to read in from config file
    spectral_reference_library_files = SPECTRAL_REFERENCE_LIBRARY
    libraries = {}
    for key, item in spectral_reference_library_files.items():
        library = envi.open(item + '.hdr', item)
        library_reflectance = library.spectra.copy()
        library_records = [int(q) for q in library.metadata['record']]
        wavelengths = np.array([float(q) for q in library.metadata['wavelength']])

        if ';;;' in key:
            key = key.replace(';;;', ',')
            logging.debug(f'found comma replacement, now: {key}')

        libraries[key] = {'reflectance': library_reflectance,
                          'library_records': library_records, 'wavelengths': wavelengths}

        band_depths = np.zeros(fractions.shape[0])
        for _f, (filename, library_name, record) in enumerate(zip(unique_file_names, library_names.tolist(), records.tolist())):
            if library_name == key:
                band_depths[_f] = calculate_band_depth(wavelengths,
                                                       library_reflectance[library_records.index(record), :],
                                                       decoded_expert[filepath_to_key(filename)]['features'][0]['continuum'])
        libraries[key]['band_depths'] = band_depths

    logging.info('Begin actual aggregation')
    if args.detailed_outputs:
        detailed_outputs = {}
        for name in mineral_fractions.keys():
            detailed_outputs[name] = {}
            detailed_outputs[name]['files'] = []
            detailed_outputs[name]['values'] = []

    for _c, constituent_file in enumerate(unique_file_names):

        ref_lib = libraries[library_names[_c]]

        fullpath_constituent_file = os.path.join(args.tetracorder_output_base, constituent_file)

        # read band depth
        with open(fullpath_constituent_file, 'rb') as fin:
            compressed = fin.read()
        decompressed = gzip.decompress(compressed)

        band_depth_header = envi.read_envi_header(fullpath_constituent_file + '.hdr')
        offs = max(50,int(band_depth_header['header offset']))
        vicar = decompressed[:offs].decode('ascii').split(' ')[0]
        if vicar[:7] != 'LBLSIZE':
            raise AttributeError(f'Incorrect file format {fullpath_constituent_file},'
                                 'no LBLSIZE found in VICAR header')
        # Read the header size from the VICAR header
        header_size = int(vicar.split('=')[-1])

        # Now pull out the rest of the binary file and reshape
        band_depth = np.frombuffer(decompressed, dtype=np.uint8, count=(rows * cols), offset=header_size)
        band_depth = band_depth.reshape((rows, cols))

        # convert data type
        band_depth = band_depth.astype(dtype=np.float32) / 255.0 * scaling[_c]

        # normalize to the depth of the library spectrum, translating to aerial fractions
        library_normalized_band_depth = band_depth / ref_lib['band_depths'][_c]

        # convert values < 0, > 1, or bad (nan/inf) to 0
        library_normalized_band_depth[np.logical_not(
            np.isfinite(library_normalized_band_depth))] = 0
        library_normalized_band_depth[library_normalized_band_depth < 0] = 0
        library_normalized_band_depth[library_normalized_band_depth > 1] = 1

        # determine the mix of EMIT minerals
        current_mixture_fractions = fractions[_c, :].reshape(1, 1, -1)
        unmixed_outputs = library_normalized_band_depth.reshape((rows, cols, 1)) @ current_mixture_fractions
        out_data = out_data + unmixed_outputs

        if args.detailed_outputs:
            for _b in range(unmixed_outputs.shape[-1]):
                if np.sum(unmixed_outputs[...,_b]) > 0:
                    detailed_outputs[list(mineral_fractions.keys())[_b]]['files'].append(constituent_file)
                    detailed_outputs[list(mineral_fractions.keys())[_b]]['values'].append(unmixed_outputs[...,_b])

        # Calculate uncertainty
        if args.calculate_uncertainty and np.sum(library_normalized_band_depth != 0) > 0:
            mixture_uncertainty = calculate_uncertainty(ref_lib['wavelengths'], observed_reflectance,
                                                        observed_reflectance_uncertainty,
                                                        ref_lib['reflectance'][ref_lib['records'].index(records[_c]), :],
                                                        decoded_expert[constituent_file]['features'][0]['continuum'])
            out_uncertainty = out_uncertainty + \
                              mixture_uncertainty.reshape((rows, cols, 1)) @ current_mixture_fractions

    if args.detailed_outputs:
        for mineral in mineral_fractions.keys():
            num_outputs = len(detailed_outputs[mineral]['files'])
            if num_outputs > 0:
                output_header['bands'] = num_outputs
                output_header['band names'] = detailed_outputs[mineral]['files']
                envi.write_envi_header(f'{args.output_base}_{mineral}_details.hdr', output_header)
                towrite = np.stack(detailed_outputs[mineral]['values'],axis=1)
                with open(f'{args.output_base}_{mineral}_details', 'wb') as fout:
                    fout.write(towrite.astype(dtype=np.float32).tobytes())

    detailed_header = output_header.copy()
    # write as BIL interleave
    out_data = np.transpose(out_data, (0, 2, 1))
    with open(args.output_base, 'wb') as fout:
        fout.write(out_data.astype(dtype=np.float32).tobytes())

    if args.calculate_uncertainty:
        out_uncertainty = np.transpose(out_uncertainty, (0, 2, 1))
        with open(f'{args.output_base}_uncert', 'wb') as fout:
            fout.write(out_uncertainty.astype(dtype=np.float32).tobytes())
    emit_utils.common_logs.logtime()


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

    right_term = np.sum(np.power(len(inds)*a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis] - np.sum(a[:,np.newaxis,:]*library_reflectance[np.newaxis,inds,np.newaxis],axis=1)[:,np.newaxis,:],2) * np.power(observed_reflectance_uncertainty[:,inds,:],2),axis=1)

    psi = np.sqrt(left_term*right_term)

    # Safegaurd against NANs that might occur for non-fits
    psi[np.isnan(psi)] = 0

    return psi


def filepath_to_key(value: str):
    """ Convert a filepath to a dictionary key in a systematic way.
    Args:
        value: filepath to convert

    Returns:
        string key to dictionary

    """
    return os.path.basename(value).split('.depth.gz')[0]


def unique_file_fractions(fraction_dict: OrderedDict, decoded_expert: OrderedDict):
    """

    Args:
        fraction_dict: mineral fractions dictionary
        decoded_expert: tetracorder expert system decoded to a dictionary

    Returns:
        unique_file_names: the set of files that must be scanned based on the mineral fractions specified
        fractions: fraction of each of the minerals included in the fraction_dict as a matrix with rows corresponding
            to unique_file_names
        scaling: scaling values for read in band depths, rows correspond to unique_file_names
        library: the reference library used, corresponding to rows of unique_file_names
        record: the reference library record number, corresponding to rows of unique_file_names
    """
    file_names = []
    for key, item in fraction_dict.items():
        file_names = file_names + [x['file'] for x in item]
    unique_file_names = np.unique(file_names).tolist()

    mineral_names = list(fraction_dict.keys())
    fractions = np.zeros((len(unique_file_names), len(fraction_dict)))
    scaling = np.zeros(len(unique_file_names))
    record = np.zeros(len(unique_file_names), dtype=int)
    library = np.empty(len(unique_file_names), dtype="<U10")

    for key, item in fraction_dict.items():
        for constituent in item:
            idx = unique_file_names.index(constituent['file'])
            fractions[idx, mineral_names.index(key)] = constituent['BD_factor']
            scaling[idx] = constituent['DN_scale']
            library[idx] = constituent['spectral_library']
            if 'record' in constituent.keys():
                record[idx] = constituent['record']
            else:
                logging.debug('scanning for record in decoded expert system')
                record[idx] = decoded_expert[filepath_to_key(constituent['file'])]['record']
            logging.debug(f'file: {unique_file_names[idx]}, DN_scale: {scaling[idx]}, library: {library[-1]}, record: {record[idx]}')

    return unique_file_names, fractions, scaling, library, record


if __name__ == "__main__":
    main()
