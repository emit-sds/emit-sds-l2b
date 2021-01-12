# Philip G. Brodrick


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
from aggregator import calculate_band_depth

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
    parser.add_argument('-log_file', type=str, default=None)
    parser.add_argument('-log_level', type=str, default='INFO')
    parser.add_argument('-condensed', type=int, choices=[0,1], default=1)
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    emit_utils.common_logs.logtime()

    args.condensed = args.condensed == 1

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

    import glob
    potential_files = glob.glob(f'{args.tetracorder_output_base}/group.1um/*depth.gz') + glob.glob(f'{args.tetracorder_output_base}/group.2um/*depth.gz')
    #potential_files = glob.glob(f'{args.tetracorder_output_base}/group*/*depth.gz')
    actual_files = []
    for fi in potential_files:
        if fi.replace(args.tetracorder_output_base,'') not in unique_file_names:
            if np.any([x in fi for x in ['vegetation','veg','snow','ice','algae','chlorophyll','plastic','organic']]):
                print(fi)
                continue
            elif filepath_to_key(fi) not in decoded_expert.keys():
                print(fi)
                continue
            else:
                actual_files.append(fi.replace(args.tetracorder_output_base,''))
                #print(fi)

    records = []
    for fi in actual_files:
        records.append(decoded_expert[filepath_to_key(fi)]['record'])

    dn_scales = []
    for fi in actual_files:
        dn_scales.append(decoded_expert[filepath_to_key(fi)]['data_type_scaling'])

    library_names = []
    for fi in actual_files:
        library_names.append(decoded_expert[filepath_to_key(fi)]['spectral_library'])


    input_header = envi.read_envi_header(args.tetracorder_output_base + '/' + actual_files[0] + '.hdr')
    output_header = input_header.copy()
    if 'file compression' in output_header.keys():
        del output_header['file compression']
    if 'wavelengths' in output_header.keys():
        del output_header['wavelengths']

    output_header['interleave'] = 'bil'
    output_header['data type'] = 4
    if args.condensed:
        output_header['bands'] = 1
    else:
        output_header['bands'] = len(actual_files)
    output_header['header offset'] = 0
    output_header['band names'] =  actual_files
    cols = int(input_header['samples'])
    rows = int(input_header['lines'])
    envi.write_envi_header(f'{args.output_base}.hdr', output_header)

    if args.condensed:
        out_data = np.zeros((rows, cols, 1), dtype=np.float32)
    else:
        out_data = np.zeros((rows, cols, len(actual_files)), dtype=np.float32)


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

        band_depths = np.zeros(len(actual_files))
        for _f, (filename, library_name, record) in enumerate(zip(actual_files, library_names, records)):
            if library_name == key:
                band_depths[_f] = calculate_band_depth(wavelengths,
                                                       library_reflectance[library_records.index(record), :],
                                                       decoded_expert[filepath_to_key(filename)]['features'][0]['continuum'])
        libraries[key]['band_depths'] = band_depths


    for _c, constituent_file in enumerate(actual_files):

        ref_lib = libraries[library_names[_c]]
        fullpath_constituent_file = args.tetracorder_output_base + '/' + constituent_file
        logging.debug(fullpath_constituent_file)

        # read band depth
        with open(fullpath_constituent_file, 'rb') as fin:
            compressed = fin.read()
        decompressed = gzip.decompress(compressed)

        band_depth_header = envi.read_envi_header(fullpath_constituent_file + '.hdr')
        offs = int(band_depth_header['header offset'])
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
        band_depth = band_depth.astype(dtype=np.float32) / 255.0 * dn_scales[_c]
        library_normalized_band_depth = band_depth / ref_lib['band_depths'][_c]

        if args.condensed:
            out_data[...,0] += library_normalized_band_depth
        else:
            out_data[...,_c] = library_normalized_band_depth


    out_data = np.transpose(out_data, (0, 2, 1))
    with open(args.output_base, 'wb') as fout:
        fout.write(out_data.astype(dtype=np.float32).tobytes())


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
