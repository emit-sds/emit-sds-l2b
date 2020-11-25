# David Thompson


import argparse
import yaml
import json
import numpy as np
import os
import logging
from typing import List
from collections import OrderedDict
import emit_utils.common_logs
import re
import pandas as pd

DEFAULT_GROUPS = [1,2]
DEFAULT_LOGFILE = None
DEFAULT_LOGLEVEL = 'INFO'
DEFAULT_SORT_KEYS = False


def recast_globals(to_recast: List, globals: dict):
    for _i in range(len(to_recast)):
        for key, item in globals.items():
            to_recast[_i] = to_recast[_i].replace(key, item)
    return to_recast


def decode_expert_system(tetra_expert_file, groups=DEFAULT_GROUPS, log_file=DEFAULT_LOGFILE,
                         log_level=DEFAULT_LOGLEVEL):
    """ Convert a tetracorder .txt 'expert system file' into a dictionary.  Not all
    parameters are preserved, only those necessary for EMIT calculation.

    Args:
        tetra_expert_file: tetracorder expert system file to read
        groups: tetracorder groups to use (all others will be ignored)
        log_file: file to write logging to
        log_level: logging level to use

    Returns:
        Coverted dictionary of tetracorder file
    """
    if log_file is None:
        logging.basicConfig(format='%(message)s', level=log_level)
    else:
        logging.basicConfig(format='%(message)s', level=log_level, filename=log_file)

    #emit_utils.common_logs.logtime()

    decoded_expert = {}
    # read expert system file and strip comments
    with open(tetra_expert_file, 'r') as fin:
        expert_file_commented = fin.readlines()

    expert_file_text, orig_lineno = [], []
    for line_index, line in enumerate(expert_file_commented):
        if not line.strip().startswith('\#') or 'TITLE=' in line:
            orig_lineno.append(line_index)
            expert_file_text.append(line)
    del expert_file_commented

    # Go through expert system file one line at a time, after initializing key variables
    expert_line_index, group, spectrum, output_data, header, out_hdr, rows, cols = \
        0, None, None, None, True, None, 0, 0
    constituent_constraints = {}

    globals = {
        '[GLBLFITALL]': '[GLBLFITALL]',
        '[GLBLFDFIT]': '[GLBLFDFIT]',
        '[GLBLDPFIT]': '[GLBLDPFIT]',
        '[GLBLDPFITg2]': '[GLBLDPFITg2]',
        '[G2UMRBDA]': '[G2UMRBDA]',
    }

    while expert_line_index < len(expert_file_text):

        # The Header flag excludes the definitions at the start
        if expert_file_text[expert_line_index].startswith('BEGIN SETUP'):
            header = False
        elif header:
            # Check for globals
            if expert_file_text[expert_line_index].startswith('=='):
                for key in globals.keys():
                    if key in expert_file_text[expert_line_index]:
                        line_remainder = expert_file_text[expert_line_index][len(key)+2].strip()
                        if '\#' in line_remainder:
                            line_remainder = line_remainder[:line_remainder.index('\#')]
                        globals[key] = line_remainder

            expert_line_index = expert_line_index + 1
            continue

        # if we've gotten to the end of the record, time to pull everything together and write our output
        if expert_file_text[expert_line_index].startswith('endaction'):
            if group in groups:
                # Populate the entry
                entry = {}
                constituent_constraints = {}
                entry['longname'] = longname
                entry['group'] = group
                entry['record'] = record
                entry['name'] = name
                entry['data_type_scaling'] = data_type_scaling
                entry['features'] = features
                entry['constituent_constraints'] = constituent_constraints
                decoded_expert[tetra_filename] = entry

        if 'TITLE=' in expert_file_text[expert_line_index]:
            toks = expert_file_text[expert_line_index].strip().split()
            name = toks[1].strip().split('=')[1]
            longname = name
            for _t in range(2, len(toks)):
                longname += ' ' + toks[_t]

        # if keyword 'group' appears, define the current group name
        if expert_file_text[expert_line_index].startswith('group'):
            group = int(expert_file_text[expert_line_index].strip().split()[1])

        # SMALL keyword tells us to find the library record number
        if 'SMALL' in expert_file_text[expert_line_index]:
            record = int(expert_file_text[expert_line_index].strip().split()[3])

        # 'define output' keyword tells us to get the 8 DN 255 scaling factor
        if 'define output' in expert_file_text[expert_line_index]:
            tetra_filename = expert_file_text[expert_line_index+2].strip().split()[0]
            data_type_scaling = float(expert_file_text[expert_line_index+3].strip().split()[4])

        # 'define features' means we've found the location to get the critical feature elements:
        #  the requisite wavelengths for now.  currently continuum removal threshold ct and lct/rct ignored
        valid_feature_constraints = ['ct','lct','rct','lct/rct>','rct/lct>','rcbblc>','rcbblc<','lcbbrc>','lcbbrc<','r*bd>']
        if 'define features' in expert_file_text[expert_line_index]:
            expert_line_index += 1
            features = []
            while ('endfeatures' not in expert_file_text[expert_line_index]):
                toks = expert_file_text[expert_line_index].strip().split()
                if len(toks) > 5 and toks[0].startswith('f') and toks[1] == 'DLw':
                    local_feature = {}
                    local_feature['continuum'] = [float(f) for f in toks[2:6]]

                    last_valid = len(toks)
                    for _t in range(len(toks)-1, 5, -1):
                        if '\#' in toks[_t]:
                            last_valid = _t
                        elif toks[_t] in valid_feature_constraints:
                            local_feature[toks[_t]] = recast_globals([toks[_local_tok] for _local_tok in range(_t+1, last_valid)], globals)
                            last_valid = _t

                    features.append(local_feature)
                expert_line_index += 1

        valid_constituent_constraints = ['FIT', 'FITALL', 'DEPTH', 'DEPTHALL', 'DEPTH-FIT', 'FD', 'FDALL', 'FD-FIT', 'FD-DEPTH']
        if 'define constraints' in expert_file_text[expert_line_index]:
            expert_line_index += 1
            while ('endconstraint' not in expert_file_text[expert_line_index]):
                toks = expert_file_text[expert_line_index].strip().split()
                if toks[0] == 'constraint:':
                    last_valid = len(toks)
                    for _t in range(len(toks)-1, 0, -1):
                        if '\#' in toks[_t]:
                            last_valid = _t
                        elif np.any([vcc in toks[_t] for vcc in valid_constituent_constraints]):
                            if '<' in toks[_t]:
                                print('less than found in constraints, revise expert reader')
                            constraint_values = toks[_t].strip().split('>')
                            constraint_name = constraint_values.pop(0)
                            for _local_tok in range(_t+1, last_valid):
                                constraint_values.append(toks[_local_tok])
                            constituent_constraints[constraint_name] = recast_globals(constraint_values, globals)

                expert_line_index += 1
        expert_line_index += 1

    return decoded_expert


def read_mineral_fractions(file_list: List):
    """
    Read in a series of mineral fractions files from tetracorder, converting to dictionaries
    Args:
        file_list: list of files to read

    Returns:
        Dictionary keyed with unique file identifiers corresponding to expert system file
    """
    mineral_fractions = OrderedDict()
    mineral_names = [re.split('\.|-', os.path.basename(x))[0] for x in file_list]
    for _f, filename in enumerate(file_list):
        with open(filename, 'r') as fin:
            fractions_file_commented = fin.readlines()
        df = None
        for _line, line in enumerate(fractions_file_commented):
            if not line.strip().startswith('#'):
                df = pd.read_fwf(filename, skiprows=_line, header=None)
                break

        ## Hard coded due to inconsistent multi-line abundance file heading
        header = ['file','DN_scale','BD_factor','title','spectral_library','record']
        fraction_list = []
        if df is None:
            continue
        else:
            for idx in range(len(df[0])):
                local_entry = {}
                for headname, tok in zip(header, df.iloc[idx,:len(header)].tolist()):
                    local_entry[headname] = tok
                if type(local_entry[header[0]]) == float and np.isnan(local_entry[header[0]]):
                    continue
                fraction_list.append(local_entry)
                #if local_entry['spectral_library'] != 'splib06':
                #    print(local_entry)
        mineral_fractions[mineral_names[_f]] = fraction_list
    return mineral_fractions


def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('tetra_expert_file', type=str, metavar='TETRA_EXPERT_SYSTEM')
    parser.add_argument('converted_file', type=str, metavar='output_converted_tetrafile')
    parser.add_argument('-groups', type=int, default=DEFAULT_GROUPS, nargs='*')
    parser.add_argument('-sort_keys', type=int, default=int(DEFAULT_SORT_KEYS), choices=[0,1])
    parser.add_argument('-log_file', type=str, default=DEFAULT_LOGLEVEL)
    parser.add_argument('-log_level', type=str, default=DEFAULT_LOGLEVEL)
    args = parser.parse_args()

    args.sort_keys = args.sort_keys == 1

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    decoded_expert = decode_expert_system(args.tetra_expert_file, args.groups, args.log_file, log_level=args.log_level)

    #emit_utils.common_logs.logtime()

    #emit_utils.common_logs.logtime()
    with open(args.converted_file, 'w') as file:
        if os.path.splitext(args.converted_file)[-1] == '.yaml':
            outstring = yaml.dump(decoded_expert, sort_keys=args.sort_keys)
        elif os.path.splitext(args.converted_file)[-1] == '.json':
            outstring = json.dumps(decoded_expert, sort_keys=args.sort_keys, indent=4)
        file.write(outstring)


if __name__ == "__main__":
    main()
