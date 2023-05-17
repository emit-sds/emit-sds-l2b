"""
This code contains support code for formatting L2B products for the LP DAAC.

Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
"""

import argparse
from netCDF4 import Dataset
from emit_utils.daac_converter import add_variable, makeDims, makeGlobalAttr, add_loc, add_glt
from emit_utils.file_checks import netcdf_ext, envi_header
from spectral.io import envi
import logging
import numpy as np
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L2B MIN PGE outputs to DAAC compatable formats, with supporting metadata''', add_help=True)

    parser.add_argument('output_abun_filename', type=str, help="Output abundance netcdf filename")
    parser.add_argument('output_abununcert_filename', type=str, help="Output abundance uncertainty netcdf filename")
    parser.add_argument('abun_file', type=str, help="EMIT L2B spectral abundance ENVI file")
    parser.add_argument('abun_unc_file', type=str, help="EMIT L2B spectral abundance uncertainty ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('glt_file', type=str, help="EMIT L1B glt ENVI file")
    parser.add_argument('version', type=str, help="3 digit (with leading V) version number")
    parser.add_argument('software_delivery_version', type=str, help="The extended build number at delivery time")
    parser.add_argument('--ummg_file', type=str, help="Output UMMG filename")
    parser.add_argument('--log_file', type=str, default=None, help="Logging file to write to")
    parser.add_argument('--log_level', type=str, default="INFO", help="Logging level")
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=args.log_level, filename=args.log_file)

    abun_ds = envi.open(envi_header(args.abun_file))
    abun_unc_ds = envi.open(envi_header(args.abun_unc_file))

    # make the netCDF4 file
    logging.info(f'Creating netCDF4 file: {args.output_abun_filename}')
    nc_ds = Dataset(args.output_abun_filename, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    makeGlobalAttr(nc_ds, args.abun_file, args.software_delivery_version, glt_envi_file=args.glt_file)

    nc_ds.title = "EMIT L2B Estimated Mineral Spectral Abundance 60 m " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis collection contains L2B band depth and geologic identification data. Band depth \
is estimated through linear feature matching - see ATBD for \
details. This collection includes band depth for both \'Group 1\' and \'Group 2\' minerals, which frequently co-occur. \
The band depth reported is that of the given mineral identified, which is also reported in a separate band. \
Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included."
    nc_ds.sync()

    logging.debug('Creating dimensions')
    makeDims(nc_ds, args.abun_file, args.glt_file)

    logging.debug('Creating and writing reflectance metadata')
    add_variable(nc_ds, "sensor_band_parameters/minerals", str, "Minerals", None,
                 abun_ds.metadata['band names'], {"dimensions": ("bands",)})

    logging.debug('Creating and writing location data')
    add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    add_glt(nc_ds, args.glt_file)

    logging.debug('Write spectral abundance data')
    add_variable(nc_ds, 'group_1_band_depth', "f4", "Group 1 Band Depth", "unitless", abun_ds.open_memmap(interleave='bip')[...,0].copy(),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_1_mineral_id', "i2", "Group 1 Mineral ID", "unitless", abun_ds.open_memmap(interleave='bip')[...,1].copy().astype(np.int16),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_band_depth', "f4", "Group 2 Band Depth", "unitless", abun_ds.open_memmap(interleave='bip')[...,2].copy(),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_mineral_id', "i2", "Group 2 Mineral ID", "unitless", abun_ds.open_memmap(interleave='bip')[...,3].copy().astype(np.int16),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    nc_ds.sync()
    logging.debug(f'Successfully created {args.output_abun_filename}')

    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'mineral_grouping_matrix_20230503.csv'))
    embed_keys = ['Index','Record','Name','URL','Group','Library']
    keytype = ['u4','u4',str,str,'u4',str]

    nc_ds.createDimension('minerals', len(df))

    for ek, ekt in zip(embed_keys, keytype):
        if ekt == str:
            converted_dat = np.array(df[ek]).astype('S')
        else:
            converted_dat = np.array(df[ek])
        add_variable(nc_ds, f'mineral_metadata/{ek.lower()}', ekt, ek, None, converted_dat, {"dimensions": ("minerals",)})
    nc_ds.sync()
    nc_ds.close()

    # make the netCDF4 file
    logging.info(f'Creating netCDF4 file: {args.output_abununcert_filename}')
    nc_ds = Dataset(args.output_abununcert_filename, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    makeGlobalAttr(nc_ds, args.abun_unc_file, args.software_delivery_version, glt_envi_file=args.glt_file)

    nc_ds.title = "EMIT L2B Estimated Mineral Spectral Abundance Uncertainty 60 m " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis collection contains L2B band depth uncertainty estimates of surface minerals, the fit quality of each mineral, \
and geolocation data. Band depth uncertainty is estimated by propogating reflectance uncertainty \
through linear feature matching used for abundance mapping - see ATBD for  details. \
Band depth uncertainty and fit qualities are provided for both \'Group 1\' and \'Group 2\' minerals, which frequently co-occur. \
Fit quality indicates how well the library-normalized observed spectra match the selected library spectra. \
Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included."
    nc_ds.sync()

    logging.debug('Creating dimensions')
    makeDims(nc_ds, args.abun_unc_file, args.glt_file)

    logging.debug('Creating and writing reflectance metadata')
    add_variable(nc_ds, "sensor_band_parameters/minerals", str, "Minerals", None,
                 abun_ds.metadata['band names'], {"dimensions": ("bands",)})

    logging.debug('Creating and writing location data')
    add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    add_glt(nc_ds, args.glt_file)

    add_variable(nc_ds, 'group_1_band_depth_unc', "f4", "Group 1 Band Depth Uncertainty", "unitless", abun_unc_ds.open_memmap(interleave='bip')[...,0].copy(),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_1_fit', "f4", "Group 1 Fit", "unitless", abun_unc_ds.open_memmap(interleave='bip')[...,1].copy(),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_band_depth_unc', "f4", "Group 2 Band Depth Uncertainty", "unitless", abun_unc_ds.open_memmap(interleave='bip')[...,2].copy(),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_fit', "f4", "Group 2 Fit", "unitless", abun_unc_ds.open_memmap(interleave='bip')[...,3].copy(),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})

    nc_ds.sync()
    nc_ds.close()
    logging.debug(f'Successfully created {args.output_abununcert_filename}')


    return


if __name__ == '__main__':
    main()
