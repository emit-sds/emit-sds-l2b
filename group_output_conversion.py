"""
This code contains support code for formatting L2B products for the LP DAAC.

Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
"""

import argparse
from netCDF4 import Dataset
from emit_utils.daac_converter import add_variable, makeDims, makeGlobalAttrBase, add_loc, add_glt, get_spatial_extent_res
from emit_utils.file_checks import netcdf_ext, envi_header
from osgeo import gdal
from spectral.io import envi
import logging
import numpy as np
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L2B MIN PGE outputs to DAAC compatable formats, with supporting metadata''', add_help=True)

    parser.add_argument('output_abun_file', type=str, help="Output abundance netcdf filename")
    parser.add_argument('output_abununcert_file', type=str, help="Output abundance uncertainty netcdf filename")
    parser.add_argument('abun_file', type=str, help="EMIT L2B spectral abundance NetCDF file")
    parser.add_argument('abununcert_file', type=str, help="EMIT L2B spectral abundance uncertainty NetCDF file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('glt_file', type=str, help="EMIT L1B glt ENVI file")
    parser.add_argument('mineral_grouping_file', type=str, help="Path to mineral grouping file")
    parser.add_argument('--start_time', type=str, help="Start time of the acquisition", default="2020-01-01T00:00:00")
    parser.add_argument('--stop_time', type=str, help="Stop time of the acquisition", default="2020-01-01T00:00:00")
    parser.add_argument('--version', type=str, help="3 digit (with leading V) version number", default="V000")
    parser.add_argument('--software_build_version', type=str, help="The extended build number when the product was created", default="000000")
    parser.add_argument('--history_file', type=str, help="File containing NetCDF history attribute - concatenated run command and input files list", default="")
    parser.add_argument('--daynight', type=str, help="Value of Day/Night flag", default="Day")
    parser.add_argument('--ummg_file', type=str, help="Output UMMG filename")
    parser.add_argument('--rfl_file', type=str, help="Reflectance file for masking only (any file with -9999 in the correct positions will work)")
    parser.add_argument('--log_file', type=str, default=None, help="Logging file to write to")
    parser.add_argument('--log_level', type=str, default="INFO", help="Logging level")
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=args.log_level, filename=args.log_file)

    abun_ds = Dataset(args.abun_file, "r")
    abununcert_ds = Dataset(args.abununcert_file, "r")

    # make the netCDF4 file
    logging.info(f'Creating netCDF4 file: {args.output_abun_file}')
    nc_ds = Dataset(args.output_abun_file, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    makeGlobalAttrBase(nc_ds)

    logging.debug('Fetch reflectance-based mask')
    rfl_mask = envi.open(envi_header(args.rfl_file)).open_memmap(interleave='bip')[:, :, 0] == -9999

    # Add scene specific attributes
    nc_ds.flight_line = os.path.basename(args.output_abun_file)[:19]
    nc_ds.time_coverage_start = args.start_time
    nc_ds.time_coverage_end = args.stop_time
    nc_ds.software_build_version = args.software_build_version
    nc_ds.product_version = args.version
    history = "None specified"
    if len(args.history_file) > 0:
        with open(args.history_file, "r") as f:
            history = f.read()
    nc_ds.history = history

    # Add spatial extent
    ul_lr, res = get_spatial_extent_res(args.glt_file)
    nc_ds.easternmost_longitude = ul_lr[2]
    nc_ds.northernmost_latitude = ul_lr[1]
    nc_ds.westernmost_longitude = ul_lr[0]
    nc_ds.southernmost_latitude = ul_lr[3]
    nc_ds.spatialResolution = res

    glt_ds = gdal.Open(args.glt_file)
    nc_ds.spatial_ref = glt_ds.GetProjection()
    nc_ds.geotransform = glt_ds.GetGeoTransform()

    nc_ds.day_night_flag = args.daynight

    nc_ds.title = "EMIT L2B Estimated Mineral Identification and Band Depth 60 m " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis collection contains L2B band depth and geologic identification data. Band depth \
is estimated through linear feature matching - see ATBD for \
details. This collection includes band depth for both \'Group 1\' and \'Group 2\' minerals, which frequently co-occur. \
The band depth reported is that of the given mineral identified, which is also reported in a separate band. \
Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included."
    nc_ds.sync()

    logging.debug('Creating dimensions')
    nc_ds.createDimension('downtrack', len(abun_ds.dimensions['downtrack']))
    nc_ds.createDimension('crosstrack', len(abun_ds.dimensions['crosstrack']))
    nc_ds.createDimension('ortho_y', glt_ds.RasterYSize)
    nc_ds.createDimension('ortho_x', glt_ds.RasterXSize)

    logging.debug('Creating and writing location data')
    add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    add_glt(nc_ds, args.glt_file)

    logging.debug('Load and mask abundance')
    def mask_key(ds, key):
        data = ds.variables[key][:].copy()
        data[rfl_mask,:] = -9999
        return data

    logging.debug('Write spectral abundance data')
    add_variable(nc_ds, 'group_1_band_depth', "f4", "Group 1 Band Depth", "unitless", mask_key(abun_ds, 'group_1_band_depth'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_1_mineral_id', "i2", "Group 1 Mineral ID", "unitless", mask_key(abun_ds, 'group_1_mineral_id'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_band_depth', "f4", "Group 2 Band Depth", "unitless", mask_key(abun_ds, 'group_2_band_depth'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_mineral_id', "i2", "Group 2 Mineral ID", "unitless", mask_key(abun_ds, 'group_2_mineral_id'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    nc_ds.sync()
    logging.debug(f'Successfully created {args.output_abun_file}')

    logging.debug("Embedding mineral metadata")
    ref_df = pd.read_csv(args.mineral_grouping_file)
    ref_df['path_name'] = ref_df['path'].apply(os.path.basename)
    nc_ds.createDimension("minerals", len(ref_df))
    METADATA_COLUMNS = [
        ("index", "index", "u4"),
        ("record", "record", "u4"),
        ("name", "path_name", str),
        ("sample_name", "title", str),
        ("url", "url", str),
        ("group", "group", "u4"),
        ("library", "library", str),
    ]
    for name, column, dtype in METADATA_COLUMNS:
        if column not in ref_df.columns:
            logging.warning(f"Reference matrix has no '{column}' column, skipping mineral_metadata/{name}")
            continue
        if dtype is str:
            data = np.array(ref_df[column]).astype("S")
        else:
            data = np.array(ref_df[column])
        add_variable(nc_ds, f"mineral_metadata/{name}", dtype, column, None, data, {"dimensions": ("minerals",)})

    # df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'mineral_grouping_matrix_20230503.csv'))
    # embed_keys = ['Index','Record','Name','URL','Group','Library']
    # keytype = ['u4','u4',str,str,'u4',str]

    # nc_ds.createDimension('minerals', len(df))

    # for ek, ekt in zip(embed_keys, keytype):
    #     if ekt == str:
    #         converted_dat = np.array(df[ek]).astype('S')
    #     else:
    #         converted_dat = np.array(df[ek])
    #     add_variable(nc_ds, f'mineral_metadata/{ek.lower()}', ekt, ek, None, converted_dat, {"dimensions": ("minerals",)})
    nc_ds.sync()
    nc_ds.close()

    # make the netCDF4 uncertainty file
    logging.info(f'Creating netCDF4 file: {args.output_abununcert_file}')
    nc_ds = Dataset(args.output_abununcert_file, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    makeGlobalAttrBase(nc_ds)

    # Add scene specific attributes
    nc_ds.flight_line = os.path.basename(args.output_abun_file)[:19]
    nc_ds.time_coverage_start = args.start_time
    nc_ds.time_coverage_end = args.stop_time
    nc_ds.software_build_version = args.software_build_version
    nc_ds.product_version = args.version
    nc_ds.history = history
    
    # Add spatial extent
    # ul_lr, res = get_spatial_extent_res(args.glt_file)
    nc_ds.easternmost_longitude = ul_lr[2]
    nc_ds.northernmost_latitude = ul_lr[1]
    nc_ds.westernmost_longitude = ul_lr[0]
    nc_ds.southernmost_latitude = ul_lr[3]
    nc_ds.spatialResolution = res
    nc_ds.spatial_ref = glt_ds.GetProjection()
    nc_ds.geotransform = glt_ds.GetGeoTransform()
    nc_ds.day_night_flag = args.daynight

    nc_ds.title = "EMIT L2B Estimated Mineral Identification and Band Depth Uncertainty 60 m " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis collection contains L2B band depth uncertainty estimates of surface minerals, the fit quality of each mineral, \
and geolocation data. Band depth uncertainty is estimated by propogating reflectance uncertainty \
through linear feature matching used for abundance mapping - see ATBD for  details. \
Band depth uncertainty and fit qualities are provided for both \'Group 1\' and \'Group 2\' minerals, which frequently co-occur. \
Fit quality indicates how well the library-normalized observed spectra match the selected library spectra. \
Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included."
    nc_ds.sync()

    logging.debug('Creating dimensions')
    nc_ds.createDimension('downtrack', len(abununcert_ds.dimensions['downtrack']))
    nc_ds.createDimension('crosstrack', len(abununcert_ds.dimensions['crosstrack']))
    nc_ds.createDimension('ortho_y', glt_ds.RasterYSize)
    nc_ds.createDimension('ortho_x', glt_ds.RasterXSize)

    logging.debug('Creating and writing location data')
    add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    add_glt(nc_ds, args.glt_file)

    add_variable(nc_ds, 'group_1_band_depth_unc', "f4", "Group 1 Band Depth Uncertainty", "unitless", mask_key(abununcert_ds, 'group_1_band_depth_unc'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_1_fit', "f4", "Group 1 Fit", "unitless", mask_key(abununcert_ds, 'group_1_fit'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_band_depth_unc', "f4", "Group 2 Band Depth Uncertainty", "unitless", mask_key(abununcert_ds, 'group_2_band_depth_unc'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})
    add_variable(nc_ds, 'group_2_fit', "f4", "Group 2 Fit", "unitless", mask_key(abununcert_ds, 'group_2_fit'),
                 {"dimensions":("downtrack", "crosstrack"), "zlib": True, "complevel": 9})

    nc_ds.sync()
    nc_ds.close()
    logging.debug(f'Successfully created {args.output_abununcert_file}')


    return


if __name__ == '__main__':
    main()
