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


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L2B PGE outputs to DAAC compatable formats, with supporting metadata''', add_help=True)

    parser.add_argument('output_filename', type=str, help="Output netcdf filename")
    parser.add_argument('abun_file', type=str, help="EMIT L2B spectral abundance ENVI file")
    parser.add_argument('abun_unc_file', type=str, help="EMIT L2B spectral abundance uncertainty ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('glt_file', type=str, help="EMIT L1B glt ENVI file")
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
    logging.info(f'Creating netCDF4 file: {args.output_filename}')
    nc_ds = Dataset(args.output_filename, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    makeGlobalAttr(nc_ds, args.abun_file, args.glt_file)

    nc_ds.title = "EMIT L2B Estimated Mineral Spectral Abundance and Uncertainty 60 m V001"
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis collection contains L2B spectral abundance estimates of surface mineralology \
        and geolocation data. Spectral abundance is estimated through linear feature matching - see ATBD for  \
        details. Spectral abundance values are reported as fractions (relative to 1). "
    nc_ds.sync()

    logging.debug('Creating dimensions')
    makeDims(nc_ds, args.abun_file, args.glt_file)

    logging.debug('Creating and writing reflectance metadata')
    add_variable(nc_ds, "sensor_band_parameters/minerals", str, "Minerals", None,
                 abun_ds.metadata['band names'], {"dimensions": ("number_of_bands",)})

    logging.debug('Creating and writing location data')
    add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    add_glt(nc_ds, args.glt_file)

    logging.debug('Write spectral abundance data')
    add_variable(nc_ds, 'spectral_abundance', "f4", "Spectral Abundance", "unitless", abun_ds.open_memmap(interleave='bip')[...].copy(),
                 {"dimensions":("number_of_scans", "pixels_per_scan", "number_of_bands")})
    nc_ds.sync()
    add_variable(nc_ds, 'spectral_abundance_uncertainty', "f4", "Spectral Abundance Uncertainty", "unitless",
                 abun_unc_ds.open_memmap(interleave='bip')[...].copy(),
                 {"dimensions":("number_of_scans", "pixels_per_scan", "number_of_bands")})
    nc_ds.sync()

    nc_ds.close()
    logging.debug(f'Successfully created {args.output_filename}')

    return


if __name__ == '__main__':
    main()
