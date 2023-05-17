"""
Quicklook generation from L2bmin products.

Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
"""

from spectral.io import envi
import numpy as np
from osgeo import gdal
import argparse
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from emit_utils.file_checks import envi_header

plt.switch_backend('agg')


def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('input_file', type=str, metavar='l2b file')
    parser.add_argument('output_file', type=str, metavar='output file to write')
    parser.add_argument('--unc_file', type=str, metavar='uncertainty file')
    args = parser.parse_args()

    ds = envi.open(envi_header(args.input_file))
    dat = ds.open_memmap(interleave='bip')
    dat[dat == -9999] = np.nan
    if args.unc_file is not None:
        unc_ds = envi.open(envi_header(args.unc_file))
        unc = unc_ds.open_memmap(interleave='bip')
        unc[unc == -9999] = np.nan

    fig = plt.figure(figsize=(20,20)) 
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

    ax = plt.subplot(gs[0,0])
    im = plt.imshow(dat[...,0], vmin=0, vmax=0.25)
    plt.xlabel('Crosstrack')
    plt.ylabel('Downtrack')
    plt.title('Group 1 Band Depth')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    ax = plt.subplot(gs[0,1])
    im = plt.imshow(dat[...,2], vmin=0, vmax=0.25)
    plt.xlabel('Crosstrack')
    plt.ylabel('Downtrack')
    plt.title('Group 2 Band Depth')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)


    if args.unc_file is not None:
        ax = plt.subplot(gs[1,0])
        im = plt.imshow(unc[...,0], vmin=0, vmax=0.025)
        plt.xlabel('Crosstrack')
        plt.ylabel('Downtrack')
        plt.title('Group 1 Band Depth Uncertainty')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)

        ax = plt.subplot(gs[1,1])
        im = plt.imshow(unc[...,2], vmin=0, vmax=0.025)
        plt.xlabel('Crosstrack')
        plt.ylabel('Downtrack')
        plt.title('Group 2 Band Depth Uncertainty')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)

    plt.savefig(args.output_file, dpi=400, bbox_inches='tight')


if __name__ == "__main__":
    main()
