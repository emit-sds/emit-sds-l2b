




import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import gdal
import gzip
import numpy as np
from matplotlib.patches import Patch
import math

import spectral.io.envi as envi
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, os
from spectral import envi
import argparse
import subprocess
from scipy.signal import savgol_filter, medfilt

from matplotlib import rcParams
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"

parser = argparse.ArgumentParser(description='DEMO L3 aggregation')
parser.add_argument('aggregated_file')
parser.add_argument('tetra_base_dir')
parser.add_argument('refl_file')
parser.add_argument('-mineral',default='calcite')
parser.add_argument('-uncertainty_file',default=None)
parser.add_argument('-simple_names',default=0,type=int)
args = parser.parse_args()



tetra_file_list = 'file_lists/{}_{}.txt'.format(os.path.basename(args.tetra_base_dir),args.mineral)
cmd_str = 'find {} -iname "*{}*" | grep gz | grep depth | grep -v .hdr > {}'.format(args.tetra_base_dir, args.mineral, tetra_file_list)
subprocess.call(cmd_str,shell=True)


files = np.genfromtxt(tetra_file_list,dtype='str').tolist()

agg_hdr = envi.read_envi_header(args.aggregated_file + '.hdr')
mineral_band_names = agg_hdr['band names']

aggregated_dataset = gdal.Open(args.aggregated_file,gdal.GA_ReadOnly)
rows = aggregated_dataset.RasterYSize
cols = aggregated_dataset.RasterXSize
aggregated = np.memmap(args.aggregated_file, mode='r', shape=(rows, aggregated_dataset.RasterCount, cols), dtype=np.float32)

refl_dataset = gdal.Open(args.refl_file,gdal.GA_ReadOnly)
#refl = np.memmap(args.refl_file, mode='r', shape=(rows, refl_dataset.RasterCount, cols), dtype=np.int16)
refl = envi.open(args.refl_file + '.hdr').open_memmap()
#rgb = refl[:,np.array([12,21,30]),:].copy().transpose((0,2,1)).astype(float)/20000.
rgb = refl[:,:,np.array([54 ,36 ,18])].copy().astype(float )


use_rows = cols
#use_rows = [0,-1]
#use_rows = [0,5000]
#use_rows = [16000, 19500]
use_rows = [13500, 17000]
aggregated = aggregated[use_rows[0]:use_rows[1], :, :]
rgb = rgb[use_rows[0]:use_rows[1], :, :]

image_stack = []
names = []
fraction = []

for file in files:

    hdr = envi.read_envi_header(file + '.hdr')
    offs = int(hdr['header offset'])

    with open(file, 'rb') as fin:
        compressed = fin.read()
    decompressed = gzip.decompress(compressed)
    

    data_type_scaling = 0.5 #Not global, will work for viz

    band_depth = np.frombuffer(decompressed, dtype=np.uint8, count=(rows*cols)+offs)
    band_depth = band_depth[offs:]  # one line offset by convention?
    band_depth = band_depth.reshape((rows, cols))
    band_depth = band_depth.astype(dtype=np.float32) / 255.0 * data_type_scaling

    band_depth = band_depth[use_rows[0]:use_rows[1],:]
    image_stack.append(band_depth)
    names.append(os.path.splitext(os.path.basename(file))[0])
    fraction.append(np.sum(band_depth != 0) / float(np.product(band_depth.shape)))

top_fractions = np.argsort(fraction)[::-1]
print('Top Fractions: {}'.format(top_fractions[:3]))
image_stack = [image_stack[x] for x in top_fractions[:3]]
names = [names[x] for x in top_fractions]
print(names)

bd_ind = -1
for ind, name in enumerate(mineral_band_names):
    if args.mineral in name.lower():
        bd_ind = ind
if bd_ind == -1:
    print('Could not find mineral: {} in: {}'.format(args.mineral, mineral_band_names))
    quit()
bd = aggregated[:,bd_ind,:].copy()



figsize = (11,8)
fig = plt.figure(figsize=figsize)
spec = gridspec.GridSpec(ncols=len(image_stack)+4, nrows=2, figure=fig, left=0.05, wspace=0.05)

ax = fig.add_subplot(spec[:, 0])
#rgb = rgb[:,:,::-1]
rgb -= np.percentile(rgb,2,axis=(0,1))[np.newaxis,np.newaxis,:]
rgb /= np.percentile(rgb,98,axis=(0,1))[np.newaxis,np.newaxis,:]
ax.imshow(rgb)
plt.axis('off')
ax.set_title('RGB')

for _i in range(len(image_stack)):
    bd = image_stack[_i]
    ax = fig.add_subplot(spec[:,_i+1])

    ax.imshow(bd != 0,cmap='inferno',vmin=0,vmax=1)
    plt.xticks([],[])
    plt.yticks([],[])
    tit = ''
    print(names[_i])
    for subtit in names[_i].split('+'):
        tit += subtit + '\n'
    if args.simple_names == 1:
        ax.set_title(mineral_band_names[bd_ind] + '\nMixture {}'.format(_i+1))
    else:
        if _i % 2 == 1:
            ax.set_title(tit[:-2])
        else:
            ax.set_xlabel(tit[:-2])

ax = fig.add_subplot(spec[:,-3])
ag_im = ax.imshow(aggregated[:,bd_ind,:], cmap='inferno', vmin=0, vmax=np.percentile(aggregated[...,bd_ind],99))
plt.axis('off')
ax.set_title('{}\nSpectral\nAbundance'.format(mineral_band_names[bd_ind]))


ax = fig.add_subplot(spec[:,-2])
most_bands = np.nansum(aggregated!=0,axis=(0,2))
most_bands = np.argsort(most_bands)[::-1]
bd = np.transpose(aggregated[:,most_bands[:3],:],(0,2,1))
for _b in range(bd.shape[-1]):
    bd[...,_b] -= np.percentile(bd[bd[...,_b] != 0,_b],2)
    bd[...,_b] /= np.percentile(bd[bd[...,_b] != 0,_b],98)
plt.imshow(bd)
plt.axis('off')

ax.set_title('3-Mineral\nComposite')
mineral_leg_handles = [Patch(facecolor='red', edgecolor='black',label=mineral_band_names[most_bands[0]]),
                       Patch(facecolor='green', edgecolor='black',label=mineral_band_names[most_bands[1]]),
                       Patch(facecolor='blue', edgecolor='black',label=mineral_band_names[most_bands[2]])]
ax = fig.add_subplot(spec[0,-1])
ax.legend(handles=mineral_leg_handles, loc='center', ncol=1, frameon=False)
plt.axis('off')

ax = fig.add_subplot(spec[1,-1])
plt.axis('off')

divider = make_axes_locatable(ax)
cax = divider.append_axes("left", size="30%")
cbar = plt.colorbar(ag_im, cax=cax)

cbar.set_label('{}\nSpectral Abundance'.format(mineral_band_names[bd_ind]))

plt.savefig('figs/aggregation_example_{}_{}'.format(os.path.basename(args.tetra_base_dir),args.mineral),dpi=300,bbox_inches='tight')
#plt.show()
plt.clf()


if args.uncertainty_file is not None:
    fig = plt.figure(figsize=figsize)
    spec = gridspec.GridSpec(ncols=len(image_stack) + 5, nrows=10, figure=fig, left=0.05, wspace=0.15)

    uncertainty_dataset = gdal.Open(args.uncertainty_file, gdal.GA_ReadOnly)
    rows = uncertainty_dataset.RasterYSize
    cols = uncertainty_dataset.RasterXSize
    uncertainty = np.memmap(args.uncertainty_file, mode='r', shape=(rows, uncertainty_dataset.RasterCount, cols),dtype=np.float32)
    uncertainty = uncertainty[use_rows[0]:use_rows[1], :, :]

    ax = fig.add_subplot(spec[:, 0])
    rgb = rgb[:, :, ::-1]
    rgb -= np.percentile(rgb, 2, axis=(0, 1))[np.newaxis, np.newaxis, :]
    rgb /= np.percentile(rgb, 98, axis=(0, 1))[np.newaxis, np.newaxis, :]
    ax.imshow(rgb)
    plt.axis('off')
    ax.set_title('RGB')

    ax = fig.add_subplot(spec[:, 1])
    bd = aggregated[:, bd_ind, :].copy()
    ub = np.percentile(bd[bd != 0], 98)
    ag_im = ax.imshow(bd, cmap='inferno', vmin=np.percentile(bd, 0), vmax=ub)
    plt.axis('off')
    ax.set_title('{}\nSpectral Abundance'.format(mineral_band_names[bd_ind]))

    mask = bd == 0

    ax = fig.add_subplot(spec[:, 2])
    bd = uncertainty[:,bd_ind,:].copy()
    bd[mask] = 0
    uncert_im = ax.imshow(bd/4, cmap='inferno',vmin=0,vmax=ub/3.)
    plt.axis('off')
    ax.set_title('{}\nSpectral Abundance\nUncertainty'.format(mineral_band_names[bd_ind]))


    ax = fig.add_subplot(spec[:,4])
    most_bands = np.nansum(aggregated!=0,axis=(0,2))
    most_bands = np.argsort(most_bands)[::-1]
    bd = np.transpose(aggregated[:,most_bands[:3],:],(0,2,1))
    for _b in range(bd.shape[-1]):
        bd[...,_b] -= np.percentile(bd[bd[...,_b] != 0,_b],2)
        bd[...,_b] /= np.percentile(bd[bd[...,_b] != 0,_b],98)
        bd[...,_b] = medfilt(bd[...,_b],(1,9))
        bd[...,_b] = np.minimum(bd[...,_b], savgol_filter(bd[...,_b], window_length=9, polyorder=3,axis=1))
    plt.imshow(bd)
    plt.axis('off')

    ax.set_title('3-Mineral\nComposite')
    mineral_leg_handles = [Patch(facecolor='red', edgecolor='black',label=mineral_band_names[most_bands[0]]),
                           Patch(facecolor='green', edgecolor='black',label=mineral_band_names[most_bands[1]]),
                           Patch(facecolor='blue', edgecolor='black',label=mineral_band_names[most_bands[2]])]
    ax = fig.add_subplot(spec[1,3])
    ax.legend(handles=mineral_leg_handles, loc='center', ncol=1, frameon=False)
    plt.axis('off')






    ax = fig.add_subplot(spec[6:9, 3])
    plt.axis('off')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("left", size="30%")
    cbar = plt.colorbar(ag_im, cax=cax)
    cbar.set_label('{}\nSpectral Abundance'.format(mineral_band_names[bd_ind]),fontsize=12)
    #ax.set_title('{}\nSpectral\nAbundance'.format(mineral_band_names[bd_ind]),fontsize=12)

    ax = fig.add_subplot(spec[3:6, 3])
    plt.axis('off')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("left", size="30%")
    cbar = plt.colorbar(uncert_im, cax=cax)
    cbar.set_label('{} Spectral\nAbundance Uncertainty'.format(mineral_band_names[bd_ind]),fontsize=12)
    #ax.set_title('{}\nSpectral\nAbundance\nUncertainty'.format(mineral_band_names[bd_ind]),fontsize=12)


    #ax = fig.add_subplot(spec[:,-2])
    #bd = np.transpose(uncertainty[:,most_bands[:3],:],(0,2,1))
    #bd -= np.percentile(bd,0,axis=(0,1))[np.newaxis,np.newaxis,:]
    #bd /= np.percentile(bd,100,axis=(0,1))[np.newaxis,np.newaxis,:]
    #plt.imshow(bd)
    #plt.axis('off')

    #ax.set_title('3-Mineral\nComposite')
    #mineral_leg_handles = [Patch(facecolor='red', edgecolor='black',label=mineral_band_names[most_bands[0]]),
    #                       Patch(facecolor='green', edgecolor='black',label=mineral_band_names[most_bands[1]]),
    #                       Patch(facecolor='blue', edgecolor='black',label=mineral_band_names[most_bands[2]])]
    #ax = fig.add_subplot(spec[0,-1])
    #ax.legend(handles=mineral_leg_handles, loc='center', ncol=1, frameon=False)
    #plt.axis('off')
    plt.savefig('figs/aggregation_example_{}_{}_uncert'.format(os.path.basename(args.tetra_base_dir),args.mineral),dpi=300,bbox_inches='tight')
    #plt.show()

