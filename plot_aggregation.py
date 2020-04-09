




import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gdal
import gzip
import numpy as np
import math

import spectral.io.envi as envi
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, os

#import matplotlib.font_manager
#names = matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
#for name in names:
#    print(name)
#quit()

#from matplotlib import rcParams
#plt.rcParams["font.family"] = "serif"
#plt.rcParams["font.serif"] = "Times New Roman"


aggregated_file = sys.argv[1]
tetra_file_list = sys.argv[2]
files = np.genfromtxt(tetra_file_list,dtype='str').tolist()

aggregated_dataset = gdal.Open(aggregated_file,gdal.GA_ReadOnly)
rows = aggregated_dataset.RasterYSize
cols = aggregated_dataset.RasterXSize
aggregated = np.memmap(aggregated_file, mode='r', shape=(rows, aggregated_dataset.RasterCount, cols), dtype=np.float32)

use_rows = cols
use_rows = 500
aggregated = aggregated[:use_rows,-2,:]
image_stack = []
names = []

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
    
    band_depth = band_depth[:use_rows,:]
    if np.sum(band_depth != 0) / np.product(band_depth.shape) > 0.0005:
        image_stack.append(band_depth)
        print(file)
        print(np.sum(band_depth != 0) / np.product(band_depth.shape))
        names.append(os.path.splitext(os.path.basename(file))[0])


fig = plt.figure(figsize=(6,6))
spec = gridspec.GridSpec(ncols=2, nrows=3, figure=fig, left=0.05, wspace=0.05)

for _i in range(3):

    bd = image_stack[_i]
    ax = fig.add_subplot(spec[_i,0])
    print(np.min(bd[bd!=0]))
    print(np.max(bd[bd!=0]))

    ax.imshow(bd,cmap='hot')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    #tit = ''
    #for subtit in names[_i].split('+'):
    #    tit += subtit + '\n'
    #ax.set_title(tit)
    #ax.set_title(names[_i],fontsize=6)
    #ax.set_title('Reference Library {}'.format(_i+1),fontsize=8)
    

ax = fig.add_subplot(spec[:,1])
bd = aggregated.copy()
print(np.min(bd[bd!=0]))
print(np.max(bd[bd!=0]))

im = ax.imshow(bd,cmap='hot')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xticks([])
ax.set_yticks([])
#ax.set_title('Aggregated')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = ax.figure.colorbar(im,cax=cax)

    
plt.savefig('figs/l2b_aggregation_example.png',dpi=300,bbox_inches='tight')
    
    
