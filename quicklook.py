# Philip G. Brodrick

from spectral.io import envi
import numpy as np
from osgeo import gdal
import argparse
import os
import subprocess



def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('input_file', type=str, metavar='aggregated abundance file')
    parser.add_argument('output_file', type=str, metavar='output file to write')
    args = parser.parse_args()


    band_names = ['calcite',  'chlorite' ,'dolomite', 'goethite' , 'gypsum' , 'hematite', 'illite+muscovite', 'kaolinite' , 'montmorillonite', 'vermiculite']
    
    ds = gdal.Open(args.input_file, gdal.GA_ReadOnly)
    dat = ds.ReadAsArray()
    
    maxband = np.argmax(dat, axis=0)
    
    colorlist=[\
    [255, 225, 25], #yellow
    [0, 130, 200], #blue
    [245, 130, 48], #orange
    [145, 30, 180], #purple
    [240, 50, 230], #magenta
    [250, 190, 212], #pink
    [170, 110, 40], #brown
    [128, 128, 0], #olive
    [128, 128, 128], #grey
    [255, 250, 200], #beige
    ]
    #[230, 25, 75], #red
    #60, 180, 75], #green
    #[70, 240, 240], #cyan
    #[210, 245, 60], #lime
    #[0, 128, 128], #teal
    #[220, 190, 255], #lavendar
    #[128, 0, 0], #maroon
    #[170, 255, 195], #mint
    #[255, 215, 180], #apricot
    #[0, 0, 128], #navy
    #[255, 255, 255], #white
    #[0, 0, 0], #black
    
    un_vals = np.unique(maxband)
    
    output_img = np.zeros((dat.shape[1], dat.shape[2], 3), dtype=np.uint8)
    leg_handles = []
    for _v, val in enumerate(un_vals):
        print(f'{_v}/{len(un_vals)-1})')
    
        subset = np.logical_and(maxband == val, dat[val,...] > 0)
    
        for _c in range(3):
            #output_img[subset,_c] = colorlist[val][_c] * ((dat[val, subset] - minval) / maxval)
            output_img[subset,_c] = colorlist[val][_c]
    
    del dat
    
    
    tmp_envi_file = os.path.splitext(args.output_file)[0]
    driver = gdal.GetDriverByName('ENVI')
    driver.Register()
    outDataset = driver.Create(tmp_envi_file,ds.RasterXSize,ds.RasterYSize,3,gdal.GDT_Byte)
    for _b in range(output_img.shape[-1]):
        outDataset.GetRasterBand(_b+1).WriteArray(output_img[...,_b])
    del outDataset
    subprocess.call(f'gdal_translate {tmp_envi_file} {args.output_file} -co ZLEVEL=9',shell=True)


if __name__ == "__main__":
    main()
