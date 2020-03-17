import argparse
import gdal
import numpy as np
import sys 
from tqdm import tqdm


def main():

    # Read arguments
    parser = argparse.ArgumentParser(description="Scale and cast reflectance for tetracorder")
    parser.add_argument('infile_name', type=str)
    parser.add_argument('outfile_name', type=str)
    parser.add_argument('-nodata_value', type=float, default=-9999)
    parser.add_argument('-output_nodata_value', type=int, default=-1)
    parser.add_argument('-scaling_factor', type=float, default=20000)
    args = parser.parse_args()
    
    
    # Open source dataset
    dataset = gdal.Open(args.infile_name,gdal.GA_ReadOnly)
    
    # Create output dataset
    driver = gdal.GetDriverByName('ENVI')
    driver.Register()
    
    outDataset = driver.Create(args.outfile_name, dataset.RasterXSize, dataset.RasterYSize, dataset.RasterCount, gdal.GDT_Int16, options=['INTERLEAVE=BIL'])
    outDataset.SetProjection(dataset.GetProjection())
    outDataset.SetGeoTransform(dataset.GetGeoTransform())
    del outDataset
    
    
    # Run through line by line
    for _line in tqdm(range(dataset.RasterYSize),ncols=80):
        # Load line
        dat = dataset.ReadAsArray(0,_line,dataset.RasterXSize,1)
    
        # find where data are bad
        mask = dat == args.nodata_value
    
        # convert data to 16-bit int, scaled by scaling factor
        dat = np.round(dat*args.scaling_factor).astype(np.int16)
    
        # re-assign nodata value, and swap axes for BIL write
        dat[mask] = args.output_nodata_value
        dat = np.swapaxes(dat,0,1)
    
        # open memmap for output writing
        memmap = np.memmap(args.outfile_name,mode='r+',shape=(dataset.RasterYSize,dataset.RasterCount,
                           dataset.RasterXSize),dtype=np.int16)
    
        # write data and flush cache
        memmap[_line:_line+1,...] = dat 
        del memmap

if __name__ == "__main__":
    main()
