# David Thompson


import os
import sys, gzip
import argparse
from scipy import logical_and as aand
import scipy as s
import spectral
import spectral.io.envi as envi
from scipy.interpolate import interp1d

def band_depth(wl, reflectance, feature):
    """ Feature is a four-wavelength tuple defining the start and end of two
        averaging windows used for continnum interpolation"""
    left_inds = aand(wl>=feature[0], wl<=feature[1])
    left_x = wl[int(left_inds.mean())]
    left_y = reflectance[left_inds].mean()
    right_inds = aand(wl>=feature[2], wl<=feature[3])
    right_x = wl[int(right_inds.mean())]
    right_y = reflectance[right_inds].mean()
    ctm = interp1d([left_x, right_x],[left_y, right_y])(wl) 
    feature_inds = aand(wl>=feature[0], wl<=feature[3])
    depths = 1.0-rfl[feature_inds]/wl[feature_inds]
    return max(depths)


# parse the command line (perform the correction on all command line arguments)
def main():

  parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
  parser.add_argument('expert_system_file', type=str, metavar='EXPERT_SYSTEM')
  parser.add_argument('library_file',type=str, metavar='LIBRARY_FILE')
  parser.add_argument('tetracorder_output_base',type=str, metavar='TETRA_OUTPUT_DIR')
  parser.add_argument('output', type=str, metavar='OUTPUT')
  args = parser.parse_args()

  lib         = envi.open(args.library_file+'.hdr', args.library_file)
  lib_rfl     = lib.spectra.copy()
  lib_records = [int(q) for q in lib.metadata['record']]
  names       = [q.strip() for q in lib.metadata['spectra names']]
  wl          = s.array([float(q) for q in lib.metadata['wavelength']])

  nminerals = 10
  out_data = None

  with open(args.expert_system_file,'r') as fin:
    lines = fin.readlines()

  # Go through expert system file one line at a time
  i, group, spectrum, output_data, header = 0, None, None, None, True
  while i<len(lines):

      if lines[i].startswith('BEGIN SETUP'):
          header = False

      # ignore comments
      if header or lines[i].startswith('\#'):
          i = i+1
          continue

      if lines[i].startswith('group'):
          
          if group is not None:
             
             # We are finalizing this spectrum.
             # get band depth of spectrum library file
             bd_library = band_depth(wl, rfl, features[0])
             
             # get band depths from map.  
             groupdir = args.tetracorder_output_base+'/group.'+int(group)+'um/'
             hdr = envi.read_envi_header(groupdir+filename+'.depth.gz.hdr')
             cols, rows = hdr['samples'], hdr['lines']
             if out_data is None:

                out_hdr = copy(hdr)
                out_hdr['interleave'] = 'bil'
                out_hdr['data type'] = 4
                out_hdr['wavelengths'] = '{'+','.join([str(q) for q in wl])+'}'
                out_hdr['bands'] = 10
                samples = int(hdr['samples']) 
                lines = int(hdr['lines']) 
                out_data = s.zeros((rows, cols, nminerals))

             mapfile = gzip.open(groupdir+filename+'.depth.gz')
             bd = s.fromfile(mapfile,dtype=s.uint8,count=rows*cols)
             bd = bd.reshape((rows,cols))
 
             # normalize to the depth of the library spectrum,
             # translating to aerial fractions
             bd_map = bd_map / bd_library 
             
             # add to an output channel
             for key, chans in chanmap.items():
                 if key in filename:
                     my_chans = chans.copy()
                     my_chans = my_chans.reshape(1,1,nminerals)
             out_data = out_data + bd.reshape((rows,cols,1)) @ my_chans

          group = int(lines[i].strip().split()[-1])

      if 'SMALL' in lines[i]:
          record = int(lines[i].strip().split()[3])
          rfl = lib_rfl[lib_records.index(record),:]

      if 'define output' in lines[i]:
          filename = lines[i+2].strip().split()[0]
          bd_scaling = float(lines[i+3].strip().split()[-1])

      if 'define features' in lines[i]:
          j= i+1
          features = []
          while ('endfeatures' not in lines[j]):
              if lines[j].strip().startswith('f'):
                  features.append([float(f) for f in 
			lines[j].strip().split()[2:6]])
              j = j + 1
      i = i + 1
    
  out_img = envi.create_image(args.output+'.hdr', meta=out_hdr, force=True,
                    shape=(rows,cols,nminerals), dtype=np.float32, ext='')
  mm = out_img.open_memmap(writable=True)
  for r in range(rows):
      mm[r,:,:] = out_data[r,:,:]
 
if __name__ == "__main__":
  main()



