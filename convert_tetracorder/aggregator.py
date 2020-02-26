# David Thompson


import argparse
import gzip
import numpy as np
from scipy.interpolate import interp1d
import spectral.io.envi as envi

def band_depth(wl, reflectance, feature):
    """ Feature is a four-wavelength tuple defining the start and end of two
        averaging windows used for continnum interpolation"""
    left_inds = np.where(np.logical_and(wl>=feature[0], wl<=feature[1]))[0]
    left_x = wl[int(left_inds.mean())]
    left_y = reflectance[left_inds].mean()
    right_inds = np.where(np.logical_and(wl>=feature[2], wl<=feature[3]))[0]
    right_x = wl[int(right_inds.mean())]
    right_y = reflectance[right_inds].mean()
    ctm = interp1d([left_x, right_x],[left_y, right_y], 
            bounds_error=False, fill_value='extrapolate')(wl) 
    feature_inds = np.logical_and(wl>=feature[0], wl<=feature[3])
    depths = 1.0-reflectance[feature_inds]/ctm[feature_inds]
    return max(depths)


# parse the command line (perform the correction on all command line arguments)
def main():

  parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
  parser.add_argument('tetra_expert_file', type=str, metavar='TETRA_EXPERT_SYSTEM')
  parser.add_argument('tetra_library_file',type=str, metavar='TETRA_LIBRARY_FILE')
  parser.add_argument('tetra_output_base',type=str, metavar='TETRA_OUTPUT_DIR')
  parser.add_argument('translation_file',type=str, metavar='MINERAL_FRACTIONS_FILE')
  parser.add_argument('output', type=str, metavar='OUTPUT')
  args = parser.parse_args()

  lib         = envi.open(args.tetra_library_file+'.hdr', args.tetra_library_file)
  lib_rfl     = lib.spectra.copy()
  lib_records = [int(q) for q in lib.metadata['record']]
  names       = [q.strip() for q in lib.metadata['spectra names']]
  wl          = np.array([float(q) for q in lib.metadata['wavelength']])

  nminerals = 10
  out_data = None

  # Find the mapping of output file names to EMIT mineral fractions
  # The 12 columns are: record number, then the name, then 10 fractions
  with open(args.translation_file,'r') as fin:
    lines = fin.readlines()
  emit_band_names = lines[0].split(',')[2:]
  chanmap = {}
  for line in lines[1:]:
      toks = line.split(',')
      tetra_spec_record = int(toks[0].strip())
      record = tetra_spec_record
      fracs = s.array([float(q.strip()) for q in toks[2:]])
      chanmap[record] = fracs

  # read expert system file and strip comments
  with open(args.tetra_expert_file,'r') as fin:
    lines_commented = fin.readlines()
  lines, orig_lineno = [],[]
  for lineno, line in enumerate(lines_commented):
    if not line.strip().startswith('\#'):
      orig_lineno.append(lineno)
      lines.append(line)

  # Go through expert system file one line at a time
  i, group, spectrum, output_data, header, out_hdr, rows, cols = \
       0, None, None, None, True, None, 0,0
  while i<len(lines):

      # The Header flag excludes the definitions at the start
      if lines[i].startswith('BEGIN SETUP'):
          header = False
      elif header:
          i = i + 1
          continue

      if lines[i].startswith('group'):
          group = int(lines[i].strip().split()[1])

      if lines[i].startswith('endaction'):
          if group in [1,2]:

              # get band depth of spectrum library file
              try: 
                  bd_library = band_depth(wl, rfl, features[0])
              except ValueError:
                  i = i + 1
                  continue
              
              # get band depths from map.  
              groupdir = args.tetra_output_base+'/group.'+str(group)+'um/'
              hdrpath = groupdir+filename+'.depth.gz.hdr'
              datapath = groupdir+filename+'.depth.gz'
              print('loading',hdrpath)
              try:
                  # read header, arrange output data files
                  hdr = envi.read_envi_header(hdrpath)
                  offs = int(hdr['header offset'])
                  if out_data is None:
                      out_hdr = hdr.copy()
                      out_hdr['interleave'] = 'bil'
                      out_hdr['data type'] = 4
                      out_hdr['wavelengths'] = '{'+','.join([str(q) for q in wl])+'}'
                      out_hdr['bands'] = nminerals
                      out_hdr['header offset'] = 0
                      out_hdr['band names'] = emit_band_names
                      cols = int(hdr['samples']) 
                      rows = int(hdr['lines']) 
                      out_data = np.zeros((rows, cols, nminerals),
                                   dtype=np.float32)
                  
                  # read band depth
                  with open(datapath,'rb') as fin:
                    compressed = fin.read()
                  decompressed = gzip.decompress(compressed)
                  bd = np.frombuffer(decompressed, dtype=np.uint8,
                      count=(rows*cols)+offs)
                  bd = bd[offs:] # one line offset by convention? 
                  nz = np.where(bd!=0)[0]
                  bd = bd.reshape((rows,cols))
                 
                  # normalize to the depth of the library spectrum,
                  # translating to aerial fractions
                  bd = bd.astype(dtype=np.float32) / 255.0 * bd_scaling
                  bd_map = bd / bd_library 

                  bd_map[np.logical_not(np.isfinite(bd_map))] = 0
                  bd_map[bd_map<0] = 0
                  bd_map[bd_map>1] = 1
                                    
                  # determine the mix of EMIT minerals
                  my_chans = chanmap[record].copy()
                  my_chans = my_chans.reshape(1,1,nminerals)
                  out_data = out_data + bd_map.reshape((rows,cols,1)) @ my_chans

              except FileNotFoundError:
                 print(filename+'.depth.gz.hdr not found')

      if 'SMALL' in lines[i]:
          record = int(lines[i].strip().split()[3])
          rfl = lib_rfl[lib_records.index(record),:]

      if 'define output' in lines[i]:
          filename = lines[i+2].strip().split()[0]
          bd_scaling = float(lines[i+3].strip().split()[4])

      if 'define features' in lines[i]:
          j= i+1
          features = []
          while ('endfeatures' not in lines[j]):
              toks = lines[j].strip().split()
              if len(toks)>5 and toks[0].startswith('f') and toks[1] == 'DLw':
                  features.append([float(f) for f in toks[2:6]])
              j = j + 1
      i = i + 1
    
  # write as BIL interleave
  out_data = np.transpose(out_data,(0,2,1))
  with open(args.output,'wb') as fout:
    fout.write(out_data.astype(dtype=np.float32).tobytes())
  envi.write_envi_header(args.output+'.hdr', out_hdr)
 
if __name__ == "__main__":
  main()



