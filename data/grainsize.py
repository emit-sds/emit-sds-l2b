


import argparse
import numpy as np
import spectral.io.envi as envi
from scipy import interpolate
from emit_utils.file_checks import envi_header
import logging
import pickle
import os



def write_output(outdat, outfile, outuncdat, outuncfile):
    logging.info(f'Writing output file {outfile}')
    # write as BIL interleave
    outdat = np.transpose(outdat, (0, 2, 1))
    with open(outfile, 'wb') as fout:
        fout.write(outdat.astype(dtype=np.float32).tobytes())

    if outuncdat is not None:
        # write as BIL interleave
        outuncdat = np.transpose(outuncdat, (0, 2, 1))
        with open(outuncfile, 'wb') as fout:
            fout.write(outuncdat.astype(dtype=np.float32).tobytes())


def spectral_derivative(rfl):

    # Everything in between
    d = [np.array(-1.5*rfl[:,:,0]+2*rfl[:,:,1]-0.5*rfl[:,:,2])]
    d.append(np.empty(d[0].shape))
    for i in range(1, 283):
      d.append((rfl[:,:,i+1]-rfl[:,:,i-1])/2)
    d.append(np.array(1.5*rfl[:,:,-1]-2*rfl[:,:,-2]+0.5*rfl[:,:,-3]))

    d = np.stack(d,axis=-1)
    return d



def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('rfl_file', type=str, metavar='OUTPUT')
    parser.add_argument('model_file', type=str, metavar='Band Depth file.  4 bands (G1 BD, G1 Ref, G2 BD, G2 Ref)')
    parser.add_argument('output_file', type=str)
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--log_level', type=str, default='INFO')
    args = parser.parse_args()

    if os.path.isfile(args.output_file):
        print('already found output file...terminating')
        exit()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    logging.info(f'Opening {envi_header(args.rfl_file)}')
    rfl_ds = envi.open(envi_header(args.rfl_file))
    wl = np.array([float(x) for x in rfl_ds.metadata['wavelength']])
    
    output_header = rfl_ds.metadata.copy()
    for key in ['wavelength', 'fwhm']:
        if key in output_header.keys():
            del output_header[key]
    output_header['bands'] = 8
    output_header['band names'] = ['S1','S2','S3','S4','S5','TSI','Clay','Median Grainsize']
    output_header['header offset'] = 0
    output_header['data type'] = 4
    envi.write_envi_header(envi_header(args.output_file), output_header)

    logging.info(f'Loading Data...')
    rfl = rfl_ds.open_memmap(interleave='bip').copy().astype(np.float32)
    logging.info(f'...Data Loaded')

    logging.info(f'Calculate Spectral Derivative')
    rfl_d = spectral_derivative(rfl)

    valid_wl = np.array([x > 454 and not (x > 1298 and x < 1505) and not (x > 1775 and x < 1980) and x < 2301 for x in wl])

    rfl = rfl[...,valid_wl]
    rfl_d = rfl_d[...,valid_wl]

    rfl = np.append(rfl_d,rfl,axis=2)
    rfl = rfl.reshape((rfl.shape[0]*rfl.shape[1],rfl.shape[2]))
    rfl = np.nan_to_num(rfl)

    model = pickle.load(open(args.model_file, 'rb'))
    pred = model.predict(rfl)

    # normalize
    pred = pred / np.sum(pred,axis=-1)[:,np.newaxis]
    predc = np.cumsum(pred,axis=-1)

    size_classes = np.array([1500, 7500, 375, 187.5, 93.75, 30, 1])
    ifuns = [interpolate.interp1d(predc[n,:],size_classes) if rfl[n,-100] > 0 else np.array([-1]) for n in range(pred.shape[0])]
    del predc
    median_size = np.array([fun([0.5]) if fun != -1 else np.array([-1]) for fun in ifuns])

    pred = np.hstack([pred,median_size])

    # mask
    pred[rfl[...,-100] <= 0,:] = -9999
    
    pred = pred.reshape((rfl_ds.shape[0],rfl_ds.shape[1], pred.shape[-1]))

    write_output(pred, args.output_file, None, None)





if __name__ == "__main__":
    main()
