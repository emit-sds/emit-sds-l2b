// David R Thompson
//

#include <iostream>
#include <fstream>
#include <list>
#include "SpectralData.h"
#include "USGS_SpectralDataReader.h"
#define MAX_SPECTRA 10000
#define MAX_CHN 500
#define MAX_STR 40

int main(int argc, char **argv) 
{ 
    USGS_SpectralDataReader reader(argv[1]);
    reader.cLoadRecords();
    int numrecords = reader.aGetNumRecords();
    int nchan=0, tocopy,k,i,j=0,nspec=-1,state=0;
    float *spectra = (float *) malloc(sizeof(float)*MAX_CHN * MAX_SPECTRA);
    char *names = (char *) malloc(sizeof(char) * MAX_STR * MAX_SPECTRA);
    int *records = (int *) malloc(sizeof(int) * MAX_SPECTRA);
    float *wavelengths =  (float *) malloc(sizeof(float)*MAX_CHN);
    wavelengths[0] = -1;

    // Read one record at a time
    for (i=0; i<numrecords; i++) {
             
      USGS_SpectralDataReader::USGSRecord record = reader.aGetRecordUSGSFormat(i);

      // Fill in nuumber of channels
      if (record.itchan > 0) {
        nchan = record.itchan;
      }

      // New spectrum? Record its name and record number
      if ((int(record.icflag[0]) == 0) && 
          (int(record.icflag[1]) == 0) && 
          (strncmp((const char *) record.ititl,"..",2)!=0)){

        nspec++;
        j = 0;
        records[nspec] = i;
        strncpy((char *) &(names[MAX_STR*nspec]), (const char *) record.ititl, 40);
        names[MAX_STR*nspec+39] = '\0';
      } 

      // Copy spectrum data.  Test to see whether this is the first 
      // record (copy up to 256 values) or a continuation record (copy
      // up to 383 values).
      if ((record.icflag[0] == 0) && (record.icflag[1] == 0)){
        tocopy = min(nchan, 256);
        memcpy(&(spectra[MAX_CHN * nspec]), record.data, tocopy*sizeof(float)); 
        j=tocopy;
      }
      else if ((record.icflag[0] == 1) && (record.icflag[1] == 0)){
        tocopy = min(nchan-j, 383);
        memcpy(&(spectra[MAX_CHN * nspec + j]), record.cdata, tocopy*sizeof(float)); 
        j = j + tocopy;
      }

      // The first spectrum shouuld carry wavelength values
      if ((j == nchan) && (wavelengths[0] == -1))
      {
        memcpy(wavelengths, &(spectra[MAX_CHN * nspec]), nchan*sizeof(float)); 
      }
    }

    // Write binary data 
    ofstream out(argv[2], ios::out);
    for (i=0; i<nspec; i++){
      out.write((char *) &(spectra[MAX_CHN*i]), sizeof(float)*nchan);
    }
    out.close();

    // Write the header
    string hdrfile = string(argv[2]).append(".hdr");
    ofstream header(hdrfile, ios::out);
    header << "ENVI" << endl;
    header << "bands = " << nchan << endl;
    header << "samples = 1" << endl;
    header << "lines = " << nspec << endl;
    header << "wavelength = {";
    for (i=0; i<(nchan-1); i++) {
      header << wavelengths[i] << ",";
    }
    header << wavelengths[nchan-1] << "}" << endl;
    header << "record = {";
    for (i=0; i<(nspec-1); i++) {
      header << records[i] << ",";
    }
    header << records[nspec-1] << "}" << endl;
    header << "header offset = 0 " << endl << "data type = 4" << endl;
    header << "interleave = bip " << endl << "byte order = 0" << endl;
    header.close();

    // Clean up
    free(spectra);
    free(records);
    free(names);
    free(wavelengths);
    return 0; 
} 

