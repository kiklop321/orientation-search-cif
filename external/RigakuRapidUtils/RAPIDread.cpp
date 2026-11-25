/******************************************************************************
 *                                 RAPIDread.cpp                              *
 *                                                                            *
 * A Matlab MEX function for reading x-ray diffraction image plate data       *
 * from Rigaku Rapid II.                                                      *
 *                                                                            *
 * The code heavy utilizes a Rigaku code for reading raw binary image data.   *
 * See the original Rigaku code: rigaku.h, rigaku.cpp,                        *
 * and license-rigaku.txt file enclosed                                       *
 *                                                                            *
 * This is a MEX-file for MATLAB.                                             *
 * Author: Zdenek Matej, Apr 2011                                             *
 *                                                                            *
 * RAPIDread code is distributed under the Simplified BSD License.            *
 *                                                                            *
 ******************************************************************************/

/*-----------------------------------------------------------------------------

 Copyright (c) 2011, Zdenek Matej, Charles University in Prague
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  ---------------------------------------------------------------------------*/


/*=============================================================================
   
   Calling:
             xd = RAPIDread(filename)

   This should read x-ray data from a file 'filename' and return them in
   a structure 'xd'. The field 'xd.image' should contain the image intensity
   data as 4 byte unsigned integers (uint32).

   Note: Using  uint32 storage is quite memory wasting, because the image data
   have range 2^15 * 32 = 2^20 (15 bits of data + 1 bit as a sign of Hi/Lo data,
   Lo data are multiplied by Hi/Lo ratio (32), this gives a range 2^20).
   Hopefully this type of storage should be straightforward for further treatment.
   You can try >> image(dx.image); axis image;

 % ----------------------------------------------------------------------------
 
   To COMPILE:

     [1] Mex file for Matlab:
                    >> mex RAPIDread.cpp rigaku.cpp

     [2] Original Rigaku program for Linux:
                    $ g++ -o ReadRAXISImage -D__RIGAKU_MAIN__ rigaku.cpp

% -----------------------------------------------------------------------------
 
   ============================================================================*/


#include "mex.h"
#include "matrix.h"

#include <sstream>

#include "rigaku.h" // include Rigaku code

// Truncate empty characters at the end of the string
std::string strtrunc(const std::string in);

void mexFunction(
		 int          nlhs,
		 mxArray      *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
  char *filename = NULL;

  /* Check for proper number of arguments */
  if (nrhs != 1)
    mexErrMsgTxt("One input required.");
  else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");

  /* input must be a string */
  if ( mxIsChar(prhs[0]) != 1)
    mexErrMsgTxt("Input must be a string.");

  /* copy the string data from prhs[0] into a C string filename. */
  filename = mxArrayToString(prhs[0]);
  
  if ( filename == NULL ) 
    mexErrMsgTxt("Could not convert input to string.");

  // Open the file

  std::ifstream inp(filename,std::ios::in|std::ios::binary);
  if ( inp.fail() || inp.bad() ) {
    if(inp.is_open())
      inp.close();
    mexErrMsgTxt("Could not open file. No data read.");
    return;
  }

  // Read the header of the image
 
  RigakuHeader header;
  if ( !ReadHeader(inp,header) ) {
    if(inp.is_open())
      inp.close();
    mexErrMsgTxt("Could not read data header. No data read.");
    return;
  }

  std::ostringstream sout;
  OutputHeaderInfo(header,sout);

  /* create output structure */
  const char* fieldnames[12] = {"filename","sample","username","comment","date","image",
                                "time","zero","tube","Lambda","collimator","headerInfo"};
  mxArray *out_struct = mxCreateStructMatrix(1,1,12,fieldnames);
  
  /* set filename field */
  mxSetField( out_struct, 0, "filename", mxCreateString(filename) );

  /* set header information text */
  mxSetField( out_struct, 0, "headerInfo", mxCreateString(sout.str().c_str()) );
  
  // Clear the sout stream
  sout.str("");
  
  /* set sample field */
  mxSetField( out_struct, 0, "sample",
              mxCreateString( strtrunc(std::string(header.crnm,20)).c_str() ) );
  
  /* set operator name field */
  mxSetField( out_struct, 0, "username",
              mxCreateString( strtrunc(std::string(header.opnm,20)).c_str() ) );
  
  /* set comment field */
  mxSetField( out_struct, 0, "comment",
              mxCreateString( strtrunc(std::string(header.memo,80)).c_str() ) );
  
  /* set date field */
  mxSetField( out_struct, 0, "date",
              mxCreateString( strtrunc( std::string(header.date,12) ).c_str() ) );
  
  /* set time field */
  mxSetField( out_struct, 0, "time", mxCreateDoubleScalar( 60 * header.fext ) );
  
  /* set zero field */
  mxSetField( out_struct, 0, "zero", mxCreateNumericMatrix( 1, 2, mxINT32_CLASS, mxREAL ) );
  
  *((int*) mxGetPr( mxGetField( out_struct, 0, "zero" ) ) + 0) = header.drtz;
  *((int*) mxGetPr( mxGetField( out_struct, 0, "zero" ) ) + 1) = header.drtx;
  
  /* set tube field */
  mxSetField( out_struct, 0, "tube",
              mxCreateString( strtrunc(std::string(header.trgt,4)).c_str() ) );
  
  /* set time field */
  mxSetField( out_struct, 0, "Lambda", mxCreateDoubleScalar( header.wlng ) );
  
  /* set collimator field */
  mxSetField( out_struct, 0, "collimator",
              mxCreateString( strtrunc(std::string(header.colm,20)).c_str() ) );
  
  /* set image field */
  
  // Allocate space to store the entire image

  LONG size = header.rlng*header.rnum;
  unsigned short *pBuffer = new unsigned short [size];

  // Read in the image
  
  if ( !ReadImageBuffer( inp, header, pBuffer) ) {
    if(inp.is_open())
      inp.close();
    delete [] pBuffer;
    mexErrMsgTxt("Could not read image data.");
    return;
  }
  
  /* allocate memory for image data */
  mxSetField( out_struct, 0, "image",
	      mxCreateNumericMatrix( header.zpxl, header.xpxl, mxINT32_CLASS, mxREAL ) );

  /* read image data */
  unsigned short *p1 = pBuffer;
  int *p2 = (int*) mxGetPr( mxGetField( out_struct, 0, "image" ) );

  std::cout<<std::flush;

  for(LONG iz=0; iz<header.rnum; iz++)
    for(LONG ix=0; ix<header.xpxl; ix++) {
      *(p2 + ((ix+1)*header.rnum-(iz+header.ipst+1))) = UnPackPixelValue(header, *p1);
      p1++;
    }

  // Close the file

  if(inp.is_open())
    inp.close();

  /* set C-style string sout.str().c_str() to MATLAB mexFunction output*/
  plhs[0] = out_struct;

  delete [] pBuffer;
  mxFree(filename);
  return;
}

// -----------------------------------------------------------------------------
std::string strtrunc(const std::string in)
{ 
  // find the first nonempty character form the end
  int iend = in.find_last_not_of(" \0");
  
  if (iend==std::string::npos)
    iend = 0;
  else
    iend++;
      
  return std::string( in, 0, iend );
}
