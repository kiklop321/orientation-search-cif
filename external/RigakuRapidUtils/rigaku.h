//
// Copyright © 2001-2007 Rigaku Americas Corp.
//                       9009 New Trails Drive
//                       The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Rigaku Americas
//
// All rights reserved
//
// ReadRAXISImage.cpp    Initial author: T.L.Hendrixson           Aug 2001
//
//
// Description:
//
// This program reads in an R-AXIS image written using the Rigaku-style
// format and writes the values of a few of the pixels to standard output.
// It additionally writes out some of the information contained in the
// header, such as axis positions, exposure time, pixel information.
//
// This program is primarily intended for illustrative purposes, and so
// there have been no attempts made to optimize either the reads or the
// byte swapping used to get the data in the correct endianness for the
// computer being used.
//
// Although the routines used for reading the image suggest that an RAXIS
// class should be developed, that is left as an exercise for the reader.
//
// This program assumes that the image was written using either SGI Control
// software or the instrument server provided with CrystalClear.  In these
// programs, the data is written in the file so that it is big-endian.
//
//
// Calling syntax:
//
// ReadRAXISImage.exe filename
//
// where filename is the name of the image file to be read.
//
// LEGAL ISSUES
// ============
//
// In plain English:
//
// 1. We don't promise that this software works.  (But if you find any bugs,
//    please let us know!)
//
// 2. You can use this software for whatever you want.  You don't have to
//    pay us.
//
// 3. You may not pretend that you wrote this software.  If you use it in a
//    program, you must acknowledge somewhere in your documentation that
//    you've used Rigaku code.
//
// In legalese:
//
// Rigaku makes NO WARRANTY or representation, either express or implied,
// with respect to this software, its quality, accuracy, merchantability, or
// fitness for a particular purpose.  This software is provided "AS IS", and
// you, its user, assume the entire risk as to its quality and accuracy.
//
// Permission is hereby granted to use, copy, modify, and distribute this
// software (or portions thereof) for any purpose, without fee, subject to
// these conditions:
//
// 1. If any part of the source code for this software is distributed, then
//    this notice must be included, with this copyright and no-warranty
//    notice unaltered; and any additions, deletions, or changes to the
//    original files must be clearly indicated in accompanying documentation.
//
// 2. If only executable code is distributed, then the accompanying
//    documentation must state that "this software is based in part on
//    software developed at Rigaku".
//
// 3. Permission for use of this software is granted only if the user accepts
//    full responsibility for any undesirable consequences; Rigaku accepts
//    NO LIABILITY for damages of any kind.
//

#ifndef __RIGAKU_RAPID_17042011_H__

/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/

#include <string>
#include <iostream>
#include <fstream>
//#include <sstream>

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

// long should have 4 bytes (under 64-bit Linux long has 8 bytes)
#define  LONG  int

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

struct RigakuHeader {
  char  name[10];     /* type of instrument */
  char  vers[10];     /* version */
  char  crnm[20];     /* crystal name */
  char  crsy[12];     /* crystal system */
  float alng;         /* a */
  float blng;         /* b */
  float clng;         /* c */
  float alfa;         /* alpha */
  float beta;         /* beta */
  float gamm;         /* gamma */
  char  spgr[12];     /* space group symbol */
  float mosc;         /* mosaic spread */
  char  memo[80];     /* memo, comments */
  char  res1[84];     /* reserved space for future use */

  char  date[12];     /* date of measurement */
  char  opnm[20];     /* username of account collecting image */
  char  trgt[4];      /* type of X-ray target (Cu, Mo, etc.) */
  float wlng;         /* X-ray wavelength */
  char  mcro[20];     /* type of monochromator */
  float m2ta;         /* monochromator 2theta (deg) */
  char  colm[20];     /* collimator size and type */
  char  filt[4];      /* filter type (Ni, etc.) */
  float camr;         /* crystal-to-detector distance */
  float vltg;         /* generator voltage (kV) */
  float crnt;         /* generator current (mA) */
  char  focs[12];     /* focus info */
  char  optc[80];     /* xray memo */
  LONG  cyld;         /* IP shape, 0=flat,1=cylinder */
  float weis;         /* Weissenberg oscillation 1 */
  char  res2[56];     /* reserved space for future use */

  char  mnax[4];      /* crystal mount axis closest to spindle axis */
  char  bmax[4];      /* crystal mount axis closest to beam axis */
  float phi0;         /* datum phi angle (deg) */
  float phis;         /* phi oscillation start angle (deg) */
  float phie;         /* phi oscillation end angle (deg) */
  LONG  oscn;         /* frame number */
  float fext;         /* exposure time (min) */
  float drtx;         /* direct beam X position */
  float drtz;         /* direct beam Z position */
  float omga;         /* goniostat angle omega */
  float fkai;         /* goniostat angle chi */
  float thta;         /* goniostat angle 2theta */
  float mu;           /* spindle inclination angle */
  char  scan[204];    /* reserved space for future use */
                      /* This space is now used for storing the scan
                         template information - tlh, 01 Feb 1999 */

  LONG  xpxl;         /* number of pixels in X direction */
  LONG  zpxl;         /* number of pixels in Z direction */
  float xsiz;         /* size of pixel in X direction (mm) */
  float zsiz;         /* size of pixel in Z direction (mm) */
  LONG  rlng;         /* record length (bytes) */
  LONG  rnum;         /* number of records (lines) in image */
  LONG  ipst;         /* starting line number */
  LONG  ipnm;         /* IP number */
  float rato;         /* photomultiplier output hi/lo ratio */
  float ft_1;         /* fading time, end of exposure to start of read */
  float ft_2;         /* fading time, end of exposure to end of read */
  char  host[10];     /* type of computer (IRIS, VAX) => endian */
  char  ip[10];       /* type of IP */
  LONG  dr_x;         /* horizontal scanning code: 0=left->right, 1=>right->left */
  LONG  dr_z;         /* vertical scanning code: 0=down->up, 1=up->down */
  LONG  drxz;         /* front/back scanning code: 0=front, 1=back */
  float shft;         /* pixel shift, R-AXIS V */
  float ineo;         /* intensity ratio E/O R-AXIS V */
  LONG  majc;         /* magic number to indicate next values are legit */
  LONG  naxs;         /* Number of goniometer axes */
  float gvec[5][3];   /* Goniometer axis vectors */
  float gst[5];       /* Start angles for each of 5 axes */
  float gend[5];      /* End angles for each of 5 axes */
  float goff[5];      /* Offset values for each of 5 axes */
  LONG  saxs;         /* Which axis is the scan axis? */
  char  gnom[40];     /* Names of the axes (space or comma separated?) */

/*
 * Most of below is program dependent.  Different programs use
 * this part of the header for different things.  So it is essentially
 * a big "common block" area for dumping transient information.
 */
  char  file[16];     /* */
  char  cmnt[20];     /* */
  char  smpl[20];     /* */
  LONG  iext;         /* */
  LONG  reso;         /* */
  LONG  save;         /* */
  LONG  dint;         /* */
  LONG  byte;         /* */
  LONG  init;         /* */
  LONG  ipus;         /* */
  LONG  dexp;         /* */
  LONG  expn;         /* */
  LONG  posx[20];     /* */
  LONG  posy[20];     /* */
  int   xray;         /* */
  char  res5[768];    /* reserved space for future use */
};

struct axis {

  axis() : start(0.0),end(0.0),offset(0.0),name("") {}

  double start,end,offset;
  std::string name;
};

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

extern const int MinimumRigakuHeaderLength;

extern std::string ERROR_OPENING_FILE;
extern std::string ERROR_READING_HEADER;
extern std::string ERROR_READING_IMAGE;

/****************************************************************************
 *                           Function prototypes                            *
 ****************************************************************************/

  // This routine returns a string containing a description of the program
  // usage.
std::string ProgramUsage(void);

  // This routine reads a Rigaku-style header from the input stream
  // <VAR>InputFile</VAR>, and returns a boolean indicating if the header
  // could be read successfully.  If the value returned is <I>true</I>,
  // then the header was read successfully and the header information will
  // be placed in the header structure <VAR>Header</VAR>.  If the value
  // returned is <I>false</I>, then an error occurred while attempting to
  // read the header and the contents of <VAR>Header</VAR> are undefined.
bool ReadHeader(std::ifstream &InputFile, RigakuHeader &Header);

  // This routine reads a Rigaku-style image from the input stream
  // <VAR>InputFile</VAR>, and returns a boolean indicating if the image
  // could be read successfully.  The routine assumes that the header
  // structure <VAR>Header</VAR> contains valid header information for
  // the image being read.  If the value returned is <I>true</I>,
  // then the image was read successfully and the image will be
  // placed as (possibly packed) unsigned shorts at the memory location
  // indicated by pImageBuffer.  If the value returned is <I>false</I>,
  // then an error occurred while attempting to read the image and the
  // contents of <VAR>pImageBuffer</VAR> are undefined.  It is up to the
  // calling routine to ensure that <VAR>pImageBuffer</VAR> contains
  // enough space to hold the image.
bool ReadImageBuffer(std::ifstream &InputFile,
                    const RigakuHeader &Header,
                    unsigned short * const pImageBuffer);

  // This routine swaps the unsigned short <VAR>Value</VAR> so that the
  // endianness of <VAR>Value</VAR> is reversed and then returns the
  // swapped value.  This routine assumes that the size of a unsigned short
  // is 2 bytes.
unsigned short Swap(unsigned short Value);

  // This routine swaps the long <VAR>Value</VAR> so that the endianness
  // of <VAR>Value</VAR> is reversed and then returns the swapped value.
  // This routine assumes that the size of a long is 4 bytes.
LONG Swap(LONG Value);

  // This routine swaps the float <VAR>Value</VAR> so that the endianness
  // of <VAR>Value</VAR> is reversed and then returns the swapped value.
  // This routine assumes that the size of a float is 4 bytes and that the
  // IEEE standard format for floats is used.
float Swap(float Value);

  // This routine swaps the values in a Rigaku-style header,
  // <VAR>Header</VAR> so that the endianness of the <CODE>long</CODE>
  // and </CODE>float</CODE> is reversed.
void Swap(RigakuHeader &Header);

  // This routine returns a long containing the actual value of a packed
  // R-AXIS image pixel, <VAR>PixelValue</VAR>.  The Rigaku-style header
  // contains information concerning the packing scheme used.
LONG UnPackPixelValue(const RigakuHeader &Header,
                     const unsigned short PixelValue);

  // This routine writes out pixel values for some prefined pixel positions
  // in the image.  The values are written to standard output.  The header
  // of the image is contained in <VAR>Header</VAR> and
  // <VAR>pImageBuffer</VAR> points to the first pixel in the image.
void OutputSamplePixelValues(const RigakuHeader &Header,
                            unsigned short * const pImageBuffer);

  // This routine writes out some basic information contained in the
  // header <VAR>Header</VAR>.
void OutputHeaderInfo(const RigakuHeader &Header, std::ostream &sout = std::cout);

  // This routine writes out the axis information, <VAR>Axis</VAR>, on
  // the output stream <VAR>Stream</VAR>
std::ostream &operator<<(std::ostream &Stream, axis &Axis);

/****************************************************************************
 ****************************************************************************/

#endif /* __RIGAKU_RAPID_17042011_H__ */
