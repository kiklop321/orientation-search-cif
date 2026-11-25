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

/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/

#ifdef WIN32
#pragma warning(disable:4786) // so warning about debug term > 255 chars ignored
#endif

#include <sstream>

#include "rigaku.h"

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

const int MinimumRigakuHeaderLength = 1400;

std::string ERROR_OPENING_FILE   = "ERROR Can not open file: ";
std::string ERROR_READING_HEADER = "ERROR Can not read Rigaku-style header.";
std::string ERROR_READING_IMAGE  = "ERROR Can not read image data.";

/****************************************************************************
 ****************************************************************************/

// If __RIGAKU_MAIN __ is defined, create an entry point main function
#ifdef __RIGAKU_MAIN__
int main(int argc,
        char *argv[])
{
  // Make sure that there is a filename on the command line.

  if(2 != argc){
     std::cerr << ProgramUsage() << std::endl;
     return 1;
  }

  // Open the file

  std::ifstream inp(argv[1],std::ios::in|std::ios::binary);
  if(inp.fail() || inp.bad()){
     std::cerr << ERROR_OPENING_FILE << argv[1] << std::endl;
     if(inp.is_open())
        inp.close();
     return 2;
  }

  // Read the header of the image

  RigakuHeader header;
  if(!ReadHeader(inp,header)){
     std::cerr << ERROR_READING_HEADER << std::endl;
     if(inp.is_open())
        inp.close();
     return 3;
  }

  OutputHeaderInfo(header);

  // Allocate space to store the entire image

  LONG size = header.rlng*header.rnum;
  unsigned short *pBuffer = new unsigned short [size/2];

  // Read in the image

  if(!ReadImageBuffer(inp,header ,pBuffer)){
     std::cerr << ERROR_READING_IMAGE << std::endl;
     if(inp.is_open())
        inp.close();
     delete [] pBuffer;
     return 4;
  }

  if(inp.is_open())
     inp.close();

  OutputSamplePixelValues(header ,pBuffer);

  delete [] pBuffer;

  return 0;
}
#endif /* __RIGAKU_MAIN __ */

/****************************************************************************
 ****************************************************************************/
std::string ProgramUsage(void)
{
  std::string s = "ReadRAXISImage.exe image_filename\n\n";
  s += "image_filename : Filename of image to be read.\n";

  return s;
}

/****************************************************************************
 ****************************************************************************/
bool ReadHeader(std::ifstream &InputFile,
               RigakuHeader &Header)
{
  // Read in the minimal amount of header information.

  InputFile.read(reinterpret_cast<char *>(&Header),MinimumRigakuHeaderLength);

  if(InputFile.bad() || InputFile.fail())
     return false;

  // Check to see if endianness of the header matches endianness of this
  // architecture.  If it doesn't, then we have to swap the non-character
  // elements of the header.

  unsigned short one = 1;
  char *p = reinterpret_cast<char *>(&one);
  if(1 == static_cast<int>(*p))
     Swap(Header);

  // Read the additional padding in the header

  LONG size = Header.rlng-MinimumRigakuHeaderLength;
  p = new char [size];
  InputFile.read(p,size);
  delete [] p;

  return !(InputFile.bad() || InputFile.fail());
}

/****************************************************************************
 ****************************************************************************/
bool ReadImageBuffer(std::ifstream &InputFile,
                    const RigakuHeader &Header,
                    unsigned short * const pImageBuffer)
{
  // Figure out how much of the image is contained in the file.  Note that
  // if the image file contains only a partial image, we will NOT construct
  // an entire image.  The header element ipst should be used to determine
  // the line in the entire image corresponding to a line in the image file.
  //
  // Note that rlng is the size of a record in bytes, and we are assuming
  // 2 bytes/pixel.

  LONG size = Header.rlng*Header.rnum/2;

  // Read in the data.

  // modified 17.4.2011 - But we are reading 2*size bytes
  InputFile.read(reinterpret_cast<char *>(pImageBuffer),2*size);

  if(InputFile.bad() || InputFile.fail())
     return false;

  // Do we need to swap the data?

  unsigned short one = 1;
  char *p = reinterpret_cast<char *>(&one);

  if(1 == static_cast<int>(*p)){

     // Swap each pixel value in the image buffer

     unsigned short *p = const_cast<unsigned short *>(pImageBuffer);
     unsigned short *pEnd = p+size;
     for(NULL; p < pEnd; p++)
        *p = Swap(*p);

  }

  return true;

}

/****************************************************************************
 ****************************************************************************/
unsigned short Swap(const unsigned short Value)
{
  unsigned short swap = Value;
  char *p = reinterpret_cast<char *>(&swap);

  char tmp = p[0];
  p[0] = p[1];
  p[1] = tmp;

  return swap;
}

/****************************************************************************
 ****************************************************************************/
LONG Swap(const LONG Value)
{
  LONG swap = Value;
  char *p = reinterpret_cast<char *>(&swap);

  char tmp = p[0];
  p[0] = p[3];
  p[3] = tmp;

  tmp = p[1];
  p[1] = p[2];
  p[2] = tmp;

  return swap;
}

/****************************************************************************
 ****************************************************************************/
float Swap(const float Value)
{
  float swap = Value;
  char *p = reinterpret_cast<char *>(&swap);

  char tmp = p[0];
  p[0] = p[3];
  p[3] = tmp;

  tmp = p[1];
  p[1] = p[2];
  p[2] = tmp;

  return swap;
}

/****************************************************************************
 ****************************************************************************/
void Swap(RigakuHeader &Header)
{

// Swap the bytes in the long elements of the header.

  Header.cyld = Swap(Header.cyld);
  Header.oscn = Swap(Header.oscn);
  Header.xpxl = Swap(Header.xpxl);
  Header.zpxl = Swap(Header.zpxl);
  Header.rlng = Swap(Header.rlng);
  Header.rnum = Swap(Header.rnum);
  Header.ipst = Swap(Header.ipst);
  Header.ipnm = Swap(Header.ipnm);
  Header.dr_x = Swap(Header.dr_x);
  Header.dr_z = Swap(Header.dr_z);
  Header.drxz = Swap(Header.drxz);

// Swap the bytes in the float elements of the header.

  Header.alng = Swap(Header.alng);
  Header.blng = Swap(Header.blng);
  Header.clng = Swap(Header.clng);
  Header.alfa = Swap(Header.alfa);
  Header.beta = Swap(Header.beta);
  Header.gamm = Swap(Header.gamm);
  Header.mosc = Swap(Header.mosc);
  Header.wlng = Swap(Header.wlng);
  Header.m2ta = Swap(Header.m2ta);
  Header.camr = Swap(Header.camr);
  Header.vltg = Swap(Header.vltg);
  Header.crnt = Swap(Header.crnt);
  Header.weis = Swap(Header.weis);
  Header.phi0 = Swap(Header.phi0);
  Header.phis = Swap(Header.phis);
  Header.phie = Swap(Header.phie);
  Header.fext = Swap(Header.fext);
  Header.drtx = Swap(Header.drtx);
  Header.drtz = Swap(Header.drtz);
  Header.omga = Swap(Header.omga);
  Header.fkai = Swap(Header.fkai);
  Header.thta = Swap(Header.thta);
  Header.mu   = Swap(Header.mu);
  Header.xsiz = Swap(Header.xsiz);
  Header.zsiz = Swap(Header.zsiz);
  Header.rato = Swap(Header.rato);
  Header.ft_1 = Swap(Header.ft_1);
  Header.ft_2 = Swap(Header.ft_2);

// New variables for RAPID image headers

  Header.shft = Swap(Header.shft);
  Header.ineo = Swap(Header.ineo);
  Header.majc = Swap(Header.majc);
  Header.naxs = Swap(Header.naxs);
  for(int i = 0; i < 5; i++){
     for(int j = 0; j < 3; j++)
        Header.gvec[i][j] = Swap(Header.gvec[i][j]);
     Header.gst[i] = Swap(Header.gst[i]);
     Header.gend[i] = Swap(Header.gend[i]);
     Header.goff[i] = Swap(Header.goff[i]);
  }
  Header.saxs = Swap(Header.saxs);

  return;
}

/****************************************************************************
 ****************************************************************************/
LONG UnPackPixelValue(const RigakuHeader &Header,
                     const unsigned short PixelValue)
{
  // Pixels are packed so that if the most significant bit is set, then
  // the value has been bit-shifted by the amount specified in the rato
  // header element.
  //
  // Note that this form of packing is actually quite inspired as the
  // amount of data lost due to packing is much less than the error in
  // the measurement due to normal counting statistics (assuming that
  // they are Poissonian).

  if(0x7fff >= PixelValue)
     return PixelValue;

  // We only need to use the lower 15 bits of the value to uncompress it.

  LONG Value = 0x7fff&PixelValue;

  // Shift the value by the ratio in the header.

  Value *= Header.rato;

  return Value;
}

/****************************************************************************
 ****************************************************************************/
void OutputSamplePixelValues(const RigakuHeader &Header,
                            unsigned short * const pImageBuffer)
{
  // Output a few pixel values.  Note that the first line in the image
  // buffer may not correspond to the first line in a "complete" image.

  std::cout << std::endl << "Sample pixel information: " << std::endl;

  std::cout << "   First pixel (0," << Header.ipst << ") in image = " <<
               UnPackPixelValue(Header, *pImageBuffer) << std::endl;

  LONG index = (Header.rnum/2-1)*Header.xpxl+(Header.xpxl/2-1);
  std::cout << "   Middle pixel (" << Header.xpxl/2 << ',' <<
                Header.rnum/2+Header.ipst << ") in image = " <<
                UnPackPixelValue(Header ,pImageBuffer[index]) << std::endl;

  index = Header.rnum*Header.xpxl-1;
  std::cout << "   Last pixel (" << Header.xpxl << ',' <<
                Header.rnum+Header.ipst << ") in image = " <<
                UnPackPixelValue(Header ,pImageBuffer[index]) << std::endl;

  return;
}

/****************************************************************************
 ****************************************************************************/
std::string fromFORTRAN(char *string, int size)
{
  std::ostringstream s;
  for(int i = 0; i < size; i++)
     s << string[i];
  return s.str();
}

/****************************************************************************
 ****************************************************************************/
void OutputHeaderInfo(const RigakuHeader &Header, std::ostream &sout)
{
  sout << "Header information: " << std::endl;

  sout << "   InstrumentType: " << fromFORTRAN((char *)Header.name,10) << std::endl;
  sout << "   Version: " << fromFORTRAN((char *)Header.vers,10) << std::endl;
  sout << "   Crystal Name: " << fromFORTRAN((char *)Header.crnm,20) << std::endl;
  sout << "   Crystal System: " << fromFORTRAN((char *)Header.crsy,12) << std::endl;
  sout << "   a: " << Header.alng << std::endl;
  sout << "   b: " << Header.blng << std::endl;
  sout << "   c: " << Header.clng << std::endl;
  sout << "   alpha: " << Header.alfa << std::endl;
  sout << "   beta: " << Header.beta << std::endl;
  sout << "   gamma: " << Header.gamm << std::endl;
  sout << "   Space Group: " << fromFORTRAN((char *)Header.spgr,12) << std::endl;
  sout << "   Mosaicity: " << Header.mosc << std::endl;
  sout << "   Memo: " << fromFORTRAN((char *)Header.memo,80) << std::endl;

  sout << "   Date: " << fromFORTRAN((char *)Header.date,12) << std::endl;
  sout << "   User: " << fromFORTRAN((char *)Header.opnm,20) << std::endl;
  sout << "   Xray Target: " << fromFORTRAN((char *)Header.trgt,4) << std::endl;
  sout << "   Wavelength: " << Header.wlng << std::endl;
  sout << "   Monochromator: " << fromFORTRAN((char *)Header.mcro,20) << std::endl;
  sout << "   Monochromator 2theta: " << Header.m2ta << std::endl;
  sout << "   Collimator: " << fromFORTRAN((char *)Header.colm,20) << std::endl;
  sout << "   Filter: " << fromFORTRAN((char *)Header.filt,4) << std::endl;
  sout << "   Generator Voltage: " << Header.vltg << std::endl;
  sout << "   Generator Current: " << Header.crnt << std::endl;
  sout << "   Focus: " << fromFORTRAN((char *)Header.focs,12) << std::endl;
  sout << "   Xray Memo: " << fromFORTRAN((char *)Header.optc,80) << std::endl;
  sout << "   IP shape: " << ((1 == Header.cyld) ? "cylindrical" :
     "flat") << std::endl;
  sout << "   OscillationType: " << ((1 == Header.weis) ? "weissenberg" : "regular") << std::endl;

  sout << "   Crystal Mount (spindle axis): " << fromFORTRAN((char *)Header.mnax,4) << std::endl;
  sout << "   Crystal Mount (beam axis): " << fromFORTRAN((char *)Header.bmax,4) << std::endl;
  sout << "   Frame Number: " << Header.oscn << std::endl;
  sout << "   Exposure time: " << Header.fext << std::endl;
  sout << "   Direct beam position: (" << Header.drtx << ',' <<
     Header.drtz << ')' << std::endl;
  sout << "   Spindle Inclication Angle: " << Header.mu << std::endl;
  if (sout==std::cout) // modified 17.4.2011
    sout << "   Image Template: " << fromFORTRAN((char *)Header.scan,204) << std::endl;
  else
    sout << "   Image Template: " << std::endl;
  sout << "   Pixel information" << std::endl <<
     "      size: (" << Header.xsiz << ',' << Header.zsiz << ')' <<
     std::endl <<
     "      number: (" << Header.xpxl << ',' << Header.zpxl << ')' <<
     std::endl;
  sout << "   Line information" << std::endl <<
     "      record length: " << Header.rlng << std::endl <<
     "      total: " << Header.rnum << std::endl <<
     "      starting: " << Header.ipst << std::endl;
  sout << "   IP Number: " << Header.ipnm << std::endl;
  sout << "   Photomultiplier Ratio: " << Header.rato << std::endl;
  sout << "   Fade Time (to start of read): " << Header.ft_1 << std::endl;
  sout << "   Fade Time (to end of read): " << Header.ft_2 << std::endl;
  sout << "   Host Type/Endian: " << fromFORTRAN((char *)Header.host,10) << std::endl;
  sout << "   IP Type: " << fromFORTRAN((char *)Header.ip,10) << std::endl;
  sout << "   Horizontal Scan: " << ((0 == Header.dr_x) ? "left to right" : "right to left") << std::endl;
  sout << "   Vertical Scan: " << ((0 == Header.dr_z) ? "down to up" : "up to down") << std::endl;
  sout << "   Front/Back Scan: " << ((0 == Header.drxz) ? "front" : "back") << std::endl;
  sout << "   Pixel Shift (RAXIS V): " << Header.shft << std::endl;
  sout << "   Even/Odd Intensity Ratio (RAXIS V): " << Header.ineo << std::endl;


  axis Axes[5];

  std::string s = Header.name;
  bool IsRapid = std::string::npos != s.find("R-AXIS-CS");
  std::string ScanAxis = "Phi";

  if(0 != Header.majc){
     std::istringstream is(Header.gnom);
     for(int i = 0; i < Header.naxs; i++){
        is >> Axes[i].name;
        Axes[i].start = Header.gst[i];
        Axes[i].end = Header.gend[i];
        Axes[i].offset = Header.goff[i];
     }
     ScanAxis = Axes[Header.saxs].name;
  }
  else{
     Axes[0].name = "Omega";
     Axes[2].name = "Phi";
     if(!IsRapid){
        Axes[0].start = Header.omga;
        Axes[0].end = Header.omga;
        Axes[0].offset = 0.0;
        Axes[2].start = Header.phis;
        Axes[2].end = Header.phie;
        Axes[2].offset = Header.phi0;
        ScanAxis = Axes[2].name;
     }
     else{
        Axes[2].start = Header.omga;
        Axes[2].end = Header.omga;
        Axes[2].offset = 0.0;
        Axes[0].start = Header.phis;
        Axes[0].end = Header.phie;
        Axes[0].offset = Header.phi0;
        ScanAxis = Axes[0].name;
     }
     Axes[1].name = "Chi";
     Axes[1].start = Header.fkai;
     Axes[1].end = Header.fkai;
     Axes[1].offset = 0.0;
  }

  sout << "   Instrument axes (using " <<
               ((0 != Header.majc) ? "extended" : "conventional") <<
               " section of header):" << std::endl;
  sout << "      Scan axis: " << ScanAxis << std::endl;
  for(int i = 0; i < 5 && "" != Axes[i].name; i++)
     sout << "      " << Axes[i] << std::endl;

  axis theta,dist;
  theta.name = "2Theta";
  theta.start = Header.thta;
  theta.end = Header.thta;
  theta.offset = 0.0;
  dist.name = "Distance";
  dist.start = Header.camr;
  dist.end = Header.camr;
  dist.offset = 0.0;
  sout << "      " << theta << std::endl;
  sout << "      " << dist << std::endl;

  return;
}

/****************************************************************************
 ****************************************************************************/
std::ostream &operator<<(std::ostream &a,
                        axis &b)
{
  a << b.name << ": start = " << b.start << ", end = " << b.end <<
       ", offset = " << b.offset;
  return a;
}
