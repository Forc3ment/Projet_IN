/****************************************************************************
 *
 * This file is part of the ViSP software.
 * Copyright (C) 2005 - 2015 by Inria. All rights reserved.
 *
 * This software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * ("GPL") version 2 as published by the Free Software Foundation.
 * See the file LICENSE.txt at the root directory of this source
 * distribution for additional information about the GNU GPL.
 *
 * For using ViSP with software that can not be combined with the GNU
 * GPL, please contact Inria about acquiring a ViSP Professional
 * Edition License.
 *
 * See http://visp.inria.fr for more information.
 *
 * This software was developed at:
 * Inria Rennes - Bretagne Atlantique
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 * France
 *
 * If you have questions regarding the use of this file, please contact
 * Inria at visp@inria.fr
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Description:
 * Luminance feature.
 *
 * Authors:
 * Eric Marchand
 *
 *****************************************************************************/


#include <visp3/core/vpMatrix.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpDisplay.h>
#include <visp3/core/vpPixelMeterConversion.h>
#include <visp3/core/vpImageConvert.h>
#include <visp3/core/vpImageFilter.h>
#include <visp3/core/vpException.h>

#include <visp3/visual_features/vpFeatureBitplane.h>


/*!
  \file vpFeatureBitplane.cpp
  \brief Class that defines the image luminance visual feature

  For more details see \cite Collewet08c.
*/



/*!
  Initialize the memory space requested for vpFeatureBitplane visual feature.
*/
void
vpFeatureBitplane::init()
{
    if (flags == NULL)
      flags = new bool[nbParameters];
    for (unsigned int i = 0; i < nbParameters; i++) flags[i] = false;

    //default value Z (1 meters)
    Z = 1;

    firstTimeIn =0 ;

    nbr = nbc = 0;
}


void
vpFeatureBitplane::init(unsigned int _nbr, unsigned int _nbc, double _Z)
{
  init() ;

  nbr = _nbr ;
  nbc = _nbc ;

  if((nbr < 2*bord) || (nbc < 2*bord)){
    throw vpException(vpException::dimensionError, "border is too important compared to number of row or column.");
  }

  // number of feature = nb column x nb lines in the images
  dim_s = (nbr-2*bord)*(nbc-2*bord) ;

  s.resize(dim_s) ;
  
  if (pixInfo != NULL)
    delete [] pixInfo;

  pixInfo = new vpBitplane[dim_s] ;

  for(unsigned int i = 0; i<dim_s; i++) 
  {
    pixInfo[i].lbp = new bool[8];
    pixInfo[i].Ix = new double[8];
    pixInfo[i].Iy = new double[8];

  }
      

  
  Z = _Z ;
}

/*! 
  Default constructor that build a visual feature.
*/
vpFeatureBitplane::vpFeatureBitplane()
  : Z(1), nbr(0), nbc(0), bord(10), pixInfo(NULL), firstTimeIn(0), cam()
{
    nbParameters = 1;
    dim_s = 0 ;
    flags = NULL;

    init() ;
}

/*!
 Copy constructor.
 */
vpFeatureBitplane::vpFeatureBitplane(const vpFeatureBitplane& f)
  : vpBasicFeature(f), Z(1), nbr(0), nbc(0), bord(10), pixInfo(NULL), firstTimeIn(0), cam()
{
  *this = f;
}

/*!
 Copy operator.
 */
vpFeatureBitplane &vpFeatureBitplane::operator=(const vpFeatureBitplane& f)
{
  Z = f.Z;
  nbr = f.nbr;
  nbc = f.nbc;
  bord = f.bord;
  firstTimeIn = f.firstTimeIn;
  cam = f.cam;
  if (pixInfo)
    delete [] pixInfo;
  pixInfo = new vpBitplane[dim_s] ;
  for(unsigned int i=0; i< dim_s; i++)
    pixInfo[i] = f.pixInfo[i];
  return (*this);
}

/*! 
  Destructor that free allocated memory.
*/
vpFeatureBitplane::~vpFeatureBitplane() 
{
  if (pixInfo != NULL)
  {
    for(unsigned int i = 0; i<dim_s; i++) 
    {
      if (pixInfo[i].lbp != NULL)
        delete [] pixInfo[i].lbp;
      if (pixInfo[i].Ix != NULL)
        delete [] pixInfo[i].Ix;
      if (pixInfo[i].Iy != NULL)
        delete [] pixInfo[i].Iy;
    }
      
    delete [] pixInfo ;
  } 
}

/*!
  Set the value of \f$ Z \f$ which represents the depth in the 3D camera frame.

  \param Z_ : \f$ Z \f$ value to set.
*/
void
vpFeatureBitplane::set_Z(const double Z_)
{
    this->Z = Z_ ;
    flags[0] = true;
}


/*!
  Get the value of \f$ Z \f$ which represents the depth in the 3D camera frame.

  \return The value of \f$ Z \f$.
*/
double
vpFeatureBitplane::get_Z() const
{
    return Z ;
}


void
vpFeatureBitplane::setCameraParameters(vpCameraParameters &_cam) 
{
  cam = _cam ;
}


/*!

  Build a luminance feature directly from the image
*/

void
vpFeatureBitplane::buildFrom(vpImage<unsigned char> &I)
{
  unsigned int l = 0;
  //double Ix,Iy;

  double px = cam.get_px() ;
  double py = cam.get_py() ;

  /*this->w = I.getWidth();
  this->h = I.getHeight();*/


  if (firstTimeIn==0)
  { 
    firstTimeIn=1 ;
    l =0 ;
    for (unsigned int i=bord; i < nbr - bord ; i++)
    {
  	  //   cout << i << endl ;
  	  for (unsigned int j = bord ; j < nbc - bord; j++)
 	    {	
        double x=0,y=0;
    	  vpPixelMeterConversion::convertPoint(cam,
    		   j,i,
    		   x,y)  ;
    	   
    	  pixInfo[l].x = x;
    	  pixInfo[l].y = y;

    	  pixInfo[l].Z   = Z ;

    	  l++;
      }
	  }
  }

  l = 0 ;
  for (unsigned int i=bord; i < nbr - bord ; i++)
  {
    for (unsigned int j = bord ; j < nbc - bord; j++)
  	{
        pixInfo[l].lbp[0] = I[i-1][j-1]   <   I[i][j];
        pixInfo[l].lbp[1] = I[i-1][j]     <   I[i][j];
        pixInfo[l].lbp[2] = I[i-1][j+1]   <   I[i][j];
        pixInfo[l].lbp[3] = I[i][j-1]     <   I[i][j];
        pixInfo[l].lbp[4] = I[i][j+1]     <   I[i][j];
        pixInfo[l].lbp[5] = I[i+1][j-1]   <   I[i][j];
        pixInfo[l].lbp[6] = I[i+1][j]     <   I[i][j];
        pixInfo[l].lbp[7] = I[i+1][j+1]   <   I[i][j];
  	  
        pixInfo[l].Ix[0] = px * ((I[i][j-1]    <   I[i+1][j]) - (I[i-2][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Ix[1] = px * ((I[i][j-1]    <   I[i+1][j]) - (I[i-2][j]     <   I[i-1][j])) ;
        pixInfo[l].Ix[2] = px * ((I[i][j-1]    <   I[i+1][j]) - (I[i-2][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Ix[3] = px * ((I[i+1][j-1]  <   I[i+1][j]) - (I[i-1][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Ix[4] = px * ((I[i+1][j-1]  <   I[i+1][j]) - (I[i-1][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Ix[5] = px * ((I[i+2][j-1]  <   I[i+1][j]) - (I[i][j-1]     <   I[i-1][j])) ;
        pixInfo[l].Ix[6] = px * ((I[i+2][j-1]  <   I[i+1][j]) - (I[i][j-1]     <   I[i-1][j])) ;
        pixInfo[l].Ix[7] = px * ((I[i+2][j-1]  <   I[i+1][j]) - (I[i][j-1]     <   I[i-1][j])) ;

        pixInfo[l].Iy[0] = py * ((I[i][j-1]    <   I[i+1][j]) - (I[i-2][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Iy[1] = py * ((I[i][j-1]    <   I[i+1][j]) - (I[i-2][j]     <   I[i-1][j])) ;
        pixInfo[l].Iy[2] = py * ((I[i][j-1]    <   I[i+1][j]) - (I[i-2][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Iy[3] = py * ((I[i+1][j-1]  <   I[i+1][j]) - (I[i-1][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Iy[4] = py * ((I[i+1][j-1]  <   I[i+1][j]) - (I[i-1][j-1]   <   I[i-1][j])) ;
        pixInfo[l].Iy[5] = py * ((I[i+2][j-1]  <   I[i+1][j]) - (I[i][j-1]     <   I[i-1][j])) ;
        pixInfo[l].Iy[6] = py * ((I[i+2][j-1]  <   I[i+1][j]) - (I[i][j-1]     <   I[i-1][j])) ;
        pixInfo[l].Iy[7] = py * ((I[i+2][j-1]  <   I[i+1][j]) - (I[i][j-1]     <   I[i-1][j])) ;
  	  
  	  l++;
  	}
  }     
}



/*!

  Compute and return the interaction matrix \f$ L_I \f$. The computation is made
  thanks to the values of the luminance features \f$ I \f$
*/
void
vpFeatureBitplane::interaction(vpMatrix &L)
{  
  L.resize(dim_s,6) ;

  for(unsigned int m = 0; m< L.getRows(); m++)
  {
    double* Ix = pixInfo[m].Ix;
    double* Iy = pixInfo[m].Iy;

    double x = pixInfo[m].x ;
    double y = pixInfo[m].y ;
    double Zinv =  1 / pixInfo[m].Z;

    for (int i = 0; i < 8; ++i)
    {
      L[m][0] +=     Ix[i] * Zinv;
      L[m][1] +=     Iy[i] * Zinv;
      L[m][2] += -(x*Ix[i] + y*Iy[i]) * Zinv;
      L[m][3] += -Ix[i]*x*y - (1+y*y)*Iy[i];
      L[m][4] += (1+x*x)*Ix[i] + Iy[i]*x*y;
      L[m][5] += Iy[i]*x - Ix[i]*y;
    }

    L[m][0] /= 8;
    L[m][1] /= 8;
    L[m][2] /= 8;
    L[m][3] /= 8;
    L[m][4] /= 8;
    L[m][5] /= 8;
  }
}

/*!
  Compute and return the interaction matrix \f$ L_I \f$. The computation is made
  thanks to the values of the luminance features \f$ I \f$
*/
vpMatrix  vpFeatureBitplane::interaction(const unsigned int /* select */)
{
  /* static */ vpMatrix L  ; // warning C4640: 'L' : construction of local static object is not thread-safe
  interaction(L) ;
  return L ;
}


/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired
 
  \param s_star : Desired visual feature.
  \param e : Error between the current and the desired features.

*/
void
vpFeatureBitplane::error(const vpBasicFeature &s_star,
			  vpColVector &e)
{
  e.resize(dim_s);
  vpFeatureBitplane *s_star_lbp = dynamic_cast<vpFeatureBitplane>(&s_star)
  if (s_star_lbp!=NULL) 
  {
    int l=0;
    for (unsigned int i=bord; i < nbr - bord ; i++)
    {
      for (unsigned int j = bord ; j < nbc - bord; j++)
      { 
        e[l] =  pixInfo[l].lbp[0] ^ s_star.pixInfo[l].lbp[0];
        e[l] += pixInfo[l].lbp[1] ^ s_star.pixInfo[l].lbp[1];
        e[l] += pixInfo[l].lbp[2] ^ s_star.pixInfo[l].lbp[2];
        e[l] += pixInfo[l].lbp[3] ^ s_star.pixInfo[l].lbp[3];
        e[l] += pixInfo[l].lbp[4] ^ s_star.pixInfo[l].lbp[4];
        e[l] += pixInfo[l].lbp[5] ^ s_star.pixInfo[l].lbp[5];
        e[l] += pixInfo[l].lbp[6] ^ s_star.pixInfo[l].lbp[6];
        e[l] += pixInfo[l].lbp[7] ^ s_star.pixInfo[l].lbp[7];

        l++;
      }
      e[l] /= 8;
    }
  }
  else std::cout << "ERROR !!";
}


/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired
 
  \param s_star : Desired visual feature.
  \param select : Not used.

*/
vpColVector
vpFeatureBitplane::error(const vpBasicFeature &s_star,
			  const unsigned int /* select */)
{
  /* static */ vpColVector e ; // warning C4640: 'e' : construction of local static object is not thread-safe
  
  error(s_star, e) ;
  
  return e ;

}




/*!

  Not implemented.

 */
void
vpFeatureBitplane::print(const unsigned int /* select */) const
{
  static int firsttime =0 ;

  if (firsttime==0)
  {
    firsttime=1 ;
    vpERROR_TRACE("not implemented") ;
    // Do not throw and error since it is not subject
    // to produce a failure
  }
}



/*!

  Not implemented.

 */
void
vpFeatureBitplane::display(const vpCameraParameters & /* cam */,
                            const vpImage<unsigned char> & /* I */,
                            const vpColor &/* color */,
                            unsigned int /* thickness */) const
{
 static int firsttime =0 ;

  if (firsttime==0)
  {
    firsttime=1 ;
    vpERROR_TRACE("not implemented") ;
    // Do not throw and error since it is not subject
    // to produce a failure
  }
}

/*!

  Not implemented.

 */
void
vpFeatureBitplane::display(const vpCameraParameters & /* cam */,
                            const vpImage<vpRGBa> & /* I */,
                            const vpColor &/* color */,
                            unsigned int /* thickness */) const
{
  static int firsttime =0 ;

  if (firsttime==0)
  {
    firsttime=1 ;
    vpERROR_TRACE("not implemented") ;
    // Do not throw and error since it is not subject
    // to produce a failure
  }
}


/*!
  Create an object with the same type.

  \code
  vpBasicFeature *s_star;
  vpFeatureBitplane s;
  s_star = s.duplicate(); // s_star is now a vpFeatureBitplane
  \endcode

*/
vpFeatureBitplane *vpFeatureBitplane::duplicate() const
{
  vpFeatureBitplane *feature = new vpFeatureBitplane ;
  return feature ;
}


/*
 * Local variables:
 * c-basic-offset: 2
 * End:
 */
