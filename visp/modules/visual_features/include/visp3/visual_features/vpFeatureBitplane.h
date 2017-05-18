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
 * Luminance based feature.
 *
 * Authors:
 * Eric Marchand
 *
 *****************************************************************************/

#ifndef vpFeatureBitplane_h
#define vpFeatureBitplane_h

#include <visp3/core/vpMatrix.h>
#include <visp3/visual_features/vpBasicFeature.h>
#include <visp3/core/vpImage.h>


/*!
  \file vpFeatureBitplane.h
  \brief Class that defines the image luminance visual feature

  For more details see \cite Collewet08c.
*/


/*!
  \class vpBitplane
  \brief Class that defines the luminance and gradient of a point

  \sa vpFeatureBitplane
*/


class VISP_EXPORT vpBitplane
{
 public:
  double x, y;   // point coordinates (in meter)
  int i, j;     // point coordinates (in pixel)
  bool* lbp ; // pixel sum of bitplanes
  double *Ix, *Iy ; // pixel gradient
  double Z; // pixel depth

};


/*!
  \class vpFeatureBitplane
  \ingroup group_visual_features
  \brief Class that defines the image luminance visual feature

  For more details see \cite Collewet08c.
*/

class VISP_EXPORT vpFeatureBitplane : public vpBasicFeature
{
 protected:
  //! FeaturePoint depth (required to compute the interaction matrix)
  //! default Z = 1m
  double Z ;

  //! Number of rows.
  unsigned int nbr ;
  //! Number of column.
  unsigned int nbc ;
  //! Border size.
  unsigned int bord ;
  
  //! Store the image (as a vector with intensity and gradient I, Ix, Iy) 
  vpBitplane *pixInfo ;
  int  firstTimeIn  ;

  int w, h; // save the w and h of the original image

 public:
  void buildFrom(vpImage<unsigned char> &I) ;

public: 

  void init() ;
  void init(unsigned int _nbr, unsigned int _nbc, double _Z) ;

  vpFeatureBitplane() ;
  vpFeatureBitplane(const vpFeatureBitplane& f) ;
  vpFeatureBitplane &operator=(const vpFeatureBitplane& f) ;

  //! Destructor.
  virtual ~vpFeatureBitplane()  ;

 public:
  vpCameraParameters cam ;
  void setCameraParameters(vpCameraParameters &_cam)  ;
  /*
    section Set/get Z
  */


  void set_Z(const double Z) ;
  double get_Z() const;

  void getAsImage(vpImage<unsigned char> & ret);

  /*
    vpBasicFeature method instantiation
  */

 
  vpMatrix  interaction(const unsigned int select = FEATURE_ALL);
  void      interaction(vpMatrix &L);

  vpColVector error(const vpBasicFeature &s_star,
                    const unsigned int select = FEATURE_ALL)  ;
  void error(const vpBasicFeature &s_star,
             vpColVector &e)  ;

  void print(const unsigned int select = FEATURE_ALL ) const ;

  vpFeatureBitplane *duplicate() const ;


  void display(const vpCameraParameters &cam,
               const vpImage<unsigned char> &I,
               const vpColor &color=vpColor::green, unsigned int thickness=1) const ;
  void display(const vpCameraParameters &cam,
               const vpImage<vpRGBa> &I,
               const vpColor &color=vpColor::green, unsigned int thickness=1) const ;

  //! Compute the error between a visual features and zero
  vpColVector error(const unsigned int select = FEATURE_ALL);

} ;



#endif
