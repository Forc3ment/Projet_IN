/****************************************************************************
 *
 * This file is part of the ViSP software.
 * Copyright (C) 2005 - 2017 by Inria. All rights reserved.
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
 * Template tracker.
 *
 * Authors:
 * Amaury Dame
 * Aurelien Yol
 * Fabien Spindler
 *
 *****************************************************************************/
/*!
 \file vpTemplateTrackerLBP.h
 \brief
*/

#ifndef vpTemplateTrackerLBP_hh
#define vpTemplateTrackerLBP_hh

#include <math.h>

#include <visp3/tt/vpTemplateTracker.h>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpDisplay.h>
#include <visp3/core/vpImageTools.h>
#include <visp3/core/vpIoTools.h>
#include <visp3/core/vpImageTools.h>
#include <visp3/core/vpImageFilter.h>
#include <visp3/core/vpMath.h>
#include <visp3/vision/vpHomography.h>
#include <visp3/tt/vpTemplateTrackerWarp.h>

/*!
  \class vpTemplateTrackerLBP
  \ingroup group_tt_tracker
*/
class VISP_EXPORT vpTemplateTrackerLBP: public vpTemplateTracker
{
  protected:
    vpRowVector     DI;
    vpRowVector     temp;

  protected:
            double  getCost(const vpImage<unsigned char> &I, const vpColVector &tp);
            double  getCost(const vpImage<unsigned char> &I){ return getCost(I,p); }
    virtual void    initHessienDesired(const vpImage<unsigned char> & I) = 0;
    virtual void    trackNoPyr(const vpImage<unsigned char> &I) = 0;

  public:
    vpTemplateTrackerLBP(vpTemplateTrackerWarp *warp);
    ~vpTemplateTrackerLBP();

    double  getLBP(const vpImage<unsigned char> &I, const vpColVector &tp);
    void    setGain(double g){ gain=g; }
    void    lbp(const vpImage<unsigned char>& im, vpImage<bool>& i1, vpImage<bool>& i2, vpImage<bool>& i3, vpImage<bool>& i4, 
                vpImage<bool>& i5, vpImage<bool>& i6, vpImage<bool>& i7, vpImage<bool>& i8);
};

#endif

