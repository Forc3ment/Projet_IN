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

#include <visp3/tt/vpTemplateTrackerLBP.h>

vpTemplateTrackerLBP::vpTemplateTrackerLBP(vpTemplateTrackerWarp *warp)
  : vpTemplateTracker(warp), DI(), temp()
{
  dW.resize(2,nbParam);
  G.resize(nbParam);
  H.resize(nbParam,nbParam);
  HLM.resize(nbParam,nbParam);

  temp.resize(nbParam);

  X1.resize(2);
  X2.resize(2);
  DI.resize(2);

}

vpTemplateTrackerLBP::~vpTemplateTrackerLBP(){}

double vpTemplateTrackerLBP::getCost(const vpImage<unsigned char> &I, const vpColVector &tp)
{
  double erreur=0;
  double IW;
  int Nbpoint=0;

  Warp->computeCoeff(tp);

  int h = I.getHeight(), w = I.getWidth();
  vpImage<bool> I1(h,w,0), I2(h,w,0), I3(h,w,0), I4(h,w,0), I5(h,w,0), I6(h,w,0), I7(h,w,0), I8(h,w,0);
  lbp(I, I1, I2, I3, I4, I5, I6, I7, I8);

  for(unsigned int point=0;point<templateSize;point++)
  {
    int i=ptTemplate[point].y;
    int j=ptTemplate[point].x;
    X1[0]=j;X1[1]=i;
    Warp->computeDenom(X1,tp);
    Warp->warpX(X1,X2,tp);

    int j2=X2[0];
    int i2=X2[1];
    if((i2>=0)&&(j2>=0)&&(i2<I.getHeight()-1)&&(j2<I.getWidth()-1))
    {
      //double Tij=ptTemplate[point].val;
      if(!blur)
        IW=I.getValue(i2,j2);
      else
        IW=BI.getValue(i2,j2);
      //IW=getSubPixBspline4(I,i2,j2);

      //erreur+=((double)Tij-IW)*((double)Tij-IW);
      erreur += I1[i][j] ^ I1[i2][j2]; // ^ XOR
      erreur += I2[i][j] ^ I2[i2][j2];
      erreur += I3[i][j] ^ I3[i2][j2];
      erreur += I4[i][j] ^ I4[i2][j2];
      erreur += I5[i][j] ^ I5[i2][j2];
      erreur += I6[i][j] ^ I6[i2][j2];
      erreur += I7[i][j] ^ I7[i2][j2];
      erreur += I8[i][j] ^ I8[i2][j2];

      Nbpoint++;
    }
  }
  ratioPixelIn=(double)Nbpoint/(double)templateSize;

  if(Nbpoint==0)return 10e10;
  return erreur/Nbpoint;
}


double vpTemplateTrackerLBP::getLBP(const vpImage<unsigned char> &I, const vpColVector &tp)
{
  double erreur=0;
  double IW;
  unsigned int Nbpoint=0;

  if(pyrInitialised)
  {
    templateSize=templateSizePyr[0];
    ptTemplate=ptTemplatePyr[0];
  }

  Warp->computeCoeff(tp);
  for(unsigned int point=0;point<templateSize;point++)
  {
    int i=ptTemplate[point].y;
    int j=ptTemplate[point].x;
    X1[0]=j;X1[1]=i;
    Warp->computeDenom(X1,tp);
    Warp->warpX(X1,X2,tp);

    double j2=X2[0];
    double i2=X2[1];
    if((j2<I.getWidth()-1)&&(i2<I.getHeight()-1)&&(i2>0)&&(j2>0))
    {
      double Tij=ptTemplate[point].val;
      IW=I.getValue(i2,j2);
      //IW=getSubPixBspline4(I,i2,j2);
      erreur+=((double)Tij-IW)*((double)Tij-IW);
      Nbpoint++;
    }
  }
  if(Nbpoint==0)return 10e10;
  return erreur/Nbpoint;
}


void vpTemplateTrackerLBP::lbp(const vpImage<unsigned char> &im, vpImage<bool>& i1, vpImage<bool>& i2, vpImage<bool>& i3, vpImage<bool>& i4, 
  vpImage<bool>& i5, vpImage<bool>& i6, vpImage<bool>& i7, vpImage<bool>& i8)
{
  for (int i = 1; i < im.getHeight()-1; i++)
    for (int j = 1; j < im.getWidth()-1; j++) {
      i1[i][j] = im[i-1][j-1] < im[i][j];
      i2[i][j] = im[i-1][j] < im[i][j];
      i3[i][j] = im[i-1][j+1] < im[i][j];
      i4[i][j] = im[i][j-1] < im[i][j];
      i5[i][j] = im[i][j+1] < im[i][j];
      i6[i][j] = im[i+1][j-1] < im[i][j];
      i7[i][j] = im[i+1][j] < im[i][j];
      i8[i][j] = im[i+1][j+1] < im[i][j];
    }
}

