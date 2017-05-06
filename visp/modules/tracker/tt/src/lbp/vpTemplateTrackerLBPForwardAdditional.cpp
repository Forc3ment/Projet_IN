#include <limits>   // numeric_limits

#include <visp3/tt/vpTemplateTrackerLBPForwardAdditional.h>
#include <visp3/core/vpImageTools.h>

vpTemplateTrackerLBPForwardAdditional::vpTemplateTrackerLBPForwardAdditional(vpTemplateTrackerWarp *warp)
  : vpTemplateTrackerLBP(warp), minimizationMethod(USE_NEWTON), p_prec(), G_prec(), KQuasiNewton()
{
  useCompositionnal=false;
}

void vpTemplateTrackerLBPForwardAdditional::trackNoPyr(const vpImage<unsigned char> &I)
{
  if(blur)
    vpImageFilter::filter(I, BI,fgG,taillef);
  vpImageFilter::getGradXGauss2D(I, dIx, fgG,fgdG,taillef);
  vpImageFilter::getGradYGauss2D(I, dIy, fgG,fgdG,taillef);

  dW=0;

  double lambda=lambdaDep;
  double IW,dIWx,dIWy;
  double Tij;
  unsigned int iteration=0;
  int i,j;
  double i2,j2;
  double alpha=2.;
  do
  {
    unsigned int Nbpoint=0;
    double erreur=0;
    G=0;
    H=0 ;
    Warp->computeCoeff(p);
    for(unsigned int point=0;point<templateSize;point++)
    {
      i=ptTemplate[point].y;
      j=ptTemplate[point].x;
      X1[0]=j;X1[1]=i;

      Warp->computeDenom(X1,p);
      Warp->warpX(X1,X2,p);

      j2=X2[0];i2=X2[1];
      if((i2>=0)&&(j2>=0)&&(i2<I.getHeight()-1)&&(j2<I.getWidth()-1))
      {
        Tij=ptTemplate[point].val;

        if(!blur)
          IW=I.getValue(i2,j2);
        else
          IW=BI.getValue(i2,j2);

        dIWx=dIx.getValue(i2,j2);
        dIWy=dIy.getValue(i2,j2);
        Nbpoint++;
        //Calcul du Hessien
        Warp->dWarp(X1,X2,p,dW);
        double *tempt=new double[nbParam];
        for(unsigned int it=0;it<nbParam;it++)
          tempt[it]=dW[0][it]*dIWx+dW[1][it]*dIWy;

        for(unsigned int it=0;it<nbParam;it++)
          for(unsigned int jt=0;jt<nbParam;jt++)
            H[it][jt]+=tempt[it]*tempt[jt];

        double er=(Tij-IW);
        for(unsigned int it=0;it<nbParam;it++)
          G[it]+=er*tempt[it];

        erreur+=(er*er);
        delete[] tempt;
      }


    }
    if(Nbpoint==0) {
      //std::cout<<"plus de point dans template suivi"<<std::endl;
      throw(vpTrackingException(vpTrackingException::notEnoughPointError, "No points in the template"));
    }

    vpMatrix::computeHLM(H,lambda,HLM);
    try
    {
      dp=1.*HLM.inverseByLU()*G;
    }
    catch(vpException &e)
    {
      throw(e);
    }

    switch(minimizationMethod)
    {
    case vpTemplateTrackerLBPForwardAdditional::USE_LMA:
    {
      vpColVector p_test_LMA(nbParam);
      p_test_LMA=p+1.*dp;
      erreur=-getCost(I,p);
      double erreur_LMA=-getCost(I,p_test_LMA);
      if(erreur_LMA<erreur)
      {
        p=p_test_LMA;
        lambda=(lambda/10.<1e-6)?lambda/10.:1e-6;
      }
      else
      {
        lambda=(lambda*10.<1e6)?1e6:lambda*10.;
      }
    }
      break;
    case vpTemplateTrackerLBPForwardAdditional::USE_GRADIENT:
    {
      dp=gain*0.000001*G;
      if(useBrent)
      {
        alpha=2.;
        computeOptimalBrentGain(I,p,erreur,dp,alpha);
        dp=alpha*dp;
      }
      p += 1. * dp;
      break;
    }

    case vpTemplateTrackerLBPForwardAdditional::USE_QUASINEWTON:
    {
      if(iterationGlobale!=0)
      {
        vpColVector s_quasi=p-p_prec;
        vpColVector y_quasi=G-G_prec;
        double s_scal_y=s_quasi.t()*y_quasi;
        //if(s_scal_y!=0)//BFGS
        //	KQuasiNewton=KQuasiNewton-(s_quasi*y_quasi.t()*KQuasiNewton+KQuasiNewton*y_quasi*s_quasi.t())/s_scal_y+(1.+y_quasi.t()*(KQuasiNewton*y_quasi)/s_scal_y)*s_quasi*s_quasi.t()/s_scal_y;
        //if(s_scal_y!=0.0)//DFP
        if(std::fabs(s_scal_y) > std::numeric_limits<double>::epsilon()) //DFP
          KQuasiNewton=KQuasiNewton+0.001*(s_quasi*s_quasi.t()/s_scal_y-KQuasiNewton*y_quasi*y_quasi.t()*KQuasiNewton/(y_quasi.t()*KQuasiNewton*y_quasi));
      }
      dp=-KQuasiNewton*G;
      p_prec=p;
      G_prec=G;
      p-=1.01*dp;
    }
      break;

    case vpTemplateTrackerLBPForwardAdditional::USE_NEWTON:
    default:
    {
      if(useBrent)
      {
        alpha=2.;
        computeOptimalBrentGain(I,p,erreur,dp,alpha);
        dp=alpha*dp;
      }

      p+=1.*dp;
      break;
    }
    }

    iteration++;
    iterationGlobale++;
  }
  while( /*( erreur_prec-erreur<50) && */(iteration < iterationMax));

  //std::cout<<"erreur "<<erreur<<std::endl;
  nbIteration=iteration;
}