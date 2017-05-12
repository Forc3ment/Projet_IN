#ifndef vpTemplateTrackerLBPForwardAdditional_hh
#define vpTemplateTrackerLBPForwardAdditional_hh

#include <visp3/tt/vpTemplateTrackerLBP.h>

/*!
  \ingroup group_tt_tracker
  The algorithm implemented in this class is described in \cite Baker04a and \cite Marchand16a.
 */
class VISP_EXPORT vpTemplateTrackerLBPForwardAdditional: public vpTemplateTrackerLBP
{
  public:
  /*! Minimization method. */
  typedef enum {
    USE_NEWTON,
    USE_LMA,
    USE_GRADIENT,
    USE_QUASINEWTON
  } vpMinimizationTypeLBPForwardAdditional;

  private:
    vpMinimizationTypeLBPForwardAdditional minimizationMethod;
    //valeur pour calculer Quasi_Newton
    vpColVector        p_prec;
    vpColVector        G_prec;
    vpMatrix           KQuasiNewton;

  protected:
    void  initHessienDesired(const vpImage<unsigned char> &/*I*/){}
    void  trackNoPyr(const vpImage<unsigned char> &I);

  public:
          vpTemplateTrackerLBPForwardAdditional(vpTemplateTrackerWarp *warp);

    void  setMinimizationMethod(vpMinimizationTypeLBPForwardAdditional method){minimizationMethod=method;}
};
#endif
