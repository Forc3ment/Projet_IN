/*! \example tutorial-template-tracker.cpp */
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/io/vpVideoReader.h>
//! [Include]
#include <visp3/tt/vpTemplateTrackerLBPForwardAdditional.h>
#include <visp3/tt/vpTemplateTrackerSSDForwardAdditional.h>
#include <visp3/tt/vpTemplateTrackerWarpHomography.h>
//! [Include]

int main(int argc, char** argv)
{
#if defined(VISP_HAVE_OPENCV) && (VISP_HAVE_OPENCV_VERSION >= 0x020100) || defined(VISP_HAVE_FFMPEG)
  std::string videoname = "../video/v2.mov";

  for (int i=0; i<argc; i++) {
    if (std::string(argv[i]) == "--videoname")
      videoname = std::string(argv[i+1]);
    else if (std::string(argv[i]) == "--help") {
      std::cout << "\nUsage: " << argv[0] << " [--name <video name>] [--help]\n" << std::endl;
      return 0;
    }
  }

  std::cout << "Video name: " << videoname << std::endl;

  vpImage<unsigned char> I;

  vpVideoReader g;
  g.setFileName(videoname);
  g.open(I);

#if defined(VISP_HAVE_X11)
  vpDisplayX display;
#elif defined(VISP_HAVE_GDI)
  vpDisplayGDI display;
#elif defined(VISP_HAVE_OPENCV)
  vpDisplayOpenCV display;
#else
  std::cout << "No image viewer is available..." << std::endl;
#endif

  display.init(I, 100, 100, "Template tracker");
  vpDisplay::display(I);
  vpDisplay::flush(I);

  //! [Construction]
  vpTemplateTrackerWarpHomography warpLBP;
  vpTemplateTrackerWarpHomography warpSSD;

  vpTemplateTrackerSSDForwardAdditional trackerSSD(&warpLBP);
  vpTemplateTrackerLBPForwardAdditional trackerLBP(&warpSSD);
  //! [Construction]

  trackerLBP.setSampling(2, 2);
  trackerLBP.setLambda(0.001);
  trackerLBP.setIterationMax(200);
  trackerLBP.setPyramidal(2, 1);

  trackerSSD.setSampling(2, 2);
  trackerSSD.setLambda(0.001);
  trackerSSD.setIterationMax(200);
  trackerSSD.setPyramidal(2, 1);

  //! [Init]
  trackerLBP.initClick(I);
  trackerSSD.initClick(I);
  //! [Init]

  while(1){
    g.acquire(I);
    vpDisplay::display(I);

    //! [Track]
    trackerLBP.track(I);
    trackerSSD.track(I);
    //! [Track]

    //! [Homography]
    vpColVector pLBP = trackerLBP.getp();
    vpHomography HLBP = warpLBP.getHomography(pLBP);
    //std::cout << "Homography: \n" << HLBP << std::endl;

    vpColVector pSSD = trackerSSD.getp();
    vpHomography HSSD = warpSSD.getHomography(pSSD);
    //std::cout << "Homography: \n" << HSSD << std::endl;
    //! [Homography]

    //! [Display]
    trackerLBP.display(I, vpColor::red);
    trackerSSD.display(I, vpColor::green);
    //! [Display]

    if (vpDisplay::getClick(I, false))
      break;

    vpDisplay::flush(I);
  }
#else
  (void)argc;
  (void)argv;
#endif
}
