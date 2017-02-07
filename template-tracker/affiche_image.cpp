
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpDisplayX.h>

void afficheImage(vpImage<unsigned char> img, int posX, int posY, const char *title)
{
    vpDisplayX d(img, posX, posY, title);
    vpDisplay::display(img);
    vpDisplay::flush(img);
    vpDisplay::getClick(img);
    vpDisplay::close(img);
}

void afficheImage(vpImage<double> & D, int posX, int posY, const char *title)
{
  vpImage<unsigned char> I; // Image to display
  vpImageConvert::convert(D, I);
  afficheImage(I, posX, posY, title);
}

void afficheImage(vpImage<int> & img, int posX, int posY, const char *title)
{
  vpImage<unsigned char> I(img.getHeight(), img.getWidth()); // Image to display
  for (int i = 0; i < img.getHeight(); i++)
		for (int j = 0; j < img.getWidth(); j++)
			I[i][j] = img[i][j];
  afficheImage(I, posX, posY, title);
}