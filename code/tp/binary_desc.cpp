#include <visp/vpConfig.h>
#include <visp/vpDebug.h>
#include "affiche_image.cpp"

using namespace std;

void lbp(vpImage<unsigned char>& im, vpImage<int>& i1, vpImage<int>& i2, vpImage<int>& i3, vpImage<int>& i4, 
	vpImage<int>& i5, vpImage<int>& i6, vpImage<int>& i7, vpImage<int>& i8)
{
	for (int i = 1; i < im.getHeight()-1; i++)
		for (int j = 1; j < im.getWidth()-1; j++) {
			i1[i][j] = (im[i-1][j-1] < im[i][j])*255;
			i2[i][j] = (im[i-1][j] < im[i][j])*255;
			i3[i][j] = (im[i-1][j+1] < im[i][j])*255;
			i4[i][j] = (im[i][j-1] < im[i][j])*255;
			i5[i][j] = (im[i][j+1] < im[i][j])*255;
			i6[i][j] = (im[i+1][j-1] < im[i][j])*255;
			i7[i][j] = (im[i+1][j] < im[i][j])*255;
			i8[i][j] = (im[i+1][j+1] < im[i][j])*255;
		}
}

int main(int argc, char **argv)
{
	// creation du menu
	
	cout << "BINARY DESCRIPTOR" << endl ;
	cout << "----" << endl ;
	
	string im;
	unsigned int h, w;	
	vpImage<unsigned char>  I1;
	
	im = "../image/baboon_grey.ppm";

	vpImageIo::read(I1, im);
	//afficheImage(I1, 100, 100, "");
	h=I1.getHeight(); w=I1.getWidth();
	vpImage<int> i1(h,w,0), i2(h,w,0), i3(h,w,0), i4(h,w,0), i5(h,w,0), i6(h,w,0), i7(h,w,0), i8(h,w,0);
	lbp(I1, i1, i2, i3, i4, i5, i6, i7, i8);

	afficheImage(i1, 100, 100, "");
	
	cout << "Fin du programme " << endl ;
	return(0);
}