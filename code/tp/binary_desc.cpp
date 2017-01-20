#include <visp/vpConfig.h>
#include <visp/vpDebug.h>

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpDisplayX.h>

using namespace std;

vpImage<vpImage<int>*>* lbp(vpImage<int>& im)
{
	vpImage<vpImage<int>*>* out = new vpImage<vpImage<int>*>(im.getHeight(), im.getWidth());
	for (int i = 1; i < im.getHeight()-1; i++)
		for (int j = 1; j < im.getWidth()-1; j++) {
			out[i][j] = new vpImage<int>(3, 3);
			out[i][j][0][0] = (im[i-1][j-1] < im[i][j]);
			out[i][j][0][1] = (im[i-1][j] < im[i][j]);
			out[i][j][0][2] = (im[i-1][j+1] < im[i][j]);
			out[i][j][1][0] = (im[i][j-1] < im[i][j]);
			out[i][j][1][2] = (im[i][j+1] < im[i][j]);
			out[i][j][2][0] = (im[i+1][j-1] < im[i][j]);
			out[i][j][0][1] = (im[i+1][j] < im[i][j]);
			out[i][j][2][2] = (im[i+1][j+1] < im[i][j]);
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
	
	im = "image.png";

	vpImageIo::read(I1, im);	
	h=I1.getHeight(); w=I1.getWidth();


	
	cout << "Lecture " << sIm << " (" << h << ", " << w << ")" << endl;
	
	vpDisplayX d1(I1,100,100) ;
	vpDisplay::setTitle(I1, "Image n.d.g.");
	vpDisplay::display(I1);
	vpDisplay::flush(I1) ;	
	vpDisplay::getClick(I1);
	
	// lecture (interactive) d'une image couleur
	vpImage<vpRGBa>  I2;
	cout << "Nom de l'image (couleur) : ";
	cin >> sIm; // Ex: ../images/lena.pgm	
	vpImageIo::read(I2,sIm) ;
	
	h=I2.getHeight(); w=I2.getWidth();
	
	cout << "Lecture " << sIm << " (" << h << ", " << w << ")" << endl;
	
	vpDisplayX d2(I2,500,100) ;
	vpDisplay::setTitle(I2, "Image couleur");
	vpDisplay::display(I2);
	vpDisplay::flush(I2) ;	
	vpDisplay::getClick(I2);
	
	// modification et sauvegarde de l'image
	for (int i=50; i<=150; i++) 
		for (int j=50; j<=150; j++){
			I2[i][j].R=0;
			I2[i][j].G=0;
			I2[i][j].B=0;
		}

	cout << "Sauvegarder l'image sous : ";
	cin >> sIm; // Ex: ../tp0_results/lena_modif.pgm	
	vpImageIo::write(I2,sIm) ;
		
	vpDisplay::close(I1) ;
	vpDisplay::close(I2) ;
	
	cout << "Fin du programme " << endl ;
	return(0);
}