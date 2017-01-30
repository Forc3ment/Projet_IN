#include <visp/vpConfig.h>
#include <visp/vpDebug.h>
#include "affiche_image.cpp"
#include <visp3/gui/vpPlot.h>

using namespace std;

void lbp(vpImage<unsigned char>& im, vpImage<bool>& i1, vpImage<bool>& i2, vpImage<bool>& i3, vpImage<bool>& i4, 
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

int patch_errors(vpImage<bool>& i11, vpImage<bool>& i12, vpImage<bool>& i13, vpImage<bool>& i14, 
	vpImage<bool>& i15, vpImage<bool>& i16, vpImage<bool>& i17, vpImage<bool>& i18, vpImage<bool>& i21, vpImage<bool>& i22, vpImage<bool>& i23, vpImage<bool>& i24, 
	vpImage<bool>& i25, vpImage<bool>& i26, vpImage<bool>& i27, vpImage<bool>& i28)
{	
	int h = i11.getHeight(), w = i11.getWidth();
	int error = 0;
	for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++) {
			error += i11[i][j] ^ i21[i][j]; // ^ XOR
			error += i12[i][j] ^ i22[i][j];
			error += i13[i][j] ^ i23[i][j];
			error += i14[i][j] ^ i24[i][j];
			error += i15[i][j] ^ i25[i][j];
			error += i16[i][j] ^ i26[i][j];
			error += i17[i][j] ^ i27[i][j];
			error += i18[i][j] ^ i28[i][j];
		}
	return error;
}

void process_graph(vpImage<unsigned char>& im, int hpatch, int wpatch, vpImage<double>& graph)
{
	int h = im.getHeight(), w = im.getWidth();

	vpImage<bool> i1(h,w,0), i2(h,w,0), i3(h,w,0), i4(h,w,0), i5(h,w,0), i6(h,w,0), i7(h,w,0), i8(h,w,0);
	lbp(im, i1, i2, i3, i4, i5, i6, i7, i8);

	vpImage<bool> fixed_p1(hpatch,wpatch,0), fixed_p2(hpatch,wpatch,0), fixed_p3(hpatch,wpatch,0), fixed_p4(hpatch,wpatch,0)
				, fixed_p5(hpatch,wpatch,0), fixed_p6(hpatch,wpatch,0), fixed_p7(hpatch,wpatch,0), fixed_p8(hpatch,wpatch,0);

	vpImage<bool> p1(hpatch,wpatch,0), p2(hpatch,wpatch,0), p3(hpatch,wpatch,0), p4(hpatch,wpatch,0)
				, p5(hpatch,wpatch,0), p6(hpatch,wpatch,0), p7(hpatch,wpatch,0), p8(hpatch,wpatch,0);

	int hpatch2 = hpatch/2, wpatch2 = wpatch/2;
	
	// create & init central fixed patch
	for (int i = 0; i < hpatch; i++)
	{
		for (int j = 0; j < wpatch; j++)
		{
			fixed_p1[i][j] = i1[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p2[i][j] = i2[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p3[i][j] = i3[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p4[i][j] = i4[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p5[i][j] = i5[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p6[i][j] = i6[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p7[i][j] = i7[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
			fixed_p8[i][j] = i8[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
		}
	}			
	

	// create moving patch
	vpImage<unsigned char> moving_patch(hpatch,wpatch,0);
	for (int i = hpatch2; i < im.getHeight()-hpatch2; i++)
	{

		for (int j = wpatch2; j < im.getWidth()-wpatch2; j++) 
		{
			// fill moving patch
			for (int ii = 0; ii < hpatch; ii++)
			{
				for (int jj = 0; jj < wpatch; jj++) 
				{
					p1[ii][jj] = i1[i-hpatch2+ii][j-wpatch2+jj];
					p2[ii][jj] = i2[i-hpatch2+ii][j-wpatch2+jj];
					p3[ii][jj] = i3[i-hpatch2+ii][j-wpatch2+jj];
					p4[ii][jj] = i4[i-hpatch2+ii][j-wpatch2+jj];
					p5[ii][jj] = i5[i-hpatch2+ii][j-wpatch2+jj];
					p6[ii][jj] = i6[i-hpatch2+ii][j-wpatch2+jj];
					p7[ii][jj] = i7[i-hpatch2+ii][j-wpatch2+jj];
					p8[ii][jj] = i8[i-hpatch2+ii][j-wpatch2+jj];
				}
				// fill graph
				graph[i][j] = patch_errors(fixed_p1,fixed_p2,fixed_p3,fixed_p4,fixed_p5,fixed_p6,fixed_p7,fixed_p8,
										p1,p2,p3,p4,p5,p6,p7,p8);
			}
		}
	}
}

int main(int argc, char **argv)
{
	// creation du menu
	
	cout << "BINARY DESCRIPTOR" << endl;
	cout << "----" << endl ;
	
	string im;
	unsigned int h, w;	
	vpImage<unsigned char>  I0;
	
	im = "../image/baboon_grey.ppm";

	vpImageIo::read(I0, im);
	vpImage<double> result(I0.getHeight(), I0.getWidth(), 0);
	vpImage<unsigned char> result_gray(I0.getHeight(), I0.getWidth(), 0);

	//afficheImage(I1, 100, 100, "");
	h=I0.getHeight(); w=I0.getWidth();

	
	process_graph(I0, 30, 30, result);

  	// vpPlot A(1, 700, 700, 100, 100, "");
  	// A.initGraph(0,1);

  	//for (int i = 0; i < h; i++)
		//std::cout<< result[i][100] << std::endl;

	std::cout<< result[128][128] << std::endl;

  	vpImageConvert::convert(result,result_gray);

	afficheImage(result_gray, 100, 100, "");

	cout << "Fin du programme " << endl;
	return(0);
}