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

int patch_errors(vpImage<unsigned char>& i1, vpImage<unsigned char>& i2)
{	
	int h = i1.getHeight(), w = i1.getWidth();
	vpImage<int> i11(h,w,0), i12(h,w,0), i13(h,w,0), i14(h,w,0), i15(h,w,0), i16(h,w,0), i17(h,w,0), i18(h,w,0);
	vpImage<int> i21(h,w,0), i22(h,w,0), i23(h,w,0), i24(h,w,0), i25(h,w,0), i26(h,w,0), i27(h,w,0), i28(h,w,0);
	lbp(i1, i11, i12, i13, i14, i15, i16, i17, i18);
	lbp(i2, i21, i22, i23, i24, i25, i26, i27, i28);
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

void process_graph(vpImage<unsigned char>& im, int hpatch, int wpatch, vpImage<unsigned char>& graph)
{
	int hpatch2 = hpatch/2, wpatch2 = wpatch/2;
	cout << "aa" << endl;
	// create & init central fixed patch
	vpImage<unsigned char> fixed_patch(hpatch,wpatch,0);
	for (int i = 0; i < hpatch; i++)
		for (int j = 0; j < wpatch; j++)
			fixed_patch[i][j] = im[(im.getHeight()-hpatch)/2 + i][(im.getWidth()-wpatch)/2 + j];
	cout << "a" << endl;
	// create moving patch
	vpImage<unsigned char> moving_patch(hpatch,wpatch,0);
	for (int i = hpatch2; i < im.getHeight()-hpatch2; i++)
		for (int j = wpatch2; j < im.getWidth()-wpatch2; j++) {
			// fill moving patch
			for (int ii = 0; ii < hpatch; ii++)
				for (int jj = 0; jj < wpatch; jj++) {
					//cout << "c" << i-hpatch2+ii << ' ' << ii << " " << jj << endl;
					moving_patch[ii][jj] = im[i-hpatch2+ii][j-wpatch2+jj];
				}
			//cout << "b" << i << endl;
			// fill graph
			graph[i][j] = patch_errors(fixed_patch, moving_patch);
		}
}

int main(int argc, char **argv)
{
	// creation du menu
	
	cout << "BINARY DESCRIPTOR" << endl;
	cout << "----" << endl ;
	
	string im;
	unsigned int h, w;	
	vpImage<unsigned char>  I1;
	
	im = "../image/baboon_grey.ppm";

	vpImageIo::read(I1, im);
	vpImage<unsigned char> I2(I1.getHeight(), I1.getWidth());
	//afficheImage(I1, 100, 100, "");
	h=I1.getHeight(); w=I1.getWidth();
	vpImage<int> i1(h,w,0), i2(h,w,0), i3(h,w,0), i4(h,w,0), i5(h,w,0), i6(h,w,0), i7(h,w,0), i8(h,w,0);
	lbp(I1, i1, i2, i3, i4, i5, i6, i7, i8);
	
	process_graph(I1, 100, 100, I2);

	afficheImage(I2, 100, 100, "");

	cout << "Fin du programme " << endl;
	return(0);
}