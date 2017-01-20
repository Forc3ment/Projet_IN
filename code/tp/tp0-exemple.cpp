/****************************************************************************
 *
 *
 * Copyright (C) 2010 Université de Rennes 1. All rights reserved.
 *
 * This software was developed at:
 * Universite de Rennes 1
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 *
 * This file uses the ViSP library.
 *
 
 *
 *****************************************************************************/


#include <iostream>

/*!
 tp0_exemple.cpp
 Ouverture, affichage et sauvegarde d'une image avec la librairie VISP
 
 */


#include <visp/vpConfig.h>
#include <visp/vpDebug.h>

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpDisplayX.h>

using namespace std;



int main(int argc, char **argv)
{
	// creation du menu
	
	cout << "BINP TP : EXEMPLE MANIPULATION IMAGE AVEC VISP " << endl ;
	cout << "--" << endl ;
	
	string sIm;
	unsigned int h, w;
	
	// lecture (interactive) d'une image en niveau de gris
	vpImage<unsigned char>  I1;
	
	cout << "Nom de l'image (n.d.g.) : ";
	cin >> sIm; // Ex: ../images/lena.ppm	
	vpImageIo::read(I1,sIm) ;
	
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















