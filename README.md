# Projet_IN
Camera tracking project

------------- How to install -------------

Clone latest VISP git repo.
Clone our git repo at the same place as VISP one.
Create a folder build next to VISP folder.
Go to this folder.
Use ccmake ../visp
Use make, wait a bit ...
and then in the build directory use ./tutorial/tracking/template-tracker/tutorial-template-tracker-test1tracker and enjoy ;).
(you need to select a polygonial patch with left clicks, the last point of it needing a right click)

```
git clone https://github.com/lagadic/visp.git
git clone https://github.com/Forc3ment/Projet_IN.git
mv visp Projet_IN/visp
cd Projet_IN
mdkir build
cd build
ccmake ../visp
make
./tutorial/tracking/template-tracker/tutorial-template-tracker-test1tracker
```
