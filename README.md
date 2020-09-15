# RhoToyMc
Toy Monte Carlo generator for UPC Rho0 analysis

This project includes all dependencies for the rho toy mc generator along with CMAKE build files


TLDR
```
git clone https://github.com/jdbrice/RhoToyMc.git
cd RhoToyMc
mkdir build
cd build
cmake ..
cmake --build .
cd ..
./build/toy gen.xml
```

This generates a file `GENV2_TEST.root` with the simulated Rho distributions  
There are two analysis scripts (I will add more)  
proj.C : Project a particular fourier component for a given pT or mass range (as a function of either pT or mass)  
plot_pairpt.C : plots the rho pair pt showing the contribution from coherent, incoherent, pure pipi and from pimu decays  
  
NOTE: the only source file that is relevant is the RhoGen.h  
  
