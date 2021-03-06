# MS2Nanotechnology

This is the simulation code used to model VLP assembly in the paper "Programmable Polymorphism of a Virus-Like Particle" (submitted to Nature Communications). This stochastic model uses a Gillespie algorithm. The code is written in Fortran 90 and can be compiled using gfortran. The -O4 flag should be used to optimise the executable output.

To compile, use: gfortran -o VLP.x VLP_assembly.f90 -O4

To run the .x file, use: ./VLP.x

This code as currently configured runs one simulation of VLP assembly, using an integer seed provided in the code. The results are written to the file "VLP.dat"

On a desktop computer (~2Ghz processor) the code should take around 10mins to run and uses <1Gb of RAM.
 
The code is available for use by anyone, however if in an academic context, we would appreciate a citation: Biela, A.P., Naskalska, A., Fatehi, F., Twarock, R., Heddle, J.G. Programmable polymorphism of a virus-like particle. Commun. Mater., 3, 7 (2022). https://doi.org/10.1038/s43246-022-00229-3
