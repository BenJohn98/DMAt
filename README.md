# DMAT
Cpp code for calculating cosmological dark matter densities of thermal relic DM candidate particles

Compilation and running instructions included in this README file. 


# Compilation
Recommend compiling using ubuntu terminal available on linux systems or through 'windows subsystem for linux' with windows 10. 
Requires installation of gsl libraries and g++ compiler.

To install compiler: $sudo apt-get install g++

To install GSL libraries: $sudo apt-get install libgsl-dev

To compile: g++ -o example DMAT.cpp

# Running Instructions
After compilation, run using bash: ./example

Code does not require input options and any variables should be changed directly within the .cpp file, keeping in mind to recompile after making any changes. 
gnuplot was used within the linux terminal for generating plots.

# Notes
Documentation with background theory and deriavation of relevant equations, as well as output for specific conditions provided in 'DMATdoc.pdf'. 
No time to generate plots for scalar singlet dark matter-- maybe try in future.
