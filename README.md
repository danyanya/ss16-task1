# ss16-task1
# Lanczos algorithm MPI realisation #
Lanczos algogithm from algowiki-project.ru

## Compile and usage ##
You can compile it with MPIC++ compiller :

`mpic++ -std=c++0x -O2 lanczos-mpi.cpp -lm -o lanczos-mpi`


## Running ##
After compillation here in repo's directory lanczos-mpi had been created
And you can run it with next command :

`mpirun -np <X> lanczos-mpi`

, where X -- number of CPUs


## Authors ##

Made by Sliusar Daniil and Grigoriev Micheil
for superprac 2016 course and tested on Lomonosov supercomputer.

Moscow 2016 
