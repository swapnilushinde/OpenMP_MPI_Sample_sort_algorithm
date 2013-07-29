														README.TXT
														
Files attached -
						number.txt 						: 64 random numbers input file.
						sampleSort.c 					: MPI Only source code
						sampleSortOpenMP.c 				: MPI OpenMP hybrid code.
						Output_64_4_4_MPI_Only.txt 		: Output of MPI only implementation
						Output_64_4_4_MPI_OprnMP.txt 	: Output of MPI OPenMP implementation
						Project2.doc 					: Project details and implementation details
						

Input file structure - 
	First line has total number of elements. (n)
	Second line has total number of threads. (t)
	All numbers to be sorted are in single column. Each element at new line.

Compilation - 
	Compiling sampleSort.c :
		mpicc -o sampleSort sampleSort.c  -std=c99 -lm
	Executing sampleSort.c :
		mpiexec -n 4 ./sampleSort
		
	Compiling sampleSortOpenMP.c :
		mpicc -o sampleSortOpenMP sampleSortOpenMP.c -fopenmp -std=c99 -lm
	Executing sampleSortOpenMP.c :
		mpiexec -n 4 ./sampleSortOpenMP

The code has been tested for various number of random inputs. (32,64,256,1024,4096,16384) and number of processors(2,4,8) and local threads(4,8). 
It worked correctly for all combination of inputs..


Swapnil U. Shinde
sshinde2@gmu.edu

		