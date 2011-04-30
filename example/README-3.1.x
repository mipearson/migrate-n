Examples
--------
Several data type examples and their corresponding parmfiles are available here.
To test whether migrate works on your machine run simply:

time ../migrate-n parmfile.testml -nomenu
time ../migrate-n parmfile.testbayes -nomenu

Timing on my Macbook Pro (Intel Core 2 Duo 2.5 GHz with 2 cores,  6 MB L2 Cache, 
4 GB memory, 800 MH bus speed) running parmfile.testml and parmfile.testbayes
using: time ../migrate-n parmfile.x -nomenu, reporting real time

Version 3.2 October 2010
--------------------------------------------------------------------------------
compilation using:                     ML                      Bayes       Notes
make                                1m21.907s                  3m38.899s   - 
make                                1m9.828s                   2m57.881s   GrandCentral 
make thread                         3m19.316s                  6m50.137s   pthreads
make mpis (3 nodes)                 0m56.998s                  2m13.603s   -
make mpis (3 nodes)                 1m5.488ss                  3m0.103s    GrandCentral
--------------------------------------------------------------------------------
Timing using the FSU HPC cluster:
make (1 node)                       1m51.763s                  4m42.985s 
make mpis (11 nodes)                0m23.215s                  0m30.281s

On Mac: compiler gcc 4.1; the threaded version uses 5 threads (suboptimal on a 2-core 
machine, the MPI version uses 3 nodes [1 master + 2 worker]), GrandCentral is the Apple
architecture for fast threading introduced with MacOS 10.6. the make thread
seems not to work too well compared to the GrandCentral on a Mac, but may be the only way 
to parallelize the heated chains on LINUX (I have not tried any of this on Windows 
[if you are a specialist on windows parallelizing talk to me please]

Compare your runs with the runs labeled:
outfile-ml-saved
outfile-bayes-saved

They should look similar, although there might be differences because 
of optimization (compiler dependent) and hardware differences.
Differences in these (too) short run are due to different 
optimization on different machines, I doubt that it is
possible to compare between different computer architectures
but if the program does finish successfully on the test data
it is likely that the program will work on your data, too.


Peter Beerli
(beerli@fsu.edu)
May 2010
