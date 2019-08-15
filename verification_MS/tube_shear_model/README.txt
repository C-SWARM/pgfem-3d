MULTI-SCALE simulation of a shear model. 

The micro model consists of a void medium (~20%). The macro model consists 2 stiff blocks. The micro solution was verified with the reference solution given by Hirschberger et al. (2008).

There will be 2 stack input folders (macro and micro) placed separately.

To run multi-scale simulation, navigate to the folder /scripts/ and run "qsub server-run.sh" 

Numbers of CPUs used for the macro: 2
Numbers of CPUs used for the micro: 3
Number of steps: 3

The simulated solution of force-time curve is stored in the folder /scripts/force.dat
