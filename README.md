# plot_forces_from_VASP
Plots the evolution of forces from VASP's OUTCAR file relative to the convergence criteria. Useful for analyzing optimization progression between or during job submissions.

By default, forces will be read from 'OUTCAR' and plotted. Alternative files can be specified with -i or --input. If multiple job submissions are contained within the OUTCAR file (i.e. the VASP run was continued without purging output files from previous runs), the convergence criteria plotted will correspond to the EDIFFG parameter specified in the most recent job.

Example of intended usage: python plot_forces or python plot_forces -i my_OUTCAR

Compatible with: 
  VASP 5.4.4
  Python 2.7.13 and 3.6.5
  VTST 3.2
