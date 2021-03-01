# monitor_VASP_convergence
plot_forces.py:
Plots the evolution of forces from VASP's OUTCAR file relative to the convergence criteria. Useful for analyzing optimization progression between or during job submissions.

By default, forces will be read from 'OUTCAR' and plotted. Alternative files can be specified with -o, --outcar [OUTCAR path] or -p, --poscar [POSCAR path]. If multiple job submissions are contained within the OUTCAR file (i.e. the VASP run was continued without purging output files from previous runs), the convergence criteria plotted will correspond to the EDIFFG parameter specified in the most recent job. If selective dynamics is employed in the POSCAR file, only forces along optimizable degrees of freedom will be plotted, since frozen degrees of freedom to not contribute to the force convergence threshold.

Example of intended usage: python plot_forces or python plot_forces -o ./my_OUTCAR -p ./my_POSCAR

plot_energies.py:
Plots the total energy, without entropy contributions, against either the optimization time or total number of optimization steps (if POTIM=0).

Example of intended usage: python plot_energies or python plot_energies -o ./my_OUTCAR

show_poscar.py
Plots the coordinates of all atoms in VASP structure file. By default, ./CONTCAR is ready if it exists and is not empty. Otherwise, ./POSCAR is read. Alternative file paths can be specified by -p, --poscar.

show_forces.py
Displays the coordinates of all atoms and their current forces. If selective dynamics is enabled, only optimizable degrees of freedom are plotted. Atom display is the same as in 'show_poscar.py'. Coordiantes and forces are read from ./CONTCAR and ./OUTCAR by default.


Compatible with VASP 5.4.4, Python 2.7.13 and 3.6.5, and VTST 3.2.
