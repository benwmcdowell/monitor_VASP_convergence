import matplotlib.pyplot as plt
import sys
import getopt
from numpy import array,dot,percentile,average
from numpy.linalg import norm

def main(outcar,poscar,**args):
    if 'quiet' in args and args['quiet']:
        quiet=True
    else:
        quiet=False

    try:
        seldyn = parse_poscar(poscar)[4]
    except ValueError or FileNotFoundError:
        seldyn='none'
    time=[0]
    forces=[[],[],[],[]]
    try:
        with open(outcar,'r') as file:
            searching=True
            while searching:
                line=file.readline()
                if not line:
                    break
                if 'EDIFFG' in line:
                    line=line.split()
                    tol=abs(float(line[line.index('EDIFFG')+2]))
                if 'POTIM' in line:
                    line=line.split()
                    potim=float(line[line.index('POTIM')+2])
                    if potim==0.0:
                        potim=-1.0
                if 'NIONS' in line:
                    line=line.split()
                    atomnum=int(line[line.index('NIONS')+2])
                    if seldyn=='none':
                        seldyn=['TTT' for i in range(atomnum)]
                elif 'TOTAL-FORCE' in line:
                    line=file.readline()
                    temp_forces=[[],[],[],[]]
                    for i in range(atomnum):
                        line=file.readline().split()
                        tempvar=[]
                        for j in range(3,6):
                            if seldyn[i][j-3]=='T':
                                temp_forces[j-3].append(abs(float(line[j])))
                                tempvar.append(abs(float(line[j])))
                        if len(tempvar)>0:
                            temp_forces[3].append(norm(array(tempvar)))
                    for i in range(4):
                        forces[i].append(temp_forces[i])
                    if len(forces[0])>1:
                        time.append(time[-1]+abs(potim))
    except:
        print('error reading OUTCAR')
        sys.exit(1)
        
    if len(time)==1:
        print('zero ionic steps read from OUTCAR')
        sys.exit()
    
    minima=[[min(j) for j in i] for i in forces]
    averages=minima=[[average(j) for j in i] for i in forces]
    maxima=minima=[[max(j) for j in i] for i in forces]
    upperq=minima=[[percentile(j,75) for j in i] for i in forces]
    lowerq=minima=[[percentile(j,25) for j in i] for i in forces]
    if not quiet:
        data_labels=['minimum','lower quartile','average','upper quartile','maximum']
        data_sets=[minima,lowerq,averages,upperq,maxima]
    else:
        data_labels=['minimum','average','maximum']
        data_sets=[minima,averages,maxima]
    #each component and the total force are plotted on their own subplot, along with the convergence criteria set by EDIFFG
    fig,axs=plt.subplots(4,1,sharex=True,figsize=(14,8))
    for i,j in zip(range(4),['_x','_y','_z','_{total}']):
        for k,l in zip(data_labels,data_sets):
            axs[i].scatter(time,l[i],label=k)
        axs[i].plot([time[0],time[-1]],[tol,tol],linestyle='dashed',label='convergence')
        axs[i].set(ylabel='$F{}$'.format(j)+' / eV $\AA^{-1}$')
        max_range=max(maxima[i])-min(minima[i])
        axs[i].set_ylim(bottom=min(minima[i])-0.05*max_range,top=max(maxima[i])+0.05*max_range)
    if potim>0.0:
        axs[-1].set(xlabel='optimization time / fs')
    else:
        axs[-1].set(xlabel='optimization steps')
    handles, labels = axs[2].get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.01,0.5), loc='right')
    plt.show()

def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if lines[7].split()[0] == 'Direct':
            start=8
        else:
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        for i in range(sum(atomnums)):
            coord[i]=dot(latticevectors,coord[i])
            
    #latticevectors formatted as a 3x3 array
    #coord holds the atomic coordinates with shape ()
    try:
        return latticevectors, coord, atomtypes, atomnums, seldyn
    except NameError:
        return latticevectors, coord, atomtypes, atomnums

if __name__=='__main__':
    outcar='./OUTCAR'
    poscar='./POSCAR'
    quiet=False
    try:
        opts,args=getopt.getopt(sys.argv[1:],'ho:p:q',['help','outcar=','poscar','quiet'])
    except getopt.GetoptError:
        print('error in command line syntax')
        sys.exit(2)
    for i,j in opts:
        if i in ['-h','--help']:
            print('''
input options:
    -o, --outcar          specify a path to the OUTCAR file other than ./OUTCAR
    -p, --poscar          specify an path to the POTCAR file other than ./POSCAR
    -q, --quiet           suppresses plotting of quartiles for a less crowded output
    
help options:
    -h, --help            display this help message
                  ''')
            sys.exit()
        if i in ['-o','--outcar']:
            outcar=j
        if i in ['-p','--poscar']:
            poscar=j
        if i in ['-q','--quiet']:
            quiet=True
    main(outcar,poscar,quiet=quiet)
