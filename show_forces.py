import matplotlib.pyplot as plt
import sys
import getopt
from numpy import array,dot
from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D

def show_forces(outcar,poscar,**args):
    if 'show_all' in args and args['show_all']:
        show_all=True
    else:
        show_all=False

    try:
        lv,coord,atomtypes,atomnums,seldyn = parse_poscar(poscar)
    except IndexError or FileNotFoundError:
        print('missing POSCAR or CONTCAR file. exiting...')
        sys.exit()
        
    forces,time,tol=parse_forces(outcar,seldyn=seldyn)
    forces=[i[-1] for i in forces]
    
    origins=[[],[],[]]
    vectors=[[],[],[]]
    counter=0    
    for i in range(sum(atomnums)):
        for j in range(3):
            if seldyn[i][j]=='T':
                origins[j].append(coord[i][j])
                vectors[j].append(forces[j][counter]/tol)
            elif 'T' in seldyn[i]:
                origins[j].append(coord[i][j])
                vectors[j].append(0.0)
        if 'T' in seldyn[i]:
            counter+=1
    
    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    for i in range(len(atomtypes)):
        tempvar=[[],[],[]]
        for j in range(atomnums[i]):
            for k in range(3):
                tempvar[k].append(coord[j+sum(atomnums[:i])][k])
        ax.scatter(tempvar[0],tempvar[1],tempvar[2],s=2000/max([norm(lv[j]) for j in range(3)]),label=atomtypes[i])
    ax.quiver(origins[0],origins[1],origins[2],vectors[0],vectors[1],vectors[2])
    ax.set_xlim(0,norm(dot(array([1.0,0.0,0.0]),lv)))
    ax.set_ylim(0,norm(dot(array([0.0,1.0,0.0]),lv)))
    ax.set_zlim(0,norm(dot(array([0.0,0.0,1.0]),lv)))
    plt.show()
    
def parse_forces(ifile,**args):
    if 'seldyn' in args:
        seldyn=args['seldyn']
    else:
        seldyn='none'
    
    time=[]
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
                    if len(time)==0:
                        time.append(0.0)
                    else:
                        time.append(time[-1]+abs(potim))
    except:
        print('error reading OUTCAR')
        sys.exit(1)
        
    if len(time)==0:
        print('zero ionic steps read from OUTCAR')
        sys.exit()
        
    return forces,time,tol

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
    show_all=False
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
        if i in ['-a','--all']:
            show_all=True
    show_forces(outcar,poscar,show_all=show_all)