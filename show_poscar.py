from numpy import array,dot
from getopt import getopt
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from os.path import exists,getsize

def show_poscar(poscar):
    lv, coord, atomtypes, atomnums = parse_poscar(poscar)[:4]
    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    for i in range(len(atomtypes)):
        tempvar=[[],[],[]]
        for j in range(atomnums[i]):
            for k in range(3):
                tempvar[k].append(coord[j+sum(atomnums[:i])][k])
        ax.scatter(tempvar[0],tempvar[1],tempvar[2],s=2000/max([norm(lv[j]) for j in range(3)]),label=atomtypes[i])
    ax.set_xlabel('x / $\AA$')
    ax.set_ylabel('y / $\AA$')
    ax.set_zlabel('z / $\AA$')
    ax.set_xlim(0,norm(dot(array([1.0,0.0,0.0]),lv)))
    ax.set_ylim(0,norm(dot(array([0.0,1.0,0.0]),lv)))
    ax.set_zlim(0,norm(dot(array([0.0,0.0,1.0]),lv)))
    fig.legend()
    plt.show()
    

def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if 'Direct' in lines[7] or 'Cartesian' in lines[7]:
            start=8
            mode=lines[7].split()[0]
        else:
            mode=lines[8].split()[0]
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        if mode!='Cartesian':
            for i in range(sum(atomnums)):
                for j in range(3):
                    while coord[i][j]>1.0 or coord[i][j]<0.0:
                        if coord[i][j]>1.0:
                            coord[i][j]-=1.0
                        elif coord[i][j]<0.0:
                            coord[i][j]+=1.0
                coord[i]=dot(coord[i],latticevectors)
            
    #latticevectors formatted as a 3x3 array
    #coord holds the atomic coordinates with shape ()
    try:
        return latticevectors, coord, atomtypes, atomnums, seldyn
    except NameError:
        return latticevectors, coord, atomtypes, atomnums
    
if __name__ == '__main__':
    if exists('./CONTCAR'):
        if getsize('./CONTCAR')>0:
            ifile='./CONTCAR'
        else:
            ifile='./POSCAR'
    else:
        ifile='./POSCAR'
    short_opts='hi:'
    long_opts=['help','input=']
    try:
        opts,args=getopt(sys.argv[1:],short_opts,long_opts)
    except getopt.GetoptError:
        print('error in command line syntax')
        sys.exit(2)
    for i,j in opts:
        if i in ['-h','--help']:
            print('''
    Note: for options with multiple values, seperate the values with commas
                  
    these options take a value:
    -i, --input                    specify an input other than ./CONTCAR
    help options:
    -h, --help                    display this error message
    ''')
            sys.exit()
        if i in ['-i','--input']:
            ifile=j
    
    show_poscar(ifile)