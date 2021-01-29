import matplotlib.pyplot as plt
import sys
import getopt

def main(outcar):
    time=[0.0]
    energies=[]
    with open(outcar,'r') as file:
        searching=True
        while searching:
            line=file.readline()
            if not line:
                break
            if 'POTIM' in line:
                line=line.split()
                potim=float(line[line.index('POTIM')+2])
                if potim==0.0:
                    potim=-1.0
            elif 'energy(sigma->0)' in line:
                energies.append(float(line.split()[-1]))
                if len(energies)>1:
                    time.append(time[-1]+abs(potim))
    
    plt.figure()
    plt.scatter(time,energies)
    plt.ylabel('total energy / eV')
    if potim>0.0:
        plt.xlabel('optimization time / fs')
    else:
        plt.xlabel('optimization steps')
    plt.show()

if __name__=='__main__':
    outcar='./OUTCAR'
    try:
        opts,args=getopt.getopt(sys.argv[1:],'ho:',['help','outcar='])
    except getopt.GetoptError:
        print('error in command line syntax')
        sys.exit(2)
    for i,j in opts:
        if i in ['-h','--help']:
            print('input options:\n\t-i, --input\t\tspecify an input file name other than ''OUTCAR''\n\nhelp options:\n\t-h, --help\t\tdisplay this help message')
            sys.exit()
        if i in ['-o','--outcar']:
            outcar=j
    main(outcar)