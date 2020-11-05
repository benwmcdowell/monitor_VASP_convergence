# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 14:16:06 2020

@author: Ben
"""

import matplotlib.pyplot as plt
import sys
import getopt

def main(outcar):
    tempx=[0]
    max_force=[[],[],[]]
    avg_force=[[],[],[]]
    min_force=[[],[],[]]
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
                if 'NIONS' in line:
                    line=line.split()
                    atomnum=int(line[line.index('NIONS')+2])
                       
                elif 'TOTAL-FORCE' in line:
                    line=file.readline()
                    temp_min=[1,1,1]
                    temp_max=[0,0,0]
                    temp_avg=[0,0,0]
                    for i in range(atomnum):
                        line=file.readline().split()
                        for j in range(3,6):
                            tempvar=abs(float(line[j]))
                            if tempvar<temp_min[j-3]:
                                temp_min[j-3]=tempvar
                            if tempvar>temp_max[j-3]:
                                temp_max[j-3]=tempvar
                            temp_avg[j-3]+=tempvar/atomnum
                    for i in range(3):
                        max_force[i].append(temp_max[i])
                        min_force[i].append(temp_min[i])
                        avg_force[i].append(temp_avg[i])
                    if len(avg_force[0])>1:
                        tempx.append(tempx[-1]+potim)
        fig,ax=plt.subplots(3,1,sharex=True)
        for i,j in zip(range(3),['x','y','z']):
            for k,l in zip([max_force[i],min_force[i],avg_force[i]],['max force','min force','avg force']):
                ax[i].scatter(tempx,k,label=l)
            ax[i].plot([tempx[0],tempx[-1]],[tol,tol],linestyle='dashed',label='convergence')
            ax[i].set(ylabel='$F_{}$'.format(j)+' / eV $\AA^{-1}$')
            ax[i].set(ylim=(0-max(max_force[i])*0.05,max(max_force[i])*1.05))
        ax[2].set(xlabel='optimization time / fs')
        handles, labels = ax[2].get_legend_handles_labels()
        fig.legend(handles, labels, loc='right')
        plt.show()
    except:
        print('error reading OUTCAR')
        sys.exit(1)

if __name__=='__main__':
    inputfile='./OUTCAR'
    try:
        opts,args=getopt.getopt(sys.argv[1:],'hi:',['help','input='])
    except getopt.GetoptError:
        print('error in command line syntax')
        sys.exit(2)
    for i,j in opts:
        if i in ['-h','--help']:
            print('input options:\n\t-i, --input\t\tspecify an input file name other than ''OUTCAR''\n\nhelp options:\n\t-h, --help\t\tdisplay this help message')
            sys.exit()
        elif i in ['-i','--input']:
            inputfile=j
    main(inputfile)