#!/usr/bin/python
############
# Name:
# SetupREMD-VERSION.py
#
###########
############
# Author:
# Marie BLUNTZER
#
###########
# Version :
# 0.4 (developpment)
#
###########
# Description:
# This script setup and makes a 100 ns test-run Gromacs simulation
# using AMBER99SB-ILDN forcefield
# parametrized TFE molecule solvent if required
# and TIP4P-Ewald water molecule type
#
#
###########


import numpy as np
import sys
import subprocess
import os
import argparse
import shutil
from shutil import copyfile
import fileinput
import re
import multiprocessing
from multiprocessing import Process ,Queue ,Pool
from collections import Counter
import string
from math import *


templdir='/home/marie/Stapled_peptide_git/'

def check_arguments() :
 if args.pdbfile==True and args.pdbfile[-3:]!='pdb' and os.path.isfile(args.pdbfile)!=True:
   sys.exit('specify a valid input file')
 if os.path.exists(templdir)==False:
   sys.exit('check your templates folder')


def get_nTFE():
    volfile=open('volume','r')
    Vfree=float(volfile.readline().split(':')[1].split()[0])*0.99
    if args.conc==True:
        print args.conc
        # nTFE=C*Vfree*avogadro
        # hyp1 : linear small peptide --> volume occupied max 2 % Vfree=Vbox*0.98--> wont work for a normal protein ..

        return ceil(args.conc*Vfree*6.002)
    if args.frac!=0 :
        return  ceil(args.frac*Vfree*6.022*1.39/100.02)
    else : return False


def setupsimulation() :
 if not os.path.exists('inputs'):
   os.mkdir('inputs')
#shutil.copyfile(args.pdbfile,'inputs/peptide.pdb' )
 os.chdir('inputs')
 log=open('gmxcmd.bsh', 'w')
 aa={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN', \
     'P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR','s':'GLS',  'r' : 'GLR'  }
 sequence = 'ACE'

 indices_stapled_residues = [i for i, x in enumerate(args.seq) if x == "s" or x == "r" ]
 bond=''
 if args.stapled :
    lib= 'loadoff   %sRES_lib/GLSStapled.lib \n loadoff  %sRES_lib/GLRStapled.lib  \n'  %(templdir ,templdir)
    while indices_stapled_residues : bond= 'bond m.%s.CE  m.%s.CE' %( indices_stapled_residues.pop()+2 ,indices_stapled_residues.pop() +2)
 else  : lib= 'loadoff   %sRES_lib/GLS.lib \n loadoff  %sRES_lib/GLR.lib  \n'  %(templdir ,templdir)



 for a in args.seq[:-1]:
  sequence=' '.join([sequence , aa[a] ])
 Cterminal = 'C' + aa[args.seq[-1].upper()]
 sequence=' '.join([sequence ,Cterminal ])
 #impose an helical start( up to 20 residues )
 helical=''
 if args.helical : helical = 'impose m { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 } { { "N" "CA" "C" "N" -40.0 } { "C" "N" "CA" "C" -60.0 } }'

 tleap=open('tleap.in' , 'w')
 tleap.writelines(['source leaprc.protein.ff14SB \n' , \
                  'loadamberparams frcmod.ionsjc_tip4pew \n' ,\
                  'loadamberparams %sRES_lib/parameters-atoms-bonds-angles.frcmod \n' %(templdir), \
                  'loadamberparams %sRES_lib/parameters-dihedrals-gaff2.frcmod \n' %(templdir), \
                  '%s\n'  %(lib) ,\
                  'set default PBradii mbondi2  \n'  ,\
                  'source leaprc.water.tip4pew \n' ,\
                  'm = sequence { %s } \n'%(sequence) ,\
                  '%s \n' %(helical) ,\
		  '%s \n' %(bond) ,\
                  'relax m  \n' , \
                  'savepdb m peptide.pdb \n',\
                  'charge m \n' ,\
                  '# solvatebox m prot TIP4PBOX 10.0 iso \n'  , \
                  '# addIonsRand prot Na+ 0 Cl- 0 \n' ,\
                  'saveamberparm m peptide.prmtop peptide.inpcrd \n' , \
                  'quit' ] )
 tleap.close()
 subprocess.call('tleap -f tleap.in' , shell=True)!=0
 subprocess.call('acpype -x peptide.inpcrd -p peptide.prmtop '  , shell=True )

 os.chdir('..')

 if not os.path.exists('simprep'):
  os.mkdir('simprep')
 shutil.copyfile('%s/GROMACS_SIM/mini-peptide.mdp' %(templdir), 'simprep/prep.mdp')

 os.chdir('simprep')
 shutil.copyfile('../inputs/peptide_GMX.top', 'topol.top')
 ## 1 --> TIP3P, 2 --> TIP4P, 3 --> TIP4P-Ewald
 #cmd='echo 3 | pdb2gmx  -quiet -f ../inputs/peptide.pdb -ignh -o peptide.gro -ff amber99sb-ildn '
 #log.write(cmd+'\n')
 #if subprocess.check_output(cmd)!=0 : sys.exit( ' ######## \n !!!  ERROR !!! \n ########## \n DEBUG : try to run : \n cd %s \n %s' %(os.getcwd(), cmd ))
 namegro='peptide_GMX'
 cmd='gmx editconf  -quiet -f ../inputs/%s.gro -o %s+box.gro -bt octahedron  -box 5.13714   4.84335   4.19446  | grep "new box volume" > volume' %(namegro , namegro)
 log.write(cmd+'\n')
 subprocess.check_output(cmd  , shell=True)
 namegro+='+box'


 subprocess.check_output('mv topol.top topol_old.top' , shell=True  )
 topolfile=open('topol_old.top','r' )
 lines=topolfile.readlines()
 topolfile.close()
 topolfile=open('topol.top','a' )
 i=0
 print((len(lines)))

 while i<len(lines) and 'moleculetype' not in lines[i]:
           topolfile.write(lines[i])
           i+=1
 topolfile.write(';   Include forcefield parameters \n')
 topolfile.write('#include "/home/marie/Stapled_peptide_git/TFE_lib/fffield.itp" \n \n' )

 if args.conc!=False or args.frac!=0:
    #subprocess.check_output('genbox   -cp %s.gro -ci %s/TFE.acpype/TFE_GMX.gro -nmol %s -o %s+TFE%s.gro >> /dev\null' %(namegro,templdir,conc,namegro,str(conc)))
    nTFE=int(get_nTFE())

    cmd='gmx insert-molecules  -quiet  -f %s.gro -ci %s/TFE_lib/TFE_GMX.gro -nmol %s -o %s+TFE.gro ' %(namegro,templdir,nTFE,namegro)
    log.write(cmd+'\n')
    subprocess.check_output(cmd + ' &> /dev/null'  , shell=True )
    namegro+='+TFE'
#    topolfile.write('; Include TFE parameters \n')
#    topolfile.write('#include "%sTFE_lib/TFE_GMX_B15.itp" \n \n' %(templdir) )


# while i<len(lines) and 'system' not in lines[i]:
#         topolfile.write(lines[i])
#         i+=1
# topolfile.write('; Include Water Parameters \n')
# topolfile.write('#include "/usr/local/gromacs/share/gromacs/top/amber99sb-ildn.ff/tip4pew.itp" \n \n')
#    topolfile.write('; Include Ions Parameters \n')
#    topolfile.write('#include "/usr/local/gromacs/share/gromacs/top/amber99sb-ildn.ff/ions.itp" \n \n')

    # topolfile.write('#include "%sTFE.acpype/TFE_GMX.itp" \n' %(templdir) )
 topolfile.writelines(lines[i:])
 topolfile.flush()
 if args.conc!=False or args.frac!=0: topolfile.write('TFE	           %s\n' %(str(nTFE)))
 topolfile.flush()
 topolfile.close

 cmd='gmx solvate  -quiet -cp %s.gro -cs tip4p.gro -p topol.top -o %s+solv.gro'%(namegro,namegro)
 log.write(cmd+'\n')

 subprocess.check_output(cmd + ' &> /dev/null' , shell=True )
 namegro+='+solv'


 #print 'gmx grompp -f prep.mdp -c %s.gro -p topol.top -o topol.tpr'%(namegro)
 cmd='gmx grompp  -quiet -f prep.mdp -c %s.gro -p topol.top -o topol.tpr -maxwarn 10 ' %(namegro)
 log.write(cmd + '\n')
 subprocess.check_output(cmd  , shell=True)!=0
 namegro+='+ions'
 cmd='echo sol | gmx genion  -quiet -s topol.tpr  -o %s.gro -p topol.top -pname NA -nname CL -neutral '%(namegro)
 log.write(cmd+'\n')
 subprocess.check_output(cmd + ' &> /dev/null'   , shell=True)!=0
 shutil.copyfile('%s.gro' %(namegro), 'conf-in.gro')
 os.chdir(workdir)


def minimization():
  if not os.path.exists('minimization'):
   os.mkdir('minimization')
  shutil.copyfile('%s/GROMACS_SIM/mini-peptide.mdp' %(templdir), 'minimization/min.mdp')
  os.chdir('minimization')
  cmd = 'gmx grompp  -f min.mdp -c ../simprep/conf-in.gro -p ../simprep/topol.top -o topol.tpr -maxwarn 5 '
  subprocess.check_output(cmd + ' &> /dev/null'  , shell=True )
  subprocess.check_output('gmx mdrun -gpu_id 0 -nt 3', shell=True)
  os.chdir(workdir)



def equilibration(T=300):
     if not os.path.exists('equilibration'):
          os.mkdir('equilibration')
     direquil2='equilibration/equilibration-'+str(T)
     if not os.path.exists(direquil2):
           os.mkdir(direquil2)

     conffile=open('%s//GROMACS_SIM/equi-peptide.mdp' %(templdir),'r')
     lines=conffile.readlines()
     if len(lines)==0:
             sys.exit()
     conffile.close()
     conffile_new=open('%s/equi.mdp' %(direquil2),'w')
     for line in lines:
           conffile_new.write(line.replace('$TEMP',str(T)))
     conffile_new.close()
     os.chdir(direquil2)
     cmd='gmx grompp -maxwarn 1 -quiet -f equi.mdp -c ../../minimization/confout.gro -p ../../simprep/topol.top -o topol.tpr'
     subprocess.check_output(cmd + ' &> /dev/null ', shell=True)
     subprocess.check_output('gmx mdrun -quiet -gpu_id 0 -nt 3  ', shell=True )
     os.chdir(workdir)
     conffile.close()

def equilibration2(T=300):
 direquil2='equilibration2/equilibration2-'+str(T)
 if not os.path.exists('equilibration2'):
   os.mkdir('equilibration2')
 if not os.path.exists(direquil2):
  os.mkdir(direquil2)
 conffile=open('%s//GROMACS_SIM/equi2-peptide.mdp' %(templdir),'r')
 lines=conffile.readlines()
 if len(lines)==0:
     sys.exit()
 conffile.close()
 conffile_new=open('%s/equi2.mdp' %(direquil2),'w')
 for line in lines:
   conffile_new.write(line.replace('$TEMP',str(T)))
 conffile_new.close()
 os.chdir(direquil2)
 subprocess.check_output('gmx grompp  -quiet -f equi2.mdp -c ../../equilibration/equilibration-%s/confout.gro -p ../../simprep/topol.top -o topol.tpr -maxwarn 5' %(str(T)), shell=True)
 subprocess.check_output('gmx mdrun  -gpu_id 1 -nt 3 ' , shell=True)
 os.chdir(workdir)
 conffile.close()


def REMDtempcalc():
### Setup a linear temperature distribution according to how many cores are available
### Thi has to be changed
#temperatures=np.arange(Tmin,Tmax,int(ceil(((Tmax-Tmin)/(n_replicas-1)))))

    #New Set of Temperaratures generated by http://folding.bmc.uu.se/remd/tgenerator.php
    # 'aranged' (by changing the number of atoms) to fit the 12 cores required by ARCHER

    return temperatures[:108]


def prodREMDsetup(i,T):

     proddir='prodREMD/'
     filesetup='prod-'+str(i)
     if not os.path.exists('prodREMD'):
      os.mkdir('prodREMD')

     filemdptemplate=open('%s/GROMACS_SIM/prodTEMP-peptide.mdp' %(templdir),'r')
     filemdp=open('prodREMD/prod_%s.mdp' %('{:03d}'.format(i)), 'w')
     lines=filemdptemplate.readlines()
     filemdptemplate.close()
     for line in lines:

      filemdp.write(line.replace('$TEMP',str(T)))
     filemdp.close()

     shutil.copyfile('./equilibration2/equilibration2-%s/confout.gro' %(T) ,  './prodREMD/conf_%s.gro' %('{:03d}'.format(i)))

############################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='# This script setup an REMD Gromacs simulation \
    \n  using a concentration of TFE (default pure water) \
    \n# using AMBER99SB-ILDN forcefield \
    \n # and TIP4P-Ewald water molecule type')
    parser.add_argument('--conc', type=float,  default=False,  help=' concentration of TFA in M (Seems that maximum is a little less than 1 M, more investigation to come...)')
    parser.add_argument('--frac', default=0.25 , type=float,   help=' concentration of TFA in M (Seems that maximum is a little less than 1 M, more investigation to come...)')
    parser.add_argument('--pdbfile', default='peptide.pdb' , type=str , help='input PDB file ')
    parser.add_argument('--seq'  , required=True, type=str , help='input sequence')
    parser.add_argument('--start_from' ,  default='prep' , help='step to restart from (default system preparation values prep, min, equ1,  equ2, prodtest  )' )
    #parser.add_argument('--prodtest' , default=True )
    parser.add_argument('--helical' , action='store_true' )
    parser.add_argument('--noREMD' , action='store_true'  )
    parser.add_argument('--stapled' , action='store_true'  )
    #parser.add_argument('--n_cores', type=int , default=multiprocessing.cpu_count() , help='number of cores to use, default maximum available'  )
    args = parser.parse_args()

    steps={'prep':1, 'min':2, 'equ1':3,  'equ2':4, 'prod':5}

    workdir=os.getcwd()
    if args.noREMD==True : REMD = False
    else : REMD=True
    if steps.get(args.start_from)<2:
        check_arguments()
        setupsimulation()
        print('run minimization')
    if steps.get(args.start_from)<3:
        minimization()
        equilibration(300)
        equilibration2(300)
        shutil.copyfile( 'simprep/topol.top' , 'prodREMD/topol.top' )
'''
    if steps.get(args.start_from)<4:


        temperatures=REMDtempcalc()
        print(temperatures)

        for T in temperatures:
              print( 'equilibration2 at', T)
              equilibration(T)





    if steps.get(args.start_from)<5 and REMD==0:
        equilibration2()



    if steps.get(args.start_from)<6 and REMD==1:

     temperatures=REMDtempcalc()
     print(temperatures)

     for T in temperatures:
      print( 'equilibration2 at', T)

      equilibration2(T)
     '
     ### paralezitian 'works', but also makes gromacs crashe even with 2 processes...'
     pool = Pool(processes=2)
     print (pool.map(equilibration2,temperatures) )


     #### Other way same issue
     processes = [multiprocessing.Process(target=equilibration2, args=temperatures)]
     # Run processes
     for p in processes:
         p.start()
     # Exit the completed processes
     for p in processes:
         p.join()
     '''


    # prodREMDrun(temperatures)
