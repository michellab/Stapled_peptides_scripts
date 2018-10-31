#! /usr/bin/env python

###########################
##### This is an ugly script to generate conformations spaced from 10 degrees using tleap and
##### parse the generated mol2 files to gaussian input files.
##### please pick your atoms names from the lists at the end of this script and put the right residus name ( I said it was ugly ... )
##### This script would normally crach if anything is wrong but preview your mol2 files before running long gaussians calculations...
############################
import os
import subprocess
import argparse



def main(list_of_torsions):

    list_of_key={}

    tleap_input_file=open('tleap.in', 'w')
    inputfile=open(inputmol2, 'r')
    inputfiliecontent=inputfile.readlines()
    for l in range(len(inputfiliecontent)):
        if res_name in inputfiliecontent[l] :
            list_of_key.update({inputfiliecontent[l].split()[1]: inputfiliecontent[l].split()[0]})
    print (list_of_key)
    list_of_atoms=[]
    for torsion in list_of_torsions:
        list_of_atoms.append([list_of_key[torsion.replace("\"","").split()[0]],list_of_key[torsion.replace("\"","").split()[1]], \
                        list_of_key[torsion.replace("\"","").split()[2]],list_of_key[torsion.replace("\"","").split()[3]]])


    print(list_of_atoms)
    inputfile.close()

    tleap_input_file.writelines(['source leaprc.protein.ff14SB \n' , 'm = loadmol2 %s  \n' %(inputmol2) , 'relax m \n' ])
    print ('Shall we begin ?')
    for bond in range(len(list_of_torsions)) :
      for angle in range(-180,+180,10) :
        tleap_input_file.writelines( [ 'impose m  { 1 } {  { %s %s } }\n' %(list_of_torsions[bond] , angle) ,  'savemol2 m %s_%s.mol2 1\n' %(torsion_names[bond] , angle)] )

    tleap_input_file.writelines(['quit' ])
    tleap_input_file.close()
    print (' tleap -f tleap.in')
    subprocess.call('tleap -f tleap.in' ,shell=True  )
    for finp in os.listdir('./') :
             if finp[-4:]=='mol2' and finp[:3] in ['PHI', 'PSI' , 'CHI']:
                fi=open(finp, 'r')
                fo=open('%s' %(finp[:-4]+'gau'), 'w')
                fo.write('--Link1-- \n%%RWF=%s\n%%NoSave\n%%chk=%s \n#P b3lyp/6-31G* Opt=(ModRedundant,MaxCycle=400) Freq SCF=(Conver=6,MaxCycle=400) Pop=NoMBS\n \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],finp[:-5],charge,multiplicity))

                while 'ATOM' not in fi.readline()  : continue
                i=True
                while i==True :
                        line=fi.readline().split()
                        if 'TRIPOS' not  in line[0] :
                            #lineout=line[1][0]+'-' + line[5] + '-' + line[8]  + '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'

                            lineout=line[1][0]+ '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'
                            fo.write(lineout)
                        else:
                            i=False
                            fo.write('\n')
                numb=torsion_names.index(finp[0:4].replace('_',''))
                listdih=list_of_atoms[numb]

                #fo.write('F'+str(list_of_atoms[numb][0])+'\nF'+str(list_of_atoms[numb][1])+'\nF'+str(list_of_atoms[numb][2])+'\nF'+str(list_of_atoms[numb][3])+'\n')
                fo.write('D ' + str(listdih[0]) +' '+ str(listdih[1])+ ' '+ str(listdih[2])+ ' '+str(listdih[3])+ ' F\n')
                fo.close()
                fi.close


if __name__ =='__main__':

        PHI_ATOMS = '"C10" "N" "CA" "C"'
        PSI_ATOMS = '"N" "CA" "C" "N1"'

        CHI1_ATOMS = ['"N" "CA" "CB" "CC"',
                      '"N" "CA" "CB" "CG1"',
                      '"N" "CA" "CB" "S"',
                      '"N" "CA" "CB" "OG"',
                      '"N" "CA" "CB" "OG1"']

        CHI2_ATOMS = ['"CA" "CB" "CC" "CD"',
                      '"CA" "CB" "S" "CC"',
                      '"CA" "CB" "CG1" "CD1"',
                      '"CA" "CB" "CG" "OD1"',
                      '"CA" "CB" "CG" "ND1"',
                      '"CA" "CB" "CG" "SD"']

        CHI3_ATOMS = ['"CB" "S" "CC" "C"',
                      '"CB" "CC" "CD" "CE"',
                      '"CB" "CG" "CD" "OE1"',
                      '"CB" "CG" "SD" "CE"']

        CHI4_ATOMS = ['"CG" "CD" "NE" "CZ"',
                      '"S" "CC" "CD" "CE"',
		      '"CC" "CD" "CE" "NZ"']

        CHI5_ATOMS = '"CD" "CE" "NZ" "NY"'

      
        res_name=' NXS'
        list_of_torsions = [PHI_ATOMS, PSI_ATOMS,  CHI1_ATOMS[2], CHI2_ATOMS[1] ,CHI3_ATOMS[0] ,CHI4_ATOMS[1] ]
        torsion_names =[ 'PHI' , 'PSI' , 'CHI1', 'CHI2' , 'CHI3' ,'CHI4' ]
        '''
        res_name=' PGR'
        lis 2_-20.gau-runonphase9.bsh t_of_torsions = [PHI_ATOMS, PSI_ATOMS,  CHI1_ATOMS[0] ]
        torsion_names =[ 'PHI' , 'PSI' , 'CHI1']
        res_name=' LNR'
        list_of_torsions = [PHI_ATOMS, PSI_ATOMS,  CHI1_ATOMS[0], CHI2_ATOMS[0] ,CHI3_ATOMS[1] ,CHI4_ATOMS[2] , CHI5_ATOMS]
        torsion_names =[ 'PHI' , 'PSI' , 'CHI1', 'CHI2' , 'CHI3', 'CHI4',  'CHI5' ]
        res_name=' CLI'
        list_of_torsions = [PHI_ATOMS, PSI_ATOMS,  CHI1_ATOMS[0], CHI2_ATOMS[0] ,CHI3_ATOMS[1] ,CHI4_ATOMS[2] ]
        torsion_names =[ 'PHI' , 'PSI' , 'CHI1', 'CHI2' , 'CHI3', 'CHI4' ]
        res_name=' CPS'
        list_of_torsions = [PHI_ATOMS, PSI_ATOMS,  CHI1_ATOMS[0], CHI2_ATOMS[0] ,CHI3_ATOMS[1]  ]
        torsion_names =[ 'PHI' , 'PSI' , 'CHI1', 'CHI2' , 'CHI3' ]
	'''
        charge=0
        multiplicity=1
        inputmol2='Mol-sm_m1-c1.mol2'
        print(list_of_torsions)
        main(list_of_torsions)
