jobs -l import os
import sys
import mdtraj
import numpy as np
from  itertools import *
import subprocess
import argparse
from  matplotlib import pyplot as  plt
import pandas as pd
import seaborn as sns



def process_traj(topol, traj):

    print( '.... processing traj .... ')

    process_traj='\
    echo 1 1 | gmx trjconv -f %s -pbc cluster -o traj_cluster.xtc -s md.tpr ; \n \
    echo 1 1 | gmx trjconv -f traj_cluster.xtc  -o traj_cluster_alg.xtc -s md.tpr -fit rot+trans  ;  \n \
    echo 1 1 | gmx trjconv -f %s -pbc cluster -o traj_cluster.pdb -s md.tpr ; \n'  %(traj, topol)

    process = subprocess.check_output(process_traj,shell=True )
    print ('.... processing traj .... DONE')

def helicity(traj,peptide_chain):
     print( '.... computing helicity ....')
     trajpep= traj.atom_slice(  traj.topology.select(peptide_chain))
     dssp=mdtraj.compute_dssp(trajpep, simplified=True)
#     residues=[residue.index for residue in traj.topology.chain(peptide_chain).residues ]
     residues=peptide_chain
#     unique, counts = np.unique(dssp[:,residues[0]:residues[-1]], return_counts=True)
     unique, counts = np.unique(dssp[:,:], return_counts=True)
     hel = dict(zip(unique, counts)).get('H',0)/np.sum(counts)*10
     print('helicity ' + str(hel))

     print( '.... computing helicity .... DONE')
     return hel

def find_centroid(traj):
    print( '.... finding centroid ....')
    atom_indices = [a.index for a in traj.topology.atoms if a.name == 'CA']
    distances = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
      distances[i] = mdtraj.rmsd(traj, traj, i, atom_indices=atom_indices)
    beta = 1
    index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()
    print( '.... finding centroid .... DONE')
    return index

def compute_rmsd(traj, peptide_chain, centroid):
    print( '.... computing RMSD ....')
#  atoms_indices = traj.topology.chain("%s" %(peptide_chain))
#    atoms_indices = traj.topology.select("chainid  %s" %(peptide_chain))
    atoms_indices = traj.topology.select("%s" %(peptide_chain))
    rmsd = mdtraj.rmsd(traj, traj, frame=centroid, atom_indices=atoms_indices, parallel=True, precentered=True)

    print( '.... computing RMSD .... DONE')

    return rmsd

def compute_rmsd_sidechains(traj, peptide_chain, centroid):
    print( '.... computing RMSD side chains ....')
#  atoms_indices = traj.topology.chain("%s" %(peptide_chain))
    #atoms_indices = traj.topology.select("chainid  3 to 5 and sidechain" )
#    atoms_indices = traj.topology.select("chainid  3 and sidechain" )
    atoms_indices = traj.topology.select("%s and sidechain" %(peptide_chain))
    rmsd = mdtraj.rmsd(traj, traj, frame=centroid, atom_indices=atoms_indices, parallel=True, precentered=True)

    print( '.... computing RMSD sidechains.... DONE')

    return rmsd
def contacts_bonds(traj, peptide_chain ):
    print ('... computing peptide protein contacts ... ')
#  print(traj.topology.select( ' %s ' %(peptide_chain) ))
#  print(traj.topology.residues.index(0-10))
#  print(traj.topology.residue( ' %s ' %(peptide_chain)))
#  group_1 = [residue.index for residue in traj.topology.chain(peptide_chain).residues ]
#  group_1 = [residue.index for residue in  traj.topology.select( ' %s ' %(peptide_chain) ).residues ]
    group_1 = [residue.index for residue in traj.topology.residues if residue.index in peptide_chain ]
    group_2 = [residue.index for residue in traj.topology.residues if residue.index in chainAlist or residue.index in  chainBlist or residue.index in chainClist ]
  #group_2 = [residue.index for residue in   traj.topology.chain(0).residues or traj.topology.chain(1).residues or   traj.topology.chain(2).residues ]
 # group_2 = [residue.index for residue in   traj.topology.select( ' %s or %s or %s' %(chainA, chainB,chainC ) ).residues   ]
    pairs = list(product(group_1, group_2))

    contacts_bonds= mdtraj.compute_contacts(traj,pairs , scheme='closest-heavy', ignore_nonprotein=True, periodic=True, soft_min=False, soft_min_beta=20)


def hydrogen_bonds(traj):
  print( '.... computing hbonds (can take few minutes) ....')

  trajwtwater=traj.atom_slice(  traj.topology.select('protein'))
#  hbonds = mdtraj.baker_hubbard(traj, freq =0.4)
  hbonds = mdtraj.baker_hubbard(trajwtwater, freq =0.4)
  return hbonds

def get_peptides_hbonds(hbonds ,peptide_chain ) :
  list_peptide_hbonds=[]
  for hbond in hbonds:
      #if ( hbond[0] in  traj.topology.select( 'chainid  %s to %s' % (3,5) ) and   hbond[2] not in  traj.topology.select( 'chainid  %s to %s' % (3,5)) ) :
#      if ( hbond[0] in  traj.topology.select( 'chainid  %s ' % (3) ) and   hbond[2] not in  traj.topology.select( 'chainid  %s' % (3)) ) :
        if ( hbond[0] in  traj.topology.select( '%s' %(peptide_chain)) and   hbond[2] not in  traj.topology.select( '%s' %(peptide_chain) ) ) :
          list_peptide_hbonds.append( hbond[[0,2]] )
          print(' hbond : %s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])))
          print( 'hbond : %s -- %s' % (hbond[0], hbond[2]))
      #elif ( hbond[0] not in  traj.topology.select( 'chainid  %s to %s' % (3,5) ) and   hbond[2]  in  traj.topology.select( 'chainid  %s to %s' % (3,5)) ) :
#      elif ( hbond[0] not in  traj.topology.select( 'chainid  %s ' % (3) ) and   hbond[2]  in  traj.topology.select( 'chainid  %s ' % (3)) ) :
        elif ( hbond[0] not in  traj.topology.select( ' %s '%(peptide_chain)) and   hbond[2]  in  traj.topology.select( '%s ' %(peptide_chain) ) ) :
          list_peptide_hbonds.append( hbond[[0,2]] )
          print(' hbond : %s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])))
          print('hbond : %s -- %s' % (hbond[0], hbond[2]))
  '''
  da_distances = mdtraj.compute_distances(traj,  np.array(list_peptide_hbonds), periodic=False)

  color = cycle(['r', 'b', 'gold'])
  print(len(np.array(list_peptide_hbonds)))
  for i in  range(len(np.array(list_peptide_hbonds))):
      plt.hist(da_distances[:, i], color=next(color), label=label(hbonds[i]), alpha=0.5)
  plt.legend()
  plt.ylabel('Freq');
  plt.xlabel('Donor-acceptor distance [nm]')
  plt.savefig('hbonds.png', bbox_inches='tight')
  '''
  return np.array(list_peptide_hbonds)

def get_specific_peptides_hbonds(hbonds ,peptide_chain ) :
  list_peptide_hbonds=[]
  for hbond in hbonds:
      #if ( hbond[0] in  traj.topology.select( 'chainid  %s to %s' % (3,5) ) and   hbond[2] not in  traj.topology.select( 'chainid  %s to %s' % (3,5)) ) :
#      if ( hbond[0] in  traj.topology.select( 'chainid  %s ' % (3) ) and   hbond[2] not in  traj.topology.select( 'chainid  %s' % (3)) ) :
        if ( hbond[0] in  traj.topology.select( '%s and not (resname  ACE or resname   ALA)  ' %(peptide_chain)) and   hbond[2] not in  traj.topology.select( '%s and not (resname  ACE or resname   ALA)  ' %(peptide_chain) ) ) :
          list_peptide_hbonds.append( hbond[[0,2]] )
          print(' hbond : %s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])))
          print( 'hbond : %s -- %s' % (hbond[0], hbond[2]))
      #elif ( hbond[0] not in  traj.topology.select( 'chainid  %s to %s' % (3,5) ) and   hbond[2]  in  traj.topology.select( 'chainid  %s to %s' % (3,5)) ) :
#      elif ( hbond[0] not in  traj.topology.select( 'chainid  %s ' % (3) ) and   hbond[2]  in  traj.topology.select( 'chainid  %s ' % (3)) ) :
        elif ( hbond[0] not in  traj.topology.select(  '%s and not (resname  ACE or resname   ALA)  '  %(peptide_chain)) and   hbond[2]  in  traj.topology.select( '%s and not (resname  ACE or resname   ALA)  ' %(peptide_chain) ) ) :
          list_peptide_hbonds.append( hbond[[0,2]] )
          print(' hbond : %s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])))
          print('hbond : %s -- %s' % (hbond[0], hbond[2]))
  return np.array(list_peptide_hbonds)




def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

#------------------------------------------------------------
'''
def draw_radar(results):
    N = len(results)
#    theta = radar_factory(N, frame='polygon')
    r= np.array(results)
    theta = np.arange(0, 2 * np.pi , 2 * np.pi /len(results) ,dtype=float )
    area = 200
    fig, axes = plt.subplots(figsize=(N, N), nrows=2, ncols=2,
                             subplot_kw=dict(projection='radar'))
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.set_xticklabels(['helicity', 'rmsd' , 'hbond'])
    c = ax.scatter(theta, r, c=colors, s=area, cmap='hsv', alpha=0.75)
    fig.savefig('score.png', bbox_inches='tight')
'''
def draw_radar(resultsA, resultsB , resultsC , labels):

    fig=plt.figure('results.png')
    ax = fig.add_subplot(111, polar=True)
    angles=np.linspace(0, 2*np.pi, len(labels), endpoint=False)
    ax.set_thetagrids(angles * 180/np.pi, labels)
    angles= np.concatenate( (  angles ,[angles[0]]))
    for results in [resultsA, resultsB , resultsC]:
# close the plot
        stats=results
        stats=np.concatenate((stats,[stats[0]]))
    #    angles= np.concatenate( (  [angles[0]]*3 ,[angles[1]]*3 , [angles[2]]*3 ,[angles[3]]*3 ,[angles[4]]*3 ,[angles[0]]))


        ax.plot(angles, stats, 'o-', linewidth=2)
        ax.fill(angles, stats, alpha=0.25)

    ax.set_title('Scoring ')
    ax.grid(True)
    fig.savefig('results.png', bbox_inches='tight')
###


parser = argparse.ArgumentParser(description='Process MDtraj')
parser.add_argument('--peptide_chains', metavar='N', type=str, nargs='+',
                    help='peptide_chains')
parser.add_argument('--topol', default='topol.prmtop', type=str,
                    help='topology file (.gro or .pdb)')
parser.add_argument('--traj', default='traj.nc', type=str,
                    help='md.xtc')
parser.add_argument('--process_traj',  type=str2bool, nargs='?',
                        const=True, default=True, help='cluster and align protein')




args = parser.parse_args()

chainA= ' resid 1 to 152'
chainB= ' resid 153 to 304'
chainC= ' resid 305 to 456'
chainD= ' resid 458 to 468'
chainE= ' resid 470 to 480'
chainF= ' resid 482 to 492'
chainAlist = range ( 0  ,151)
chainBlist = range (152 , 303)
chainClist = range (304 , 455)
chainDlist = range (457 , 467)
chainElist = range (469 , 479)
chainFlist = range (481 , 491)




'''
if  args.process_traj ==True:
    process_traj(args.topol, args.traj )
    traj = mdtraj.load('traj_cluster_alg.xtc', top='traj_cluster.pdb')
'''

print('Loading trajectory' )
traj = mdtraj.load(args.traj, top=args.topol)
print('centering coodinates' )
traj.center_coordinates()
results=[]
'''  ... 1
for i in [chainDlist , chainElist , chainFlist ] :
    contact = contacts_bonds(traj, i )
'''
for i in [chainD , chainE , chainF] :
    hel = helicity(traj,i)
    results.append(hel)



plt.figure()
#centroid=find_centroid(traj)
centroid=0
for i in [chainD , chainE , chainF]:
        rmsd=compute_rmsd(traj, i, centroid)
        results.append(1/np.average(rmsd))
        plt.plot(range(traj.n_frames), rmsd, '-')


plt.title('rmsd')
plt.xlabel('time')
plt.ylabel('rmsd')
plt.savefig('rmsd.png', bbox_inches='tight')
plt.figure()

for i in [chainD , chainE , chainF]:
    side = compute_rmsd_sidechains(traj, i , centroid)
    results.append((1/np.average(side)))


hydrogenBonds= hydrogen_bonds(traj)
for i in [chainD , chainE , chainF]:
    hbs= get_peptides_hbonds(hydrogenBonds,i)
    results.append(hbs.size)



for i in [chainD , chainE , chainF]:
    hbs= get_specific_peptides_hbonds(hydrogenBonds,i)
    results.append(hbs.size)






labels=np.array(['helicity','rmsd', 'sidechain rmsd',  'hydrogen bonds' , 'specific hydrogen bonds' ])

draw_radar([res for res in  results if res % 3 == 0 ] ,[res for res in  results if res % 3 == 1 ] , [res for res in  results if res % 3 == 2 ] , labels)
