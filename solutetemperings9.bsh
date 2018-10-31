#! /bin/bash


mkdir solutetemp
cd solutetemp
echo '' > plumed.dat
cp ../simprep/topol.top .
# six replicas
nrep=6
# "effective" temperature range
tmin=290
tmax=500

# build geometric progression
list=$(
awk -v n=$nrep \
    -v tmin=$tmin \
    -v tmax=$tmax \
  'BEGIN{for(i=0;i<n;i++){
    t=tmin*exp(i*log(tmax/tmin)/(n-1));
    printf(t); if(i<n-1)printf(",");
  }
}'
)

# clean directory
rm -fr \#*


for((i=0;i<nrep;i++))
do

# choose lambda as T[0]/T[i]
# remember that high temperature is equivalent to low lambda
  lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i+1))';}')
# process topology
# (if you are curious, try "diff topol0.top topol1.top" to see the changes)
  plumed partial_tempering $lambda < topol.top > topol$i.top
# prepare tpr file
# -maxwarn is often needed because box could be charged
done

cp /home/marie/Stapled_peptide_git/GROMACS_SIM/prod-300K-peptide.mdp prod.mdp
cp ../equilibration2/equilibration2-300/confout.gro   confout.gro  


echo '#! /bin/bash

#SBATCH --job-name=SoluteTemp
#SBATCH  --ntasks=12
#SBATCH --gres=gpu:6
#SBATCH --partition=GTX
#SBATCH --time 48:00:00


for((i=0;i<nrep;i++))
do
  gmx grompp -p  topol$i.top -o topol$i.tpr -c  confout.gro -f  prod.mdp -maxwarn 1 
done
           
mpirun -n 6 ~/gromacs-2016.5/bui/bin/mdrun  -s topol.tpr -multi 6  -plumed plumed.dat -hrex -replex 50000
' > runsolute.bsh



scp -r ../solutetemp marie@section9.chem.ed.ac.uk:temp



