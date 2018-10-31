import os
import shutil
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

for i in range(len(temperatures)):
    print(temperatures[i])
    prodREMDsetup(i,temperatures[i])
# subprocess.check_output('gmx grompp  -quiet -f prodREMD.mdp -c %s/equilibration2/equilibration2-%s/confout.gro -p %s/simprep/topol.top -o topol.tpr -maxwarn 5' %(workdir,T,workdir ))



templdir='/home/marie/Stapled_peptide_git/'


temperatures=[300.00, 302.14, 304.30, 306.48, 308.67, 310.88, 313.11, 315.35, 317.60, 319.88, 322.17, 324.48, 326.80, 329.14, 331.50, \
     333.88, 336.27, 338.68, 341.11, 343.56, 346.02, 348.50, 351.00, 353.52, 356.06, 358.62, 361.19, 363.77, 366.38, 369.01, 371.67, 374.34,\
      377.03, 379.74, 382.47, 385.22, 387.98, 390.78, 393.59, 396.42, 399.27, 402.14, 405.03, 407.94, 410.88, 413.83, 416.81, 419.81, 422.83 \
      , 425.88, 428.94, 432.03, 435.15, 438.28, 441.43, 444.61, 447.81, 451.04, 454.28, 457.56, 460.85, 464.17, 467.50, 470.87, 474.26, \
      477.68, 481.12, 484.58, 488.07, 491.59, 495.13, 498.70, 502.28, 505.90, 509.54, 513.22, 516.92, 520.65, 524.40, 528.18, 531.99, 535.82\
      , 539.69, 543.57, 547.48, 551.43, 555.41, 559.41, 563.44, 567.51, 571.60, 575.72, 579.87, 584.05, 588.26, 592.50, 596.77, 601.08, \
      605.41, 609.77, 614.17, 618.60, 623.06, 627.55, 632.07, 636.62, 641.21, 645.83, 650.49, 655.21, 659.93, 664.69, 669.48, 674.30, 679.17\
      , 684.07, 689.00, 693.97, 699.97, 705.01]

for i in range(len(temperatures)):
    print(temperatures[i])
    prodREMDsetup(i,temperatures[i])
