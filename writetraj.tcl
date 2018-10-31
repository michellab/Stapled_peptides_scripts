set all [atomselect 0 "all"]
animate write crd traj.crd beg 0 end -1 skip 1 waitfor all sel $all 0
quit

