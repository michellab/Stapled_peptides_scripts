



def mut_res(resid,X, Y, filein , fileout):

    f=open(filein,r)
    output=open(fileout,w)
    g=f.readlines()
    for line in g :
      newline=line
      t=True
      if line.split()[0]=='ATOM'
        if line.split()[5]==resid:
          if line.split()[2] in [N, CA, C ,O ,CB]:
            newline=line[:17]+X+line[20:]
          else : t= False
        if line.split()[5]==resid+4:
          if line.split()[2] in [N, CA, C ,O ,CB]:
            newline=line[:17]+Y+line[20:]
          else : t= False
      if t= True : output.writelines(newline)
      f.close()
      output.close()

if __name__=__main__


    for resid in [6:10]:
      for combination in [[AKS,AKR] , [AKS,AKS], [AKR,AKR], [AKR,AKS], [CXS,CXY] ] :
      mut_res(resid,combination[0], combination[1], filein , stapled-peptide-%s-%s-%s-%s.pdb (str(resid),str(resid+4),str(combination[0]),str(combination[1])))
