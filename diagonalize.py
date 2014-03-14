
from Cluster import *

types = [ "A", "B", "C", "D" ]

# single sites
if ORDER >= 0:
  site_code = []
  for t in types:
    site_code.append("s"+t)
  print site_code
  print

  siteClus = []
  siteEval = []
  site_p = []
  site_w = []
  
  for c in range(4):
    siteClus.append( Cluster("site",1,c,site_code[c]) )
    siteEval.append( siteClus[c].Diagonalize(site_code[c],wf_switch) )

    site_p.append( siteClus[c].PropertyCSM(site_code[c],field_mag,siteEval[c]) )
    site_w.append( site_p[c] )

# single tetrahedrons
if ORDER >= 1:
  t1_code = ["t1X"]
  print t1_code
  print

  t1Clus = []
  t1Eval = []
  t1_p = []
  t1_w = []

  #print "{0:} is {1:}".format(0,t1_code[0])
  t1Clus.append( Cluster("tetra",1,None,t1_code[0]) )
  t1Eval.append( t1Clus[0].Diagonalize(t1_code[0],wf_switch) )

  t1_p.append( t1Clus[0].PropertyCSM(t1_code[0],field_mag,t1Eval[0]) )
  t1_w.append( t1Clus[0].Weight(t1_p,t1_code[0],site_w,None,None,None,t1_code,None,None,None,None) )

# two tetrahedrons
if ORDER >= 2:
  t2_code = []
  for c in types: 
    t2_code.append("t2"+c)
  print t2_code
  print

  t2Clus = []
  t2Eval = []
  t2_p = []
  t2_w = []

  for c in range(4):
    #print "{0:} is {1:}".format(c,t2_code[c])
    link2 = [ (0,1,c) ]
    t2Clus.append( Cluster("tetra",2,link2,t2_code[c]) )
    t2Eval.append( t2Clus[c].Diagonalize(t2_code[c],wf_switch) )

    t2_p.append( t2Clus[c].PropertyCSM(t2_code[c],field_mag,t2Eval[c]) )
    t2_w.append( t2Clus[c].Weight(t2_p,t2_code[c],site_w,t1_w,None,None,t1_code,t2_code,None,None,None) )

# three tetrahedrons
if ORDER >= 3:
  t3_code = []

  for i in range(len(types)):
    for j in range(i):
      label = types[j]+types[i]
      t3_code.append( "t3"+label )
  print t3_code
  print

  t3Clus = []
  t3Eval = []
  t3_p = []
  t3_w = []

  for i in range(len(types)):
    for j in range(i):
      link3 = [ (0,1,j), (1,2,i) ]
      index = len(t3Clus)
      #print
      #print "{0:} is {1:}".format(index,t3_code[index])
      t3Clus.append( Cluster("tetra",3,link3,t3_code[index]) )
      t3Eval.append( t3Clus[index].Diagonalize(t3_code[index],wf_switch) )

      t3_p.append( t3Clus[index].PropertyCSM(t3_code[index],field_mag,t3Eval[index]) )
      t3_w.append( t3Clus[index].Weight(t3_p,t3_code[index],site_w,t1_w,t2_w,None,t1_code,t2_code,t3_code,None,None) )

# four tetrahedrons
if ORDER >= 4:
  t4s_code = []

  for i in range(4):
    for j in range(4):
      for k in range(4):
        if j != i and k != j and k >= i:
          label = types[i] + types[j] + types[k]
          t4s_code.append( "t4s"+label )
  print t4s_code
  print
  t4y_code = ["t4yABC", "t4yABD", "t4yACD", "t4yBCD"]
  print t4y_code
  print


  t4sClus = []
  t4sEval = []
  t4s_p = []
  t4s_w = []

  for i in range(4):
    for j in range(4):
      for k in range(4):
        if j != i and k != j and k >= i:
          link4s = [ (0,1,i), (1,2,j), (2,3,k) ]
          index = len(t4sClus)
          #print
          #print "{0:} is {1:}. {2:}{3:}{4:}".format(index,t4s_code[index],i,j,k)
          t4sClus.append( Cluster("tetra",4,link4s,t4s_code[index]) )
          t4sEval.append( t4sClus[index].Diagonalize(t4s_code[index],wf_switch) )

          t4s_p.append( t4sClus[index].PropertyCSM(t4s_code[index],field_mag,t4sEval[index]) )
          t4s_w.append( t4sClus[index].Weight(t4s_p,t4s_code[index],site_w,t1_w,t2_w,t3_w,t1_code,t2_code,t3_code,t4s_code,t4y_code) )

  t4yClus = []
  t4yEval = []
  t4y_p = []
  t4y_w = []

  link4yABC = [ (0,1,0), (1,2,1), (1,3,2) ]
  t4yClus.append( Cluster("tetra",4,link4yABC,"t4yABC") )
  t4yEval.append( t4yClus[0].Diagonalize("t4yABC",wf_switch) )

  link4yABD = [ (0,1,0), (1,2,1), (1,3,3) ]
  t4yClus.append( Cluster("tetra",4,link4yABD,"t4yABD") )
  t4yEval.append( t4yClus[1].Diagonalize("t4yABD",wf_switch) )

  link4yACD = [ (0,1,0), (1,2,2), (1,3,3) ]
  t4yClus.append( Cluster("tetra",4,link4yACD,"t4yACD") )
  t4yEval.append( t4yClus[2].Diagonalize("t4yACD",wf_switch) )

  link4yBCD = [ (0,1,1), (1,2,2), (1,3,3) ]
  t4yClus.append( Cluster("tetra",4,link4yBCD,"t4yBCD") )
  t4yEval.append( t4yClus[3].Diagonalize("t4yBCD",wf_switch) )

  index = 0
  for Cl in t4yClus:
    t4y_p.append( Cl.PropertyCSM(t4y_code[index],field_mag,t4yEval[index]) )
    t4y_w.append( Cl.Weight(t4y_p,t4y_code[index],site_w,t1_w,t2_w,t3_w,t1_code,t2_code,t3_code,t4s_code,t4y_code) )
    index+=1


Js = "J_{0:.3f}_{1:.3f}_{2:.3f}_{3:.3f}".format(Jzz,Jpm,Jpmpm,Jzpm)
filenameC = "C"+str(ORDER)+"_"+Js+"_h_"+str(field_mag)+"_d_"+dir_str
filenameS = "S"+str(ORDER)+"_"+Js+"_h_"+str(field_mag)+"_d_"+dir_str
filenameM = "M"+str(ORDER)+"_"+Js+"_h_"+str(field_mag)+"_d_"+dir_str
fileoutC = open(filenameC,"w")
fileoutS = open(filenameS,"w")
fileoutM = open(filenameM,"w")

#summation and printing
for t in range(T_POINTS):
  for p in range(3):
    si_val = 0
    t1_val = t1_w[0][p][t]
    t2_val = 0
    t3_val = 0
    t4_aba_val = 0
    t4_abc_val = 0
    t4_y_val = 0
    for c in range(4):
      si_val += .25*site_w[c][p][t]
      t2_val += .25*t2_w[c][p][t]
    for c in range( len(t3_w) ):
      t3_val += t3_w[c][p][t]/6.0  
    l0 = si_val
    l1 = .5*t1_val
    l2 = t2_val
    l3 = 3*t3_val
    o0 = l0
    o1 = l0 + l1
    o2 = l0 + l1 + l2
    o3 = l0 + l1 + l2 + l3
    e3 = l0 + l1 + l2 + .5*l3
    if ORDER >=4:
      for c in range( len(t4s_w) ):
        left = t4s_code[c][3]
        right = t4s_code[c][5]
        if left == right:
          t4_aba_val += t4s_w[c][p][t]/12.0 # ABA types
        else:
          t4_abc_val += t4s_w[c][p][t]/12.0 # ABC types ABC == CBA
      for c in range( len(t4y_w) ):
        t4_y_val += t4y_w[c][p][t]/4.0
      l4 = 3*t4_aba_val + 6*t4_abc_val + 2*t4_y_val
      o4 = l0 + l1 + l2 + l3 + l4
      e4 = l0 + l1 + l2 + .5*l3 + .25*(l3 + l4)
      if p == 0:
        outstringC = "{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:}".format(11.6045*T[t],o0,o1,o2,o3,o4,e3,e4)
        fileoutC.write( outstringC )
        fileoutC.write("\n")
      elif p == 1:
        outstringS = "{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:}".format(11.6045*T[t],o0,o1,o2,o3,o4,e3,e4)
        fileoutS.write( outstringS )
        fileoutS.write("\n")
      else:
        outstringM = "{0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:}".format(11.6045*T[t],h_meV*o0,h_meV*o1,h_meV*o2,h_meV*o3,h_meV*o4,h_meV*e3,h_meV*e4)
        fileoutM.write( outstringM )
        fileoutM.write("\n")
    #ORDER 4
    else:
      if p == 0:
        outstringC = "{0:} {1:} {2:} {3:} {4:} {5:}".format(11.6045*T[t],o0,o1,o2,o3,e3)
        fileoutC.write( outstringC )
        fileoutC.write("\n")
      elif p == 1:
        outstringS = "{0:} {1:} {2:} {3:} {4:} {5:}".format(11.6045*T[t],o0,o1,o2,o3,e3)
        fileoutS.write( outstringS )
        fileoutS.write("\n")
      else:
        outstringM = "{0:} {1:} {2:} {3:} {4:} {5:}".format(11.6045*T[t],h_meV*o0,h_meV*o1,h_meV*o2,h_meV*o3,h_meV*e3)
        fileoutM.write( outstringM )
        fileoutM.write("\n")
    #ORDER 3 only

fileoutC.close()
fileoutS.close()
fileoutM.close()















      
