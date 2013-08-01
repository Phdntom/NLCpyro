
from constants import *
from class_lib import *
import numpy as np
#from collections import deque

class Cluster:

  def __init__(self,unit,num,links,file_ext):
    self.nodes = []
    self.edges = []
    self.sites = []
    self.bonds = []

    if unit == "site":
      self.nodes.append( Site(links) )
      self.sites.append( (self.nodes[0].site[0],links) )
      """
      tetra_map = []
      tetra_map.append( (0,0) )
      map_name = "tetra"+str(file_ext)
      map_file = open(map_name,"w")
      np.savez(map_file,tetra_map)
      map_file.close()
      """

    elif unit == "tetra":
      """
        nodes: tetrahedrons embedded onto a diamond lattice
          contains Tetra Objects
        edges: undirected paths between nodes in the diamond lattice
          contains 3tuples, (bond type, start node, end node)
        sites: points in cluster on the pyrochlore lattice
          conatins 2tuples, (site label, site type)
        bonds: undirected paths between sites in the pyrochlore lattice
          contains 2lists, [start site, end site]
      """
      print "start cluster process"
      tetra_map = []
      #build num tetrahedrons
      for n in range(num):
        self.nodes.append( Tetra(n) )

      #append the 1st one into the supercluster
      # n = 0
      for s in range(4):
        self.sites.append( (self.nodes[0].site[s],self.nodes[0].site[s]) )
        tetra_map.append( (s,s) )
      for b in range(6):
        self.bonds.append( self.nodes[0].bond[b] )

      #append the remaining via the articulation list, links
      for n in range(1,num):
        shift = 3*n+1
        maps = []
        self.edges.append( (links[n-1][2],links[n-1][0],links[n-1][1]) )
        for s in range(4):
          if s != links[n-1][2]:
            #finds the site of links, shifts the rest
            maps.append( (self.nodes[n].site[s],shift) )
            tetra_map.append( (self.nodes[n].site[s],shift) )

            self.sites.append( (shift,self.nodes[n].site[s]) )
            self.nodes[n].site[s] = shift
            shift += 1
          else:
            maps.append( (self.nodes[n].site[s], self.nodes[ links[n-1][0] ].site[s] ) )
            tetra_map.append( (self.nodes[n].site[s], self.nodes[ links[n-1][0] ].site[s] ) )
        for b in range(6):
          for s in range(4):
            if self.nodes[n].bond[b][0] == maps[s][0]:
              self.nodes[n].bond[b][0] = maps[s][1]
            if self.nodes[n].bond[b][1] == maps[s][0]:
              self.nodes[n].bond[b][1] = maps[s][1]
          self.bonds.append(self.nodes[n].bond[b])
        """
        print "map:"
        print maps
        """
      """
      map_name = "tetra"+str(file_ext)
      map_file = open(map_name,"w")
      np.savez(map_file,tetra_map)
      map_file.close()
      """
      """
      print "nodes"
      for i in range(len(self.nodes)):
        print self.nodes[i]
      print
      print "edges"
      print self.edges
      """
      print "tetra_map"
      print tetra_map

    #if it's a tetrahedron based cluster
    else:
      print 'Error'
    #else it's an error

    print "sites"
    print self.sites

    print "bonds"
    print self.bonds

  def BuildHamiltonian(self,h_field):
    bond_count = len(self.bonds)
    site_count = len(self.sites)

    N = 2**site_count
    H = np.zeros((N,N), dtype = complex)
    for k in range(N):
      for b in range(bond_count):
        s = self.bonds[b][0]
        t = self.bonds[b][1]
        s_ty = self.sites[s][1]
        t_ty = self.sites[t][1]
        spin_s = (k >> s) & 1
        spin_t = (k >> t) & 1

        H[k,k] += .25 * Jzz * (2*spin_s-1) * (2*spin_t-1)

        if spin_s ^ spin_t:
          l = k + (spin_s - spin_t) * 2**t + (spin_t - spin_s) * 2**s
          H[k,l] += -Jpm
        
        if spin_s == 0 and spin_t == 0:
          l = k + 2**s + 2**t
          H[k,l] += Jpmpm * gamma[s_ty,t_ty]
        if spin_s == 1 and spin_t == 1:
          l = k - 2**s - 2**t
          H[k,l] += Jpmpm * (-zeta[s_ty,t_ty])

        if spin_t == 0:
          l = k + 2**t
          H[k,l] += .5 * (2*spin_s-1) * Jzpm * zeta[s_ty,t_ty]
        if spin_t == 1:
          l = k - 2**t
          H[k,l] += .5 * (2*spin_s-1) * Jzpm * (-gamma[s_ty,t_ty])
        if spin_s == 0:
          l = k + 2**s
          H[k,l] += .5 * (2*spin_t-1) * Jzpm * zeta[t_ty,s_ty]
        if spin_s == 1:
          l = k - 2**s
          H[k,l] += .5 * (2*spin_t-1) * Jzpm * (-gamma[t_ty,s_ty])
      ### for b bonds

      for s in range(site_count):
        spin_s = (k >> s) & 1
        s_ty = self.sites[s][1]
        h = np.dot(R[s_ty],field_dir)
        H[k,k] -= h_field*h[2]*gz*.5*(2*spin_s-1)
        if spin_s == 0:
          l = k + 2**s
          H[k,l] -= h_field*.5*gxy*(h[0]-h[1]*1j)
        elif spin_s == 1:
          l = k - 2**s
          H[k,l] -= h_field*.5*gxy*(h[0]+h[1]*1j)
      ### for s sites

    ### for k states
    return H

  def Diagonalize(self,code_ext,wf_switch):

    field_list = [field_mag,field_mag+.0001]

    eigen_list = []
    N = 2**(len(self.sites))

    for h_ext in field_list:
      H = self.BuildHamiltonian(h_ext/h_meV)
      #print "hamiltonian constructed."
      #dir_str = str(field_dir[0])+str(field_dir[1])+str(field_dir[2])
      if wf_switch == True :
        eigen_w,eigen_v = np.linalg.eigh(H)
        #print "diagonalization completed."

        eigen_list.append( eigen_w )

        if print_eigen_switch:
          print eigen_w
          print eigen_v
          filename = "J_{0:.3f}_{1:.3f}_{2:.3f}_{3:.3f}_h_{4:.6f}_d_{5:}".format(Jzz,Jpm,Jpmpm,Jzpm,h_ext,dir_str)
          eigen_out = open(filename+"_"+code_ext+"_eval","w")
          eigenwf = open(filename+"_"+code_ext+"_wf","w")

          np.savez(eigen_out,eigen_w)
          np.savez(eigenwf,eigen_v)

          eigenwf.close()
          eigen_out.close()
        ##if print_eigen_switch == True
      ###if wf_switch == True

      else: #wf_switch == False
        eigen_w = np.linalg.eigvalsh(H)
        #print "diagonalization completed."
        #print eigen_w

        eigen_list.append( eigen_w )

        if print_eigen_switch:
          print eigen_w
          filename = "J_{0:.3f}_{1:.3f}_{2:.3f}_{3:.3f}_h_{4:.6f}_d_{5:}".format(Jzz,Jpm,Jpmpm,Jzpm,h_ext,dir_str)
          eigen_out = open(filename+"_"+code_ext+"_eval","w")
          np.savez(eigen_out,eigen_w)

          eigen_out.close()
        ##if print_eigen_switch
      ###else wf_switch == False
    ###for h_ext
    return eigen_list

  def PropertyCSM(self,code,h1,eigen_w):
    h2 = h1 + .0001
    N = 2**len(self.sites)

    """
    J = "J_{0:.3f}_{1:.3f}_{2:.3f}_{3:.3f}".format(Jzz,Jpm,Jpmpm,Jzpm)
    print code
    print h1
    print h2
    hstring1 = "_h_{0:.6f}".format(h1)
    hstring2 = "_h_{0:.6f}".format(h2)
    filename1 = J + hstring1 + "_d_" + dir_str + "_" + code+"_eval"
    filename2 = J + hstring2 + "_d_" + dir_str + "_" + code+"_eval"
    npzfile = np.load(filename1)
    eigen1_w = npzfile['arr_0']
    npzfile = np.load(filename2)
    eigen2_w = npzfile['arr_0']
    """
    #print eigen_w[0]

    min_level1 = 0
    min_level2 = 0
    for E_j in eigen_w[0]:
      if E_j.real < min_level1:
        min_level1 = E_j.real
    for E_j in eigen_w[1]:
      if E_j.real < min_level2:
        min_level2 = E_j.real
    #print min_level

    C = np.zeros( (T_POINTS), float)
    S = np.zeros( (T_POINTS), float)
    M = np.zeros( (T_POINTS), float)

    for t in range(T_POINTS):
      Z = 0
      Z_up = 0
      E = 0
      E2 = 0

      for j in range(N):
        E_j = eigen_w[0][j].real - min_level1
        if E_j / T[t] < 300:
          Boltzmann = np.exp(-E_j / T[t])
          Z += Boltzmann
          E += E_j * Boltzmann
          E2 += E_j * E_j * Boltzmann

        E_j = eigen_w[1][j].real - min_level2
        if E_j / T[t] < 300:
          Boltzmann = np.exp(-E_j / T[t])
          Z_up += Boltzmann

      ## for j
      E /= Z
      E2 /= Z
      C[t] = (E2 - E*E) / (T[t] * T[t])
      S[t] = E / T[t] + np.log(Z)
      Z *= np.exp( -min_level1 / T[t])
      Z_up *= np.exp( -min_level2 / T[t])
      M[t] = -(-T[t]*np.log(Z_up) + T[t]*np.log(Z))/(h2-h1)
      #Z *= np.exp( -min_level1 / T[t])
      #Z_up *= np.exp( -min_level2 / T[t])
      #M[t] = -(-T[t]*np.log(Z_up) + T[t]*np.log(Z))/dh
    ## for t
    return [C,S,M]

  def Weight(self,prop,code,W0,W1,W2,W3,C1,C2,C3,C4s,C4y):
    #print code
    #print sites
    #print t1
    #print t2
    #print t3
    #print t4s
    #print t4y
    WC = np.zeros( (T_POINTS), float)
    WS = np.zeros( (T_POINTS), float)
    WM = np.zeros( (T_POINTS), float)
    W = [WC,WS,WM]

    for p in range(3):
      #print "begin sub"
      if code in C1 :
        for t in range(T_POINTS):
          W[p][t] = prop[0][p][t]
          for c in range(4):
            W[p][t] -= W0[c][p][t]

      ## if 1 tetra
      elif code in C2:
        for c in range( len(C2) ):
          if code == C2[c]:
            index2 = c
        #print "loop found this index"
        #print index2
        for t in range(T_POINTS):
          W[p][t] = prop[index2][p][t]
  
          W[p][t] -= 2*W0[0][p][t]
          W[p][t] -= 2*W0[1][p][t]
          W[p][t] -= 2*W0[2][p][t]
          W[p][t] -= 2*W0[3][p][t]
          W[p][t] +=   W0[index2][p][t]

          W[p][t] -= 2*W1[0][p][t]
      ## if 2 tetra cluster
      elif code in C3:
        for c in range( len(C3) ):
          if code == C3[c]:
            index3 = c
        #print index3
        c1 = code[2]
        c2 = code[3]
        for c in range( len(C2) ):
          if "t2"+c1 == C2[c]:
            i1 = c
          if "t2"+c2 == C2[c]:
            i2 = c
        for t in range(T_POINTS):
          W[p][t] = prop[index3][p][t]
  
          W[p][t] -= 3*W0[0][p][t]
          W[p][t] -= 3*W0[1][p][t]
          W[p][t] -= 3*W0[2][p][t]
          W[p][t] -= 3*W0[3][p][t]
          W[p][t] +=   W0[i1][p][t]
          W[p][t] +=   W0[i2][p][t]

          W[p][t] -= 3*W1[0][p][t]
  
          W[p][t] -= W2[i1][p][t]
          W[p][t] -= W2[i2][p][t]
      ## if 3 tetra
      elif code in C4s:
        for c in range( len(C4s) ):
          if code == C4s[c]:
            index4 = c
        c1 = code[3]
        c2 = code[4]
        c3 = code[5]
        for c in range( len(C2) ):
          if "t2"+c1 == C2[c]:
            i1 = c
          if "t2"+c2 == C2[c]:
            i2 = c
          if "t2"+c3 == C2[c]:
            i3 = c
        """
        print c1
        print c2
        print c3
        print i1
        print i2
        print i3
        """

        #build string codes
        if c1 < c2:
          c10 = c1+c2
        else:
          c10 = c2+c1
        if c2 < c3:
          c20 = c2+c3
        else:
          c20 = c3+c2
        #print c10
        #print c20
        for c in range( len(C3) ):
          if "t3"+c10 == C3[c]:
            j10 = c
          if "t3"+c20 == C3[c]:
            j20 = c

        #print j10
        #print j20
    
        for t in range(T_POINTS):
          W[p][t] = prop[index4][p][t]

          W[p][t] -= 4*W0[0][p][t]
          W[p][t] -= 4*W0[1][p][t]
          W[p][t] -= 4*W0[2][p][t]
          W[p][t] -= 4*W0[3][p][t]
          W[p][t] +=   W0[i1][p][t]
          W[p][t] +=   W0[i2][p][t]
          W[p][t] +=   W0[i3][p][t]

          W[p][t] -= 4*W1[0][p][t]
  
          W[p][t] -= W2[i1][p][t]
          W[p][t] -= W2[i2][p][t]
          W[p][t] -= W2[i3][p][t]

          W[p][t] -= W3[j10][p][t]
          W[p][t] -= W3[j20][p][t]

     ## 4 tetrahedron string

      elif code in C4y:
        for c in range( len(C4y) ):
          if code == C4y[c]:
            index5 = c
          c1 = code[3]
          c2 = code[4]
          c3 = code[5]
          for c in range( len(C2) ):
             if "t2"+c1 == C2[c]:
               i1 = c
             if "t2"+c2 == C2[c]:
               i2 = c
             if "t2"+c3 == C2[c]:
               i3 = c
        """
        print c1
        print c2
        print c3
        print i1
        print i2
        print i3
        """
        #build string codes
        if c1 < c2:
          c10 = c1+c2
        else:
          c10 = c2+c1
        if c2 < c3:
          c20 = c2+c3
        else:
          c20 = c3+c2
        if c1 < c3:
          c30 = c1+c3
        else:
          c30 = c3+c11
        #contrast to t4s types which only have 2 three tetra subgraphs
        #print c10
        #print c20
        #print c30
        for c in range( len(C3) ):
          if "t3"+c10 == C3[c]:
            j10 = c
          if "t3"+c20 == C3[c]:
            j20 = c
          if "t3"+c30 == C3[c]:
            j30 = c

        #print j10
        #print j20
        #print j30
    
        for t in range(T_POINTS):
          W[p][t] = prop[index5][p][t]
  
          W[p][t] -= 4*W0[0][p][t]
          W[p][t] -= 4*W0[1][p][t]
          W[p][t] -= 4*W0[2][p][t]
          W[p][t] -= 4*W0[3][p][t]
          W[p][t] +=   W0[i1][p][t]
          W[p][t] +=   W0[i2][p][t]
          W[p][t] +=   W0[i3][p][t]

          W[p][t] -= 4*W1[0][p][t]

          W[p][t] -= W2[i1][p][t]
          W[p][t] -= W2[i2][p][t]
          W[p][t] -= W2[i3][p][t]

          W[p][t] -= W3[j10][p][t]
          W[p][t] -= W3[j20][p][t]
          W[p][t] -= W3[j30][p][t]

     ## 4 tetrahedron y
    return W

