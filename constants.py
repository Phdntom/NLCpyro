import numpy as np

#ORDER = 3 or 4
#3 is default, very fast and accurate
#4 is slow but is required for low field and temperature
ORDER = 3

#External FIELD
field_mag = 0
# field_mag is in units of Tesla
field_dir = (1,1,1)
#field_dir will be normalized later

#EXCHANGES----------------------

#Balents
"""
Jzz =  0.17
Jpm =  0.05
Jpmpm = 0.05
Jzpm = -0.14

gz = 1.8
gxy = 2.4*gz

"""

#Chang
Jzz = 0.0293
Jpm = -0.0234
Jpmpm = 0.0058
Jzpm = -0.0293

gz = 1.77
gxy = 4.18

# J's are in meV---------------

#These should stay false unless you want to add code that uses
#wave functions then you can print out w.f's and e-values to python
#compressed format...
wf_switch = False
print_eigen_switch = False

#Temperature is in Kelvins and is converted inside the program to meV
T_POINTS = 100
T = np.zeros((T_POINTS),float)
T[0] = .1;
for t in range(1,50):
  T[t] = T[t-1] + .1
for t in range(50,75):
  T[t] = T[t-1] + .2
for t in range(75,85):
  T[t] = T[t-1] + .5
for t in range(85,T_POINTS):
  T[t] = T[t-1] + 1

#-------------------------------------------------
#--------do not change below this line------------
#-------------------------------------------------

h_meV = 17.2762
T_meV = 11.6045
# 17.2762 T/meV

print "Jzz   = {0:}".format(Jzz)
print "Jpm   = {0:}".format(Jpm)
print "Jpmpm = {0:}".format(Jpmpm)
print "Jzpm  = {0:}".format(Jzpm)
print "gz    = {0:}".format(gz)
print "gxy   = {0:}".format(gxy)

print "H mag = {0:}".format(field_mag)
if field_mag > .00001:
  print "H dir = {0:}".format(field_dir)
print "T grid = {0:}".format(T)

dir_str = str(field_dir[0])+str(field_dir[1])+str(field_dir[2])

norm = 0
for i in range(3):
  norm += field_dir[i]*field_dir[i]
field_dir = np.array(field_dir)/np.sqrt(norm)
#print field_dir

for t in range(T_POINTS):
#  print T[t]
  T[t] = T[t]/T_meV
# 11.6045 K/meV, calculations are done in meV

z0 = np.array((+1,+1,+1))/np.sqrt(3)
z1 = np.array((+1,-1,-1))/np.sqrt(3)
z2 = np.array((-1,+1,-1))/np.sqrt(3)
z3 = np.array((-1,-1,+1))/np.sqrt(3)

x0 = np.array((-2,+1,+1))/np.sqrt(6)
x1 = np.array((-2,-1,-1))/np.sqrt(6)
x2 = np.array((+2,+1,-1))/np.sqrt(6)
x3 = np.array((+2,-1,+1))/np.sqrt(6)

y0 = np.cross(z0,x0)
y1 = np.cross(z1,x1)
y2 = np.cross(z2,x2)
y3 = np.cross(z3,x3)

R = [None for i in range(4)]
R[0] = np.vstack((x0,y0,z0))
R[1] = np.vstack((x1,y1,z1))
R[2] = np.vstack((x2,y2,z2))
R[3] = np.vstack((x3,y3,z3))

a = np.exp(+np.pi/3*1j)
b = np.exp(-np.pi/3*1j)
zeta = np.array(
               [ [  0, -1,  a,  b],
                 [ -1,  0,  b,  a],
                 [  a,  b,  0, -1],
                 [  b,  a, -1,  0] ]
                )
#print zeta
gamma = -zeta.conjugate()
#print gamma




