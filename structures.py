
class Tetra:
  def __init__(self,name):
    self.site = [ 0, 1, 2, 3]
    self.bond = [ [0,1], [0,2], [0,3], [1,2], [1,3], [2,3] ]
    self.name = name
  def __str__(self):
    return str( (self.name,"tetra") )

class HexRing:
  def __init__(self,name):
    self.site = [ 0, 1, 2, 3, 4, 5 ]
    self.bond = [ [0,1], [1,2], [2,3], [3,4], [4,5], [5,0] ]
    self.name = name
  def __str__(self):
    return str( (self.name,"ring") )

class Site:
  def __init__(self,name):
    self.site = [0,]
    self.name = name
  def __str__(self):
    return str( (self.name,"site") )

