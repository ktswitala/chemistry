
import sys
import math
from functools import partial
import copy, collections

import statistics
import matplotlib.pyplot as plt

import cchem

class NormalizeHelper(object):
  def __init__(self, f):
    self.f = f
    self.f_norm = RejectionIntegrator.Tally()

  def tally(self,x,y,z,u):
    fv = self.f(x,y,z)
    self.f_norm.tally(u, fv**2)
      
  def result(self, measure):
    return 1 / math.sqrt(self.f_norm.result() * measure)
  
class ProductHelper(object):
  def __init__(self, f1, f2):
    self.f1norm = RejectionIntegrator.Tally()
    self.f2norm = RejectionIntegrator.Tally()
    self.f1f2 = RejectionIntegrator.Tally()
    self.f1 = f1
    self.f2 = f2
  
  def max_u(self, max_u):
    self.f1norm.max_u = max_u
    self.f2norm.max_u = max_u
    self.f1f2.max_u = max_u
     
  def tally(self,x,y,z,u):
    f1v = self.f1(x,y,z)
    f2v = self.f2(x,y,z)
    self.f1norm.tally(u, f1v**2)
    self.f2norm.tally(u, f2v**2)
    self.f1f2.tally(u, f1v*f2v)

  def efficiency(self):
    print("efficiency {0} {1} {2}".format(*map(lambda t: t.efficiency(), [self.f1norm, self.f2norm, self.f1f2])))
    print("out of range {0} {1} {2}".format(*map(lambda t: t.out_of_range, [self.f1norm, self.f2norm, self.f1f2])))

  def result(self, measure):
    norm1 = 1 / math.sqrt(self.f1norm.result() * measure)
    norm2 = 1 / math.sqrt(self.f2norm.result() * measure)
    return (self.f1f2.result()*measure) * norm1 * norm2

def integrate(helper):
  results = []
  for i in range(0,10):
    region = RegionSampler.MultipleIndependent(
      RegionSampler.SphereRegion([0.0, 0.0, 0.0], [0.0,2.0]),
      RegionSampler.HyperCubeRegion([[0.0, 1.0]]))
    
    uri = RejectionIntegrator.UniformIntegrator(helper, region)
    uri.integrate(int(1e5))
    result = uri.result() 
    results.append( uri.result() )
  return results
  
def normalize(g):
  results = integrate(NormalizeHelper(g))
  #print("N1: ", g.normalize())
  mean = statistics.mean(results)
  stdev = statistics.pstdev(results, mu=mean)
  return (mean, stdev)
      
def overlap(g1, g2):
  results = integrate(ProductHelper(g1,g2))
  #print("overlap: ", g1.overlap(g2))
  #print("N1: ", g1.normalize())
  #print("N2: ", g2.normalize())
  mean = statistics.mean(results)
  stdev = statistics.pstdev(results, mu=mean)
  return (mean, stdev)

def herpplot(results):
  plt.xlim((-1,1))
  plt.hist(results, bins=[-1 + 0.02*x for x in range(0,101)])
  plt.show()

#BasisSet.create_db()
basis_library = BasisSet.load_basis_sets()
basis_library.index_all_basis()
print("libraries loaded")
basis_library.report()

fns = list(basis_library.search_basis_funcs(name="sto-2g"))

overlaps = []
for name, f1 in fns:
  overlap_row = []
  for name, f2 in fns:
    print("{0} {1} - {2} {3}".format(f1["atom"], f1["angular"], f2["atom"], f2["angular"]))
    g1 = GaussianBasis.Contraction(f1["angular"], f1["coeff"])
    g2 = GaussianBasis.Contraction(f2["angular"], f2["coeff"])
    ol = overlap(g1, g2)
    overlap_row.append(ol)
    print(ol)
  overlaps.append(overlap_row)

def test_0():
  def g(x,y,z):
    pass
  print("test_0")
  n1 = 1000000
  start_time = time.time()
  for i in range(n1):
    x = random.uniform(-2.0, 2.0)
    y = random.uniform(-2.0, 2.0)
    z = random.uniform(-2.0, 2.0)
    result = g(x,y,z)
  print((time.time() - start_time) / n1)
    
def test_1():
  g = Primitive(1.5, {"x":1, "y":2, "z":0})
  print("test_1")
  n1 = 1000000
  start_time = time.time()
  for i in range(n1):
    x = random.uniform(-2.0, 2.0)
    y = random.uniform(-2.0, 2.0)
    z = random.uniform(-2.0, 2.0)
    result = g(x,y,z)
  print((time.time() - start_time) / n1)

def test_2():
  g = Contraction("s", [ [0.5, [0.5]], [0.25, [0.25]], [0.125, [0.125]] ])
  print("test_2")
  n1 = 1000000
  start_time = time.time()
  for i in range(n1):
    x = random.uniform(-2.0, 2.0)
    y = random.uniform(-2.0, 2.0)
    z = random.uniform(-2.0, 2.0)
    result = g(x,y,z)
  print((time.time() - start_time) / n1)
    
if __name__ == "__main__":
  import random, time
  test_0()
  test_1()
  test_2()
