
import math, random

import pickle

import numpy
import scipy.integrate
import matplotlib.pyplot as plot

import cchem

main = cchem.MainRoot()

au = {}
au["length"] = 5.291772109217e-11
au["mass"] = 9.1093829140e-31
au["time"] = 2.41888432650516e-17
au["energy"] = 4.3597441775e-18

T = 273.0
#k = (1.380648813e-23) * au["length"]**(-2) * au["mass"]**(-1) * au["time"]**(2)
k = (1.380648813e-23)
B = 1 / (k * T)

def virial2_int(r, u):
  return (1 - math.exp(-B * u))*r**2

def virial_units(B2):
  # m^3 -> cm^3 / mol
  return 6.02214129e23 * 1e6 * (2*math.pi) * B2

class DiatomPotential(object):
  def __init__(self, atom):
    self.state1 = {"type":"single_atom", "theory":"ccsd", "element":atom, "charge":0, "nopen":0, "basis":"aug-cc-pvdz"}
    self.state2 = {"type":"double_atom", "theory":"ccsd", "element":atom, "charge":0, "nopen":0, "basis":"aug-cc-pvdz", "r":None}

  def create(self):
    data = {}
    data["dissoc"] = self.compute_dissoc()
    data["points"] = []
    return data
    
  def load(self, data):
    self.dissoc_energy = data["dissoc"]
    self.points = data["points"]
  
  def save(self):
    data = {}
    data["dissoc"] = self.dissoc_energy
    self.points = sorted(self.points, key=lambda p: p[0])
    data["points"] = self.points
    return data
    
  def compute_dissoc(self):
    task_output = main.pydriver.nwchem_task( cchem.single_atom(main,self.state1) )
    if task_output["returncode"] != 0:
      print(task_output["stdout"])
    print(task_output.keys())
    atomic_energy = task_output["task_output"]["task:energy"]
    return 2 * atomic_energy
    
  def sample(self, r):
    r = r / 1e-10
    self.state2["r"] = r
    task_output = main.pydriver.nwchem_task( cchem.double_atom(main,self.state2) )
    diatom_energy = task_output["task_output"]["task:energy"]
    result = (r * 1e-10, (diatom_energy - self.dissoc_energy) * au["energy"]) 
    self.points.append( result )

class LennardJonesPotential(object):
  def __init__(self, epsilon, sigma):
    self.epsilon = epsilon
    self.sigma = sigma
    self.sz = 20
    
  def sample(self, r):
    return 4 * self.epsilon * ( (self.sigma / r)**12 - (self.sigma / r)**6 )

def make_plot(name, xs, ys):
  print("plotting {0}: {1} pts".format(name, len(xs)))
  plot.plot(xs,ys)
  plot.xlim(0.0, 20e-10)
  plot.ylim(-0.5e-20, 2e-20)
  plot.savefig("pot_{0}.jpg".format(name))
  plot.cla()
  ys2 = list(map(lambda p: virial2_int(*p), zip(xs,ys)))
  plot.plot(xs,ys2)
  plot.xlim(0.0, 20e-10)
  plot.ylim(-6e-20, 6e-20)
  plot.savefig("v2int_{0}.jpg".format(name))
  plot.cla()
  ys3 = list(map(lambda y: abs(y)**(1/2.0), numpy.diff(ys2, n=1) / numpy.diff(xs, n=1)))
  plot.plot(xs[0:-1],ys3)
  plot.ylim(min(ys3), max(ys3))
  plot.savefig("v2int_deriv1_{0}.jpg".format(name))
  plot.cla()

def compute_LJ(lj_pot):
  result, uncertain = scipy.integrate.quad( lambda x: virial2_int(x, lj_pot.sample(x)), 5e-11, 20.0e-10)
  print("lennard-jones:", virial_units(result))
  xs, ys = sample_LJ(lj_pot)
  make_plot("ar_lj", xs, ys)

def integrate(xs, ys):
  N = len(xs)-1
  trap_area = 0.0
  for i in range(0, N-1):
    trap_area += (xs[i+1] - xs[i]) * ((ys[i] + ys[i+2]) / 2)
  return trap_area
  
def compute_AB(ab_init, name):
  (xs, ys) = zip(*ab_init.points)
  result = integrate(xs, list(map(lambda p: virial2_int(*p), zip(xs,ys))))
  print("ab-initio:", virial_units(result))
  make_plot(name, xs, ys)

def sample_AB_grid(ab):
  for x in numpy.linspace(1e-10, 20e-10, 20):
    print("ab sample", x)
    try: 
      ab.sample(x)
    except KeyboardInterrupt:
      raise 
    except:
      print("did not converge")
  print("sampling finished")
  
def sample_AB(ab, pts):
  for i in range(pts):
    r = random.uniform(1e-10, 20e-10)
    print("ab sample", i, r)
    try:
      ab.sample(r)
    except KeyboardInterrupt:
      raise 
    except:
      print("did not converge")
  print("sampling finished")

def sample_LJ(lj_pot, num_pts=800):
  xs = numpy.linspace(0.001e-10, 20e-10, num_pts)
  ys = []
  for x in xs:
    ys.append( lj_pot.sample(x) )
  return (xs, ys)

main.set_db("virial")
main.start()

with open("interatomic.pckl", "rb") as f:
  interatomic = pickle.load(f)
  
lj_pot_ar = LennardJonesPotential(1.65e-21, 3.4e-10)

lj_pot_ne = LennardJonesPotential(5.084e-22, 2.782e-10)
compute_LJ(lj_pot_ne)

#ab_init_ar = DiatomPotential("Ar")
#interatomic["Ar"] = ab_init_ar.create()
#ab_init_ar.load( interatomic["Ar"] )
#interatomic["Ar"] = ab_init_ar.save()

ab_init_ne = DiatomPotential("Ne")
#interatomic["Ne"] = ab_init_ne.create()
ab_init_ne.load( interatomic["Ne"] )
sample_AB(ab_init_ne, 50)
interatomic["Ne"] = ab_init_ne.save()
print(min(list(map(lambda x: x[0], ab_init_ne.points))))
compute_AB(ab_init_ne, "ab_init")

#ab_init_ne_grid = DiatomPotential("Ne")
#interatomic["Ne_grid"] = ab_init_ne_grid.create()
#ab_init_ne_grid.load( interatomic["Ne_grid"] )
#sample_AB_grid(ab_init_ne_grid)
#interatomic["Ne_grid"] = ab_init_ne_grid.save()
#compute_AB(ab_init_ne_grid, "ab_init_grid")

with open("interatomic.pckl", "wb") as f:
  pickle.dump(interatomic, f)

main.finished()
