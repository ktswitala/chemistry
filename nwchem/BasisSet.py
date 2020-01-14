
import os
import pickle
import re
from collections import defaultdict

import cchem

class BasisSet(object):
  def __init__(self):
    self.name = None
    self.has_ecp = None
    self.by_atom = defaultdict(self.new_atom)
    self.by_group = {}
    
  def new_atom(self):
    atom_data = {}
    atom_data["basis_funcs"] = {}
    atom_data["ecp_data"] = False
    atom_data["total_basis_funcs"] = 0
    return atom_data
    
  def verify(self):
    for atom, atom_data in self.by_atom.items():
      nproton = cchem.elements[atom].nproton 
      if atom_data["total_basis_funcs"] * 2 < nproton:
        raise Exception("not enough funcs for {0}, have {1}, need {2}".format(atom, atom_data["total_basis_funcs"] * 2, nproton))

  def generate_stats(self):
    pass

