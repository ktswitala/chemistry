
import cchem
from cchem.nwchem import BasisSet

angular_spherical_ct = {"s":1, "p":3, "sp":4, "d":5, "f":7, "g":9, "h":11, "i":13, "k":15, "l":17, "m":19}
angular_cartesian_ct = {"s":1, "p":3, "sp":4, "d":6, "f":10, "g":15, "h":21, "i":28, "k":36, "l":45, "m":55}
angular_symbols = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm']

class BasisSetParser(object):
  def __init__(self):
    self.valid_angular = ['s', 'p', 'sp', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm']
    
  def parse_fortran_float(self, value):
    return float(value.replace('D','E'))
    
  def move_line(self, line_amt):
    self.set_line( self.line_no + line_amt )
    while not self.finished and (len(self.line) == 0 or self.line[0] == "#"):
      self.move_line(1)

  def set_line(self, line_no):
    self.line_no = line_no
    if line_no >= len(self.lines):
      self.line = None
      self.tokens = None
      self.finished = True
      return
    self.line = self.lines[self.line_no].strip(" ")
    self.tokens = self.tokenize(self.line)
    
  def tokenize(self, line):
    pos = 0
    tokens = []
    while pos < len(line):
      start_tok = pos
      if line[pos] == '"':
        pos += 1
        while pos < len(line) and line[pos] != '"':
          pos += 1
        pos += 1
      else:
        while pos < len(line) and line[pos] not in [' ', '\t']:
          pos += 1
      end_tok = pos
      tokens.append(line[start_tok:end_tok])
      while pos < len(line) and line[pos] in [' ', '\t']:
        pos += 1
    return tokens
    
  def read_basis_contraction(self):
    block = {}
    block["atom"] = self.tokens[0]
    block["angular"] = self.tokens[1].lower()
    if block["angular"] not in self.valid_angular:
      raise Exception("bad angular function {0}".format(block["angular"]))
    block["prims"] = []
    looping = True
    while looping:
      contraction = {}
      self.move_line(1)
      try:
        contraction["exp"] = self.parse_fortran_float(self.tokens[0])
        contraction["scales"] = list(map(lambda fs: self.parse_fortran_float(fs), self.tokens[1:]))
      except:
        break
      block["prims"].append(contraction)
    if len(block["prims"]) == 0:
      raise Exception("no prims")
    return block
    
  def read_basis_group(self):
    basis = {}
    basis["name"] = self.tokens[1]
    basis["params"] = list(map(str.lower, set(self.tokens[2:])))
    basis["contractions"] = []
    self.move_line(1)
    while not self.finished and self.tokens[0] != "end":
      basis["contractions"].append(self.read_basis_contraction())
    self.move_line(1)
    return basis
  
  def read_ecp_entry(self):
    entry = {}
    entry["atom"] = self.tokens[0]
    if self.tokens[1] == "nelec":
      entry["nelec"] = int(self.tokens[2])
      self.move_line(1)
    elif self.tokens[1].lower() in self.valid_angular or self.tokens[1].lower() == "ul":
      entry["angular"] = self.tokens[1].lower()
      entry["prims"] = []
      looping = True
      params = {}
      while looping:
        self.move_line(1)
        try:
          params["n"] = int(self.tokens[0])
          params["exp"] = self.parse_fortran_float(self.tokens[1])
          params["scales"] = list(map(lambda fs: self.parse_fortran_float(fs), self.tokens[2:]))
        except:
          break
        entry["prims"].append(params)
      if len(entry["prims"]) == 0:
        raise Exception("no coeffs")
    else:
      raise Exception("unknown ecp block")
    return entry
        
  def read_ecp_group(self):
    ecp = {}
    ecp["name"] = self.tokens[1]
    ecp["blocks"] = []
    self.move_line(1)
    while not self.finished and self.tokens[0] != "end":
      block = self.read_ecp_entry()
      if "nelec" in block:
        if "nelec" in ecp:
          raise Exception("multiple nelec blocks")
        else:
          ecp["nelec"] = block
      else:
        ecp["blocks"].append(block)
    self.move_line(1)
    return ecp
  
  def read_assoc_ecp(self):
    assoc_ecp = {}
    assoc_ecp["name"] = self.tokens[1]
    self.move_line(1)
    return assoc_ecp
    
  def parse(self, data):
    self.data = data
    self.lines = data.split('\n')
    self.line_no = 0
    self.finished = False

    result = {}
    result["basis_groups"] = []
    result["ecp_groups"] = []
    result["assoc_ecps"] = []
    self.move_line(0)
    while not self.finished:
      if self.tokens[0] == "basis":
        basis = self.read_basis_group()
        result["basis_groups"].append(basis)
      elif self.tokens[0] == "ecp":
        ecp = self.read_ecp_group()
        result["ecp_groups"].append(ecp)
      elif self.tokens[0] == "ASSOCIATED_ECP":
        assoc_ecp = self.read_assoc_ecp()
        result["assoc_ecps"].append(assoc_ecp)
      else:
        raise Exception("unknown top level: {0}".format(self.tokens))
    basis_set = self.process(result)
    return basis_set

  def parse_group_name(self, full):
    temp_split = full.strip('"').split('_')
    return temp_split[0], str.join("",temp_split[1:])

  def select_one(self, params, options, default):
    found = set()
    for option in options:
      if option in params:
        found.add(option)
    if len(found) == 0:
      return default
    elif len(found) == 1:
      return found.pop()
    else:
      raise Exception("multiple selections")
      
  def process(self, data):
    basis_set = BasisSet()
    
    valid_basis_group_params = ["cartesian", "spherical", "segment", "nosegment", "print", "noprint", "rel"]
    
    if len(data["ecp_groups"]) > 0:
      basis_set.has_ecp = True
    else:
      basis_set.has_ecp = False
      
    if len(data["assoc_ecps"]) > 1:
      raise Exception("too many assoc_ecps")
    if len(data["assoc_ecps"]) > 0:
      basis_set.assoc_ecp = data["assoc_ecps"][0]["name"].strip('"')
    else:
      basis_set.assoc_ecp = None
      
    for basis_group_data in data["basis_groups"]:
      basis_group = {}
      element, basis_group["name"] = self.parse_group_name(basis_group_data["name"])
      if not hasattr(basis_set, 'first_group_name'):
        basis_set.first_group_name = basis_group["name"]
        
      for param in basis_group_data["params"]:
        if param not in valid_basis_group_params:
          raise Exception("invalid param")
        
      basis_group["harmonics"] = self.select_one(basis_group_data["params"], ["cartesian", "spherical"], "cartesian")
      basis_group["segment"] = self.select_one(basis_group_data["params"], ["segment", "nosegment"], "segment")
      basis_group["print"] = self.select_one(basis_group_data["params"], ["print", "noprint"], "print")
      if "rel" in basis_group_data["params"]:
        basis_group["rel"] = True
      else:
        basis_group["rel"] = False
        
      # process basis contractions
      for contraction_data in basis_group_data["contractions"]:
        atom = contraction_data["atom"]
        angular = contraction_data["angular"].lower()
        atom_data = basis_set.by_atom[atom]
        
    for ecp_group_data in data["ecp_groups"]:
      atom = ecp_group_data["nelec"]["atom"]
      basis_set.by_atom[atom]["ecp_data"] = True
      
    return basis_set

import os
import pickle
import useful.resources.filesystem as filesystem
def parse_all_basis(lib_path="/usr/share/nwchem/libraries"):
  parser = BasisSetParser()
  for pathinfo in filesystem.list_dir_files(lib_path):
    full_path = pathinfo["path"]
    with open(pathinfo["path"], "r") as f:
      basis_name = os.path.basename(full_path)
      print(basis_name)
      basis_set = parser.parse(f.read())
      basis_set.name = basis_name
      yield basis_set

def load_basis(path, recreate=False):
  if os.path.exists(path) and recreate is False:
    with open(path, "rb") as f:
      return pickle.load(f)
  else:
    basis_sets = {}
    for basis_set in parse_all_basis():
      basis_sets[basis_set.name] = basis_set
    with open(path, "wb") as f:
      pickle.dump(basis_sets, f)
    return basis_sets
