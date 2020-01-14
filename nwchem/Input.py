
import io
import os

def write_module(module):
  return module.text()

def write_bool(flag, t=None, f=None):
  if flag is True:
    return "{0}".format(t)
  if flag is False:
    return "{0}".format(f)

def write_option(label=None):
  return "{0}".format(label)
  
# pls escape this
def write_str(s):
  return '"{0}"'.format(s)
    
def write_arg(arg, arg_type):
  if arg_type is str:
    return write_str(arg)
  elif arg_type is bool:
    if arg is True:
      return "1"
    if arg is False:
      return "0"
  else:
    return arg
        
def write_args(args, arg_types=None, label=None):
  s = "{0}".format(label)
  for i, arg in enumerate(args):
    arg_type = arg_types[i]
    s += " {0}".format(write_arg(arg, arg_type))
  return s
        
def write_set(name, ty, datas):
  s = "set {0} {1}".format(name, ty)
  type_convert = {"integer":int, "real":float, "double":float, "logical":bool, "string":str}
  type_arg = type_convert[s["type"]]
  for data in datas:
    s += write_arg(data, type_arg)
  return s
      
def write_unset(names):
  s = "unset"
  for name in names:
    s += " {0}".format(name)
      
def write_startup_mode(mode, file_prefix=None, rtdb_filename=None):
  if not mode in ["start", "restart"]:
    raise Exception("invalid startup mode: {0}".format(mode))

  s = mode
  if "file_prefix":
    s += " {0}".format(file_prefix)
  if "rtdb_filename":
    s += " rtdb {0}".format(rtdb_filename)
  return s
  
def write_dirs(dir_type, process={}, host={}, other=[]):
  if not dir_type in ["permanent_dir", "scratch_dir"]:
    raise Exception("invalid dir type: {0}".format(dir_type))
  s += dir_type
  for i, d in process.items():
    s += " {0}:{1}".format(i, d)
  for n, d in host.items():
    s += " {0}:{1}".format(n, d)
  for d in other:
    s += " {0}".format(d)
  return s

def write_memory(amounts, units=None, verify=None, hardfail=None):
  s = "memory"
  for mem_ty, amt in mem["amounts"].items():
    s += " {0} {1} mb".format(mem_ty, amt)
  if units:
    s += " units {0}".format(units)
  if verify:
    s += " {0}".format(write_bool(verify, t="verify", f="noverify"))
  if hardfail:
    s += " {0}".format(write_bool(hardfail, t="hardfail", f="nohardfail"))
  return s
    
def write_title(title):
  return "title {0}".format(title)
    
def write_print(names, level=None):
  s = ""
  if level:
    s += "print {0}".format(level)
  else:
    s += "print"
  for name in pr["names"]:
    s += " {0}".format(name)
  return s

def write_noprint(names):
  s = "noprint"
  for name in names:
    s += " {0}".format(write_str(name))
  return s
  
def write_task(task, oper):
  valid_theory = ["scf", "dft", "sodft", "mp2", "direct_mp2", "rimp2", "ccsd", "ccsd(t)", "mcscf", "selci", "pspw", "band", "tce"]
  valid_oper = ["energy", "gradient", "optimize", "saddle", "hessian", "freq", "vscf", "property", "dynamics", "thermodynamics"]
  valid_special = ["python", "rtdbprint", "cphf", "property", "dplot"]
  if task in valid_theory:
    return "task {0} {1}".format(task, oper)
  elif task in valid_special:
    return "task {0}".format(task)
      
def write_geometry(geometry):
  s =  "geometry\n"
  s += "  symmetry c1\n"
  for nucleus in geometry.get_nuclei():
    label = nucleus.element.symbol
    s += "  {0} {1} {2} {3}\n".format(label, nucleus.location.x, nucleus.location.y, nucleus.location.z)
  s += "end"
  return s
    
def write_basis(basis):
  return basis.basis_block()

def write_ecp(basis):
  return basis.ecp_block()

def write_pydriver_exec(filename):
  s = ""
  s += "python\n"
  s += """execfile("{0}")\n""".format(filename)
  s += "end"

toplevel_directives = {}
tld = toplevel_directives
tld["startup_mode"] = write_startup_mode
tld["dirs"] = write_dirs
tld["memory"] = write_memory
tld["echo"] = lambda o: write_option(label="echo")
tld["title"] = lambda o: write_args([o], arg_types=[str], label="title")
tld["set"] = lambda o: write_set
tld["unset"] = lambda o: write_unset
tld["stop"] = lambda o: write_option(label="stop")
tld["task"] = write_task
    
tld["charge"] = lambda o: write_args([o], arg_types=[int], label="charge")
tld["geometry"] = write_geometry
tld["basis"] = write_basis
tld["ecp"] = write_ecp

tld["scf"] = write_module

class Input(object):
  def __init__(self, mode="b"):
    self.out = io.BytesIO()
    self.mode = mode

    self.top_dir = None
    
    # toplevel options
    self.directives = []

  def add_directive(self, name, *args):
    self.directives.append( [name, args] )
    
  def write_all(self):
    for name, args in self.directives:
      fn = toplevel_directives[name]
      self.write( "{0}\n".format(fn(*args)) )
    
  def write(self, s):
    if self.mode == "s":
      self.out.write(s)
    elif self.mode == "b":
      self.out.write(s.encode('ascii'))
    else:
      raise Exception("unknown mode")

  def result(self):
    self.out = io.BytesIO()
    self.write_all()
    return self.out.getvalue()

class SCFModule(object):
  def __init__(self, nopen):
    self.nopen = nopen
    
  def text(self):
    s = ""
    s += "scf\n"
    s += "  nopen {0}\n".format(self.nopen)
    s += "  sym off\n"
    s += "  adapt off\n"
    s += "end"
    return s
    
class DFTModule(object):
  def __init__(self):
    pass
    
class IdenticalBasis(object):
  def __init__(self, basis_set, format_version=4):
    self.basis_set = basis_set
    self.basis_name = self.basis_set.first_group_name.replace("_", " ")
    self.format_version = format_version
  
  def basis_block(self):
    s = ""
    s += "basis\n"
    s += '  * library "{0}"\n'.format(self.basis_name)
    s += "end"
    return s

  def ecp_block(self):
    s = ""
    s += "ecp\n"
    s += '  * library "{0}"\n'.format(self.basis_set.assoc_ecp)
    s += "end"
    return s
    
