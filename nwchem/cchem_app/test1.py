
import cchem

print("Loading basis sets")
#basis_sets = cchem.nwchem.load_basis("basis.pickle", recreate=True)
basis_sets = cchem.nwchem.load_basis("basis.pickle")
print("Done")

def prepare(do_geometry):
  nw = cchem.nwchem.Input()
  space = cchem.geometry.PhysicalSpace()
  do_geometry(space)
  theory = cchem.nwchem.SCFModule(nopen=space.get_spin_state(0))
  #basis_set_ecp = basis_sets["stuttgart_rlc_ecp"]
  basis_set = basis_sets["pc-3"]
  basis_set = basis_sets["crenbl_ecp"]
  
  basis = cchem.nwchem.IdenticalBasis(basis_set)

  nw.add_directive("charge", 0.0)
  nw.add_directive("geometry", space)
  nw.add_directive("scf", theory)
  nw.add_directive("basis", basis)
  if basis_set.assoc_ecp:
    nw.add_directive("ecp", basis)
  nw.add_directive("task", "scf", "energy")

  return nw
  
def single_atom(ps, element):
  ps.add_nucleus( element, cchem.geometry.Location(0.0, 0.0, 0.0) )

def double_atom(ps, element1, element2):
  ps.add_nucleus( element1, cchem.geometry.Location(0.0, 0.0, 0.0) )
  ps.add_nucleus( element2, cchem.geometry.Location(0.0, 0.0, 1.2) )

def doit():
  nw = prepare( lambda space: double_atom(space, 'O', 'O') )
  tr = cchem.nwchem.TaskRunner(nw, "./nwchem_tasks/task0/")
  tr.run()

def view():
  for basis_set in basis_sets.values():
    if 'O' in basis_set.by_atom:
      if basis_set.by_atom['O']["ecp_data"] is True:
        print(basis_set.name)

doit()
#view()
