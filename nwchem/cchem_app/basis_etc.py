
import random
import collections
import pickle

import cchem

def tally_results(main):
  total = 0
  format_single = collections.defaultdict(lambda: 0)
  for task in main.tasks.find():
    if task["input"]["formatter"] == 1:
      total += 1
    output = pickle.loads(task["output"])
    if output["returncode"] == 0:
      format_single[task["input"]["formatter"]] += 1
  print(format_single)
  print(total)

def assess_quality(main):
  main.set_db('test1')
  f = open("basis_quality.txt", "w")
  
  minimums = collections.defaultdict(lambda: 100)
  average_err = collections.defaultdict(lambda: Counter())

  def relevant_tasks():
    for basis_name in basis_success.keys():
      if basis_success[basis_name].average() < 0.9:
        continue
      for task in main.tasks.find( {"input.basis":basis_name} ):
        task["output"] = pickle.loads(task["output"])
        if task["output"]["returncode"] == 0:
          yield task
          
  basis_success = count_success(main)
  
  plot.hist(list(map(lambda v: v.average(), basis_success.values())), 10, facecolor='green', alpha=0.75)
  plot.xlim(0.0, 1.0)
  plot.savefig("basis_quality.jpg")
      
  print("minimums...")
  for task in relevant_tasks():
    element = task["input"]["element"]
    energy = task["output"]["task_output"]["scf:energy"]
    if energy < minimums[element]:
      minimums[element] = energy

  print("errors...")
  for task in relevant_tasks():
    basis_name = task["input"]["basis"]
    element = task["input"]["element"]
    energy = task["output"]["task_output"]["scf:energy"]
    average_err[basis_name].add_value( (minimums[element] - energy) / minimums[element] )

  for basis_name, err in sorted([(k,v.average()) for k,v in average_err.items()], key=lambda t: t[1]):
    atom_ct = len(main.basis_sets[basis_name].by_atom.keys())
    if atom_ct < 20:
      continue
    print(basis_name.ljust(32), ("%1.2e" % (err*2600)).ljust(8), atom_ct, file=f)
  f.close()
  
def compare_nmo(main):
  main.set_db('test1')
  for task in main.tasks.find():
    output = pickle.loads(task["output"])
    element = task["input"]["element"]
    basis_set = main.basis_sets[task["input"]["basis"]]
    symbol = cchem.elements[element].symbol
    predicted_nmo = basis_set.by_atom[symbol]["total_basis_funcs"]
    if output["returncode"] == 0:
      task_db = output["task_output"]
      if predicted_nmo != task_db["scf:nmo"]:
        print("prediction wrong:", basis_set.name, symbol)
        print(predicted_nmo, task_db["scf:nmo"])
      if task_db["scf:nmo"]*2 < task_db["scf:nelec"]:
        print("not enough? ran:", basis_set.name, symbol)
        print(task_db["scf:nmo"], task_db["scf:nelec"])
    else:
      print(output["stderr"].decode('ascii'))
      if predicted_nmo*2 < element:
        print("not enough? failed:", basis_set.name, symbol)
        print(basis_set.by_atom[symbol]["total_basis_funcs"], element)

def nmo_sizes(main):
  main.set_db('test1')
  f = open("nmo_sizes.txt", "w")
  for task in main.tasks.find():
    output = pickle.loads(task["output"])
    element = task["input"]["element"]
    basis_set = main.basis_sets[task["input"]["basis"]]
    symbol = cchem.elements[element].symbol
    predicted_nmo = basis_set.by_atom[symbol]["total_basis_funcs"]
    if output["returncode"] == 0:
      print("success ", basis_set.name, symbol, predicted_nmo)
    else:
      print("failure ", basis_set.name, symbol, predicted_nmo)
  f.close()      
  
def basis_errors(main):
  print("printing errors")
  main.set_db('test1')
  for basis_set in main.basis_sets.values():
    s = ""
    for task in main.tasks.find({"input.basis":basis_set.name}):
      output = pickle.loads(task["output"])
      if output["returncode"] != 0:
        s += output["stderr"].decode('ascii') + '\n'
    if s != "":
      fname = "./output/errors/{0}.txt".format(basis_set.name)
      f = open(fname, "w")
      f.write(s)
      f.close()

def determine_name_mismatch(main, basis_name, name_mismatched, name_ok):
  basis_data = list(main.basis_db.search_basis_sets(name=basis_name))[0]
  for basis_group in basis_data["basis_groups"].values():
    element, group_name = parse_group_name(basis_group["name"])
    if group_name.lower() != basis_name.replace("_", " ").lower():
      name_mismatched.add(basis_name)
    else:
      name_ok.add(basis_name)

def find_basis_mismatch(main):
  main.set_db('test1')

  basis_success = count_success(main)

  fail_mismatches = set()
  fail_ok_name = set()
  succ_mismatches = set()
  succ_ok_name = set()
  
  for basis_name in basis_success.keys():
    if basis_success[basis_name].average() < 0.1:
      determine_name_mismatch(main, basis_name, fail_mismatches, fail_ok_name)
    elif basis_success[basis_name].average() > 0.9:
      determine_name_mismatch(main, basis_name, succ_mismatches, succ_ok_name)

  print("fail mismatch", len(fail_mismatches))
  print(fail_mismatches)
  print("fail ok", len(fail_ok_name))
  print(fail_ok_name)

  print("succ mismatch", len(succ_mismatches))
  print(succ_mismatches)
  print("succ ok", len(succ_ok_name))
  print(succ_ok_name)
  
def count_success(main):
  basis_success = collections.defaultdict(lambda: Counter())
  for basis_name in main.basis_sets.keys():
    for task in main.tasks.find( {"input.basis":basis_name} ):
      output = pickle.loads(task["output"])
      if output["returncode"] == 0:
        basis_success[basis_name].add_value(1.0)
      else:
        basis_success[basis_name].add_value(0.0)
  return basis_success
