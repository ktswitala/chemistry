
def print_db_structure(main):
  main.set_db('test1')
  for task in main.tasks.find():
    output = pickle.loads(task["output"])
    if output["returncode"] == 0:
      with open("scf_db_keys.txt", "w") as f:
        dict_branches(dictify(output["task_output"]), file=f)
        break
