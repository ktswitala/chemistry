
def do_fit(main, db_name):
  main.set_db(db_name)
  terms = [[]]
  for i in range(1,5):
    terms.append( [("nmo", i)] )
    terms.append( [("nelec", i)] )
  reg = MultipleLinearRegression(terms)

  nmos = []
  nelecs = []
  times = []
  for task in main.tasks.find():
    output = pickle.loads(task["output"])
    if output["returncode"] == 0:
      result = output["task_output"]
      nmos.append( result["scf:nmo"] )
      nelecs.append( result["scf:nelec"] )
      times.append( result["task:walltime"] )
      reg.add({"nmo":result["scf:nmo"], "nelec":result["scf:nelec"]}, result["task:walltime"])

  reg.refit()
  reg.report()
  
  plt.scatter(nmos, times)
  plt.xlim(min(nmos), max(nmos))
  plt.ylim(min(times), max(times))
  plt.show()
  
  print("predictions...")
  predict_error = Counter()
  for task in main.tasks.find():
    output = pickle.loads(task["output"])
    if output["returncode"] == 0:
      result = output["task_output"]
      predict_in = {"nmo":result["scf:nmo"], "nelec":result["scf:nelec"]}
      prediction = reg.predict(predict_in)
      predict_error.add_value( abs(prediction - result["task:walltime"]) )
  print("average prediction error: {0}".format(predict_error.average()))
