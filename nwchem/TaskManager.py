
import pickle
import pymongo
import sys

class TaskManager(object):
  def __init__(self, task_collection_name, task_dir, db=None, resource_pool=None):
    self.task_collection_name = task_collection_name
    self.db = db
    self.resource_pool = resource_pool

    self.config = self.db["config"]
    self.tasks = self.db[task_collection_name]

  def run_task(self, task_input):
    resource_id = self.resource_pool.request_resource(block=True)
    
    self.resource_pool.free_resource(resource_id)
    
  def info(self):
    s = ""
    s += "ready tasks: {0}".format( self.tasks.count( {"state":"ready"} ) )
    s += "failed tasks: {0}".format( self.tasks.count( {"state":"failed"} ) )
    s += "complete tasks: {0}".format( self.tasks.count( {"state":"complete"} ) )
    
  def clear(self):
    self.db.drop_collection(task_collection_name)
    
  def create_state(self, task_input):
    entry = {}
    entry["input"] = task_input
    entry["state"] = "ready"
    self.tasks.insert(entry)
   
  def redo_tasks(self):
    for task in self.tasks.find( {"state":"complete"} ):
      output = pickle.loads(task["output"])
      if output["returncode"] != 0:
        task["state"] = "ready"
        del task["output"]
        self.tasks.update( {"_id":task["_id"]}, task )
        
  def queue_tasks(self, amount=1000):
    task_entries = self.tasks.find( {"state":"ready"} )
    for i in range(amount):
      try:
        entry = task_entries.next()
      except StopIteration:
        break
      self.runner.queue_task( self.do_task, entry )
    
  def do_task(self, entry):
    print("task ready:", entry["input"])
    input_generator = self.input_generators[entry["input"]["type"]]
    task_result = self.pydriver.nwchem_task( input_generator(self.main, entry["input"]) )
    #print( { k:len(pickle.dumps(v)) for k, v in task_result.items() } )

    if "task_output" in task_result:
      del task_result["task_output"]["basis:ao basis:contraction info"]
      del task_result["task_output"]["geometry:geometry:operators"]
      #print( { k:len(pickle.dumps(v)) for k,v in task_result["task_output"].items() } )

    entry["output"] = pickle.dumps(task_result)
    entry["state"] = "complete"
    #print(task_result["stdout"].decode('ascii'))
    #print(task_result["stderr"].decode('ascii'))
    self.tasks.update( {"_id":entry["_id"]}, entry )
    print("task finished ({0}): {1}".format(task_result["returncode"], entry["input"]))
