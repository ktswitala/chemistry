
import sys, os, shutil
import time, random
import pickle

import subprocess, threading
import socketserver

import queue

import useful
from useful.ipc import RPCHandler
from useful.tasks import ResourcePool
from useful.resources import filesystem

from cchem.nwchem.Input import Input

def cpu_count():
  if multiprocessing.cpu_count() <= 2:
    return 1
  else:
    return multiprocessing.cpu_count()-2

class PyDriver(socketserver.ThreadingMixIn, socketserver.BaseRequestHandler):
  def handle(self):
    self.rpc = RPCHandler(self.request)

    self.command_stream = self.rpc.wait_stream("command")
    self.info_stream = self.rpc.wait_stream("info")
    info_actions = {"print":print}
    self.info_thread = threading.Thread(target=self.rpc.action_loop, args=("info", info_actions))
    self.info_thread.start()
    
    command = self.rpc.read_object("command")
    if command["type"] == "helo":
      self.resource_id = command["resource_id"]
      if self.resource_id not in self.server.resources:
        raise Exception("Unknown task id")
      self.resource = self.server.resources[self.resource_id]
      self.task = self.resource["task"]

      if self.resource["status"] != "connecting":
        raise Exception("Attempt to connect to task not in connection status!")
      self.resource["status"] = "running"

      try:
        self.resource["result"]["task_output"] = self.task(self)
      except:
        return

      self.rpc.shutdown()
    else:
      self.rpc.shutdown()
      raise Exception("invalid connection attempt: bad header")
  
  def send(self, *args, **kwargs):
    return self.rpc.call(self.command_stream, *args, **kwargs)
    
class PyDriverServer(object):
  def __init__(self, top_dir, max_resource=None):
    self.top_dir = top_dir

    self.ipc_dir = "/media/puter/55DED9FF64382EC5/backup_1/Projects/science/cchem/ipc"

    self.server_port = random.randint(12000,13000)

    self.server = socketserver.TCPServer(("localhost",self.server_port), PyDriver)
    self.server_thread = threading.Thread(target=self.server_loop)
    self.server_thread.start()

  def server_loop(self):
    print("Serving pydriver at port {0}".format(self.server_port))
    self.server.serve_forever()

  def shutdown(self):
    self.server.shutdown()
    
  def create_input(self):
    conngen = Input()
    conngen.task = "pydriver"
    conngen.set_top_dir(task_dir)
    conngen.set_name("conngen")
    conngen.write(os.path.join(self.ipc_dir, "py_driver.py")) 
    filesystem.remove_dir_files(conngen.permanent_dir)
    filesystem.remove_dir_files(conngen.scratch_dir)
    with open(conngen.nw_filename, "wb") as f:
      f.write(conngen.out.getvalue())
      f.flush()
    return conngen
    
  def nwchem_task(self, task):
    resource = {}
    self.server.resources[resource_id] = resource

    resource["task"] = task
    result = {}
    resource["result"] = result
    
    task_dir = os.path.join(self.top_dir, "task{0}".format(resource_id))

    conngen = self.create_input()
    params = {"host":"localhost", "port":self.server_port, "resource_id":resource_id}
    with open(os.path.join(conngen.permanent_dir, "params.pckl"), "wb") as f:
      pickle.dump(params, f, protocol=2)
      f.flush()
 
    resource["status"] = "connecting"
    
    result = {}
    result["returncode"] = process.returncode
    result["stdout"] = stdout
    result["stderr"] = stderr
    if process.returncode == 0:
      result["task_output"] = resource["result"]["task_output"]
    return result

class SimpleIOAction(object):
  def __init__(self, task_input, theory="scf"):
    self.task_input = task_input
    self.theory = theory
  
  def __call__(self, pydriver):
    pydriver.send("input_parse", self.task_input)
    if self.theory == "ccsd":
      pydriver.send("task", "tce")
    else:
      pydriver.send(self.theory)
    return pydriver.send("get_db")
