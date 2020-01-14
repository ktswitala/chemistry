
import os, subprocess

class TaskRunner(object):
  def __init__(self, task_input, directory):
    self.task_input = task_input
    self.directory = directory
    
  def run(self):
    nw_path = os.path.join(self.directory, "task.nw")
    with open(nw_path, "wb") as f:
      f.write( self.task_input.result() )
      
    process = subprocess.Popen("nwchem task.nw", stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.directory, shell=True)
    self.stdout, self.stderr = process.communicate()
    print(self.stdout.decode('ascii'))
