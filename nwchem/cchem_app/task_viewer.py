
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *

import pymongo

class TaskSearch(QWidget):
  def __init__(self, tasks):
    self.tasks = tasks
    self.search_box = 
    
  def 
    
db_client = pymongo.MongoClient()
db = db_client["test1"]
tasks = db["tasks"]

app = QApp()
mw = MainWindow()
