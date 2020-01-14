
import sys, os
import socket
import pickle

from useful.ipc.RPCHandler import RPCHandler

def iterate_names():
  next_name = rtdb_first()
  while next_name:
    yield next_name
    try:
      next_name = rtdb_next()
    except:
      next_name = None
    
def dump_db():
  db = {}
  for name in iterate_names():
    db[name] = rtdb_get(name)
  return db

def debug_parse(value):
  sys.stderr.write(value)
  input_parse(value)

sys.stderr.write("{0}\n".format(os.getcwd()))
with open("params.pckl", "rb") as f:
  params = pickle.load(f)
  
action_dict = {}
action_dict["input_parse"] = debug_parse
action_dict["task"] = lambda *args, **kwargs: task(str(args[0]))
action_dict["get_db"] = dump_db

rpc_socket = socket.create_connection((params["host"], params["port"]))
rpc = RPCHandler(rpc_socket)
command_stream = rpc.new_stream("command")
info_stream = rpc.new_stream("info")

rpc.send_object(command_stream, {"type":"helo", "resource_id":params["resource_id"]})

looping = True
while looping:
  obj = rpc.read_object(command_stream)
  if obj is None:
    looping = False
    continue
  if obj["type"] == "action":
    rpc.handle_call(command_stream, action_dict, obj)
  else:
    raise Exception("unknown obj type")

