
def type_convert(i):
  if i == 1000:
    return str
  elif i == 1010:
    return int
  elif i == 1011:
    return bool
  else:
    raise Exception("Unknown type {0}".format(i))
    
def discover_types():
  infos = collections.defaultdict(lambda: set())
  for name in iterate_names():
    value = rtdb_get(next_name)
    if type(value) is list:
      if len(value) > 0:
        value_type = type(value[0])
      else:
        value_type = None
    else:
      value_type = type(value)
    infos[rtdb_get_info(next_name)[0]].add( value_type )
  print(infos)

def dictify(d_in):
  d_out = {}
  for key,value in d_in.items(): 
    ks = key.split(":")
    path = ks[:-1]
    var_name = ks[-1]
    current_d = d_out
    assigned = False
    for k in path:
      if k not in current_d:
        current_d[k] = {}
      current_d = current_d[k]
    if var_name in current_d:
      if type(current_d[var_name]):
        raise Exception("leaf replacing a branch")
      else:       
        raise Exception("value is set")
    current_d[var_name] = value
  return d_out
  
def dict_branches(d, file=None, depth=0):
  for k,v in d.items():
    print("{0}{1}".format(depth*" ", k), file=file)
    if type(v) is dict:
      dict_branches(d[k], file=file, depth=depth+1)
