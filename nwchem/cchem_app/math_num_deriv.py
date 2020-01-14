
import random

import math
import numpy as np
import matplotlib.pyplot as plt

def gen_coeff(M, x0, xs):
  N = len(xs)-1
  coeff = np.ones([M,N,N-1])
  return coeff

def vonde_mat(xs, order):
  N = len(xs)
  if order >= N:
    raise Exception("not enough points")
  vonde_mat = np.empty([N,order+1])
  for i in range(N):
    entry = 1
    for j in range(order+1):
      vonde_mat[i,j] = entry
      entry *= xs[i]
  return vonde_mat

xs = sorted([random.uniform(-math.pi, math.pi) for i in range(200)])
ys = list(map(math.sin, xs))

deriv_ez = np.diff(ys) / np.diff(xs)

deriv_xs = []
deriv_ys = []
for i in range(0, len(xs)-2):
  x = np.reshape(xs[i:i+3], [3,1])
  y = np.reshape(ys[i:i+3], [3,1])

  vm = vonde_mat(x, 2)
  vm_inv = np.linalg.inv(vm)
  u = np.dot(vm_inv, y)
  deriv = 2*u[2]*x[1] + u[1]
  deriv_xs.append( x[1] )
  deriv_ys.append( deriv )

plt.plot(xs[1:],deriv_ez)
plt.plot(deriv_xs, deriv_ys)
plt.show()
