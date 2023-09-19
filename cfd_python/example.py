import numpy as np
import matplotlib.pyplot as plt
from time import sleep

def T_right(T,i,BC =[0,0]):
  if(i==len(T)-1):
    return BC[1]
  return T[i+1]
def T_left(T,i,BC = [0,0]):
  if(i==0):
    return BC[0]
  return T[i-1]
def T_p(T,i,BC=[0,0]):
  return T[i]

def gamma(T):
  if (T>T_crit) and (T<T_crit2):
    return 0.001
  elif(T<T_crit):
    return 0.5
  else:
    return 0.05
  
def source(T):
  if(T>T_crit) and (T<T_crit2):
    return 1e4
  elif(T<T_crit):
      return 1e2
  else:
    return 1e2

def dsource(T):
  return 0

def solve(T,dx,coeff,s,ds,BC):
  for i in range(len(T)):
    T_p = T[i]
    T_r =T_right(T,i,BC)
    T_l =T_left(T,i,BC)
    T[i] =  - (- (s(T_p)-ds(T_p))*dx**2 - coeff(T_l) * T_l -coeff(T_r) * T_r) / (coeff(T_l) + coeff(T_r) + ds(T_p) * dx ** 2)
  return T

def reverse_solve(T,dx,coeff,s,ds,BC):
  for i in range(len(T)-1,-1,-1):
    T_p = T[i]
    T_r =T_right(T,i,BC)
    T_l =T_left(T,i,BC)
    T[i] =  - (- (s(T_p)-ds(T_p))*dx**2 - coeff(T_l) * T_l -coeff(T_r) * T_r) / (coeff(T_l) + coeff(T_r) + ds(T_p) * dx ** 2)
  return T

def plot(c,T,BC):
  x = np.insert(c,0,c[0]-(c[1]-c[0]))
  x = np.append(x,c[-1]+(c[-1]-c[-2]))
  y = np.insert(T,0,BC[0])
  y = np.append(y,BC[1])
  plt.plot(x,y)
  plt.scatter(x,y)
  plt.show()

x = np.linspace(0,1,11)
bc = [500,300]
center = (x[1:] + x[:-1]) / 2
dx = x[1]-x[0]
T_init = 300
T = np.ones_like(center) * T_init
T = T + np.array([(len(T)-i) / len(T) * 10 for i in range(len(T))])

T_crit = 400
T_crit2 = 420
iter = 10
verbose = 1

plot(center,T,bc)

for i in range(iter):
  solve(T,dx,gamma,source,dsource,bc)
  if((i+1)%verbose==0):
    plot(center,T,bc)