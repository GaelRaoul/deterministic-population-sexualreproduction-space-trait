import os 
import numpy as np
from os.path import exists
import sys

print('X,Y,dx,dy,dt,K,T,sigmax,sigmay,A,B,c,speed')

# Vector of the different climate change speeds that will be used to launch simulations.
nb=2
k_macro=(np.linspace(0,nb-1,nb)).astype(int)
c_macro=np.linspace(-0.5, -0.6, num=nb)
K_macro=100000#K_macro=10000#

# opening a file to write the results of the simulations in (without overwriting previous files)
idata_macrofile=0
while exists(f'data/data_macro{idata_macrofile}.csv')==True:
  idata_macrofile=idata_macrofile+1
f = open(f"data/data_macro{idata_macrofile}.csv", "x")
f.write('X,Y,dx,dy,dt,K,T,sigmax,sigmay,A,B,c,speed\n')
f.close()

# Launch of the simulations
for k in k_macro:
  command="python3 c07-macrod.py "+str(c_macro[k])+" "+str(K_macro)+" "+str(idata_macrofile)
  print(command)
  os.system(command)




