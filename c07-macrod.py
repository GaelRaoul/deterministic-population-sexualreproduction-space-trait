import numpy as np
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from os.path import exists
from sklearn.linear_model import LinearRegression
import sys


############################### Command line parameters reading
c=float(sys.argv[1])
K=(int)(float(sys.argv[2]))
idata_macrofile=(int)(float(sys.argv[3]))
############################### Definition of the coefficients of the model/discretisation and variables
A=1/(2*1);
B=1.2;
X=15; Y=4+X*B
dx=0.1;dy=0.1;dt=0.2;
T=2
Kt=int(T/dt)
Kx=int(X/dx);Ky=int(Y/dy)
dx=X/Kx; Nx=Kx+1;dy=Y/Ky; Ny=Ky+1;
dt=T/Kt
sigmax=1
sigmay=1*0.8**2

x=np.linspace(0,X,Nx);ix=(np.linspace(0,Nx-1,Nx)).astype(int)
y=np.linspace(0,Y,Ny);iy=(np.linspace(0,Ny-1,Ny)).astype(int)
y=y-((y[1]+y[-1])/2-B*(x[1]+x[-1])/2)
x0ini=(x[1]+x[-1])/2
t=np.linspace(0,T,Kt+1)
it=(np.linspace(0,Kt,Kt+1)).astype(int)
x0=np.copy(x)
y0=np.copy(x)

# Definition of the initial population
n00=np.zeros((Nx,Ny))
In0=np.copy(x)
for i in ix:
  for j in iy:
    n00[i][j]=np.floor(np.exp(-(x[i]-x0ini)**2/10-(y[j]-B*x0ini)**2/10)*K*dy)-50
    n0=np.where(n00 > 0, n00, 0)
  In0[i]=np.sum(n0[i])
n=np.copy(n0)
ntemp=np.copy(n)
In=np.copy(x)

############################### Initiation of the recording and the data file.
X0rec=np.copy(t)
Xtrec1=np.copy(t)
Xtrec2=np.copy(t)
Xtrec3=np.copy(t)
itrec2=np.copy(t)
x0rec=np.copy(t)
Nrec=np.copy(t)
Inonzeros0=np.flatnonzero(In0)
truc=0

Zn0=0*np.copy(x)
for i in ix[1:len(ix)]:
  for j in iy[1:len(iy)]:
    Zn0[i]=max(0,min(Y,Zn0[i]+y[j]*n[i][j]/max(10**(-6),In0[i])))

idatafile=0
while exists(f'data/results{idatafile}.csv')==True:
  idatafile=idatafile+1
f = open(f"data/results{idatafile}.csv", "x")
f.write('t,X,Xtheta,I,x1,y1\n')
f2 = open(f"data/coefficients{idatafile}.csv", "x")
f2.write('X,Y,dx,dy,dt,K,T,sigmax,sigmay,A,B,c,speed\n')
f2.write(f"{X},{Y},{dx},{dy},{dt},{K},{T},{sigmax},{sigmay},{A},{B},{c},")

############################### Definition of the dispersion kernel
alpha=1
sigmakernelx=sigmax*dt/alpha
sigmakernely=sigmay*dt/alpha
xconv=np.linspace(-X,X,2*Nx-1)
gamma_x=np.exp(-np.square(xconv)/(2*sigmakernelx))
gamma_x=gamma_x/sum(gamma_x)
yconv=np.linspace(-Y,Y,2*Ny-1)
gamma_y=np.exp(-np.square(yconv)/(2*sigmakernely))
gamma_y=gamma_y/sum(gamma_y)
nconv_x=np.copy(n)

############################### Time loop
for k in it:
  # Core of the time loop
  ntemp=np.copy(n)
  ntemp0=np.copy(n)
  ############################### Calculations for the diffusion in space and trait
  for i in ix:
    nconv_x[i]=np.convolve(n[i],gamma_y,'valid')
  temp=np.transpose(nconv_x)
  temp0=np.copy(temp)
  for j in iy:
    temp[j]=np.convolve(temp[j],gamma_x,'valid')
  n_temp=alpha*np.transpose(temp)+(1-alpha)*ntemp
  ntemp=np.copy(n_temp)
  ############################### Draws of the Poisson random variables
  for i in ix[1:len(ix)-1]:
    In[i]=np.sum(ntemp[i])
    for j in iy[1:len(iy)-1]:
      n[i][j]=np.random.binomial(ntemp[i][j], np.exp(-((A*(y[j]-B*(x[i]-c*t[k]))**2+In[i]/K))*dt), size=1)+np.random.binomial(ntemp[i][j], np.exp((1-((A*(y[j]-B*(x[i]-c*t[k]))**2+In[i]/K)))*dt)-np.exp(-((A*(y[j]-B*(x[i]-c*t[k]))**2+In[i]/K))*dt), size=1)
  # Boundary conditions
  In[0]=np.sum(ntemp[0]);
  In[ix[-1]]=np.sum(ntemp[ix[-1]])
  for j in iy[1:len(iy)-1]:
    n[0][j]=0
    n[-1][j]=0  
  # Recording of the data
  Inonzeros=np.flatnonzero(In)
  Igrandtemp=np.where(In-K/100>0,1,0)
  Igrand=np.flatnonzero(Igrandtemp)
  pos_front=Igrand[-1]
  Xtrec1[k]=x[Inonzeros[1]]
  Xtrec2[k]=x[Inonzeros[-1]]
  Xtrec3[k]=x[pos_front]
  Nrec[k]=np.sum(In)
  X0rec[k]=x[int(len(ix)/2)]
  itrec2[k]=Inonzeros[-1]
  x0rec[k]=truc
  #print(np.sum(In))
  f.write(f"{t[k]}, {Xtrec2[k]}, {Xtrec3[k]},{sum(In)*dx},{B*(x[1]+x[-1])/2},{(y[1]+y[-1])/2}\n")
  # Shift of the calculation window as the front moves
  if Xtrec3[k]> (0.3*x[1]+0.7*x[-1]):
    #print(np.max(In))
    truc=truc+1
    shift=5
    x=x+shift*dx
    ntemp=np.copy(n)
    n=np.zeros((Nx,Ny))
    for i in ix[0:len(ix)-shift]:
      for j in iy[0:len(iy)]:
        n[i][j]=ntemp[i+shift][j]
  if Xtrec3[k]< (0.5*x[1]+0.5*x[-1]):
    #print(np.max(In))
    truc=truc+1
    shift=5
    x=x-shift*dx
    ntemp=np.copy(n)
    n=np.zeros((Nx,Ny))
    for i in ix[shift:len(ix)]:
      for j in iy[0:len(iy)]:
        n[i][j]=ntemp[i-shift][j]
    for i in ix[0:shift-1]:
      for j in iy[0:len(iy)-(int)(i*B*dx/dy)]:
        n[i][j]=np.max([ntemp[shift][j+(int)(i*B*dx/dy)],ntemp[shift][j]])#n[i][j]=ntemp[shift][j+(int)(i*B*dx/dy)]
        #n[i][j]=ntemp[shift][j+(int)(i*B*dx/dy)]
        #n[i][j]=ntemp[shift][j]
  deltay=(B*((x[1]+x[-1])/2-c*t[k])-(y[-1]+y[1])/2)
  #yshift=0
  if deltay>1:
    yshift=5
    y=y+yshift*dy
    ntemp=np.copy(n)
    n=np.zeros((Nx,Ny))
    for i in ix[0:len(ix)]:
      for j in iy[0:len(iy)-yshift]:
        n[i][j]=ntemp[i][j+yshift]
  if deltay<-1:
    yshift=-5
    y=y+yshift*dy
    ntemp=np.copy(n)
    n=np.zeros((Nx,Ny))
    for i in ix[0:len(ix)]:
      for j in iy[-yshift:len(iy)]:
        n[i][j]=ntemp[i][j+yshift]


############################### Calculation and recording of the propagation speed
iinf=int(len(t)/2)
tcalc=t[iinf:len(t)]#[np.floor(len(t)/2):len(t)-1]
Xcalc=Xtrec3[iinf:len(t)]#[np.floor(len(t)/2):len(t)-1]
xlin=tcalc.reshape((-1, 1))
ylin=Xcalc
model =  LinearRegression().fit(xlin, ylin)
f2.write(f"{model.coef_[0]}")


f3 = open(f"data/data_macro{idata_macrofile}.csv", 'a')
f3.write(f"{X},{Y},{dx},{dy},{dt},{K},{T},{sigmax},{sigmay},{A},{B},{c},{model.coef_[0]}\n")
print(f"{X},{Y},{dx},{dy},{dt},{K},{T},{sigmax},{sigmay},{A},{B},{c},{model.coef_[0]}")


f.close()
f2.close()
