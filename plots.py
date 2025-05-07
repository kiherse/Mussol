import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.colors import LogNorm

FloatType = np.float64
IntType = np.int32

############################################################# Ratpenat

ratpenat = h5py.File('/storage/scratch/lv72/lv72805/D02/J3C0/5e44/J3C0RAA-00056.h5')['00001']
time = '00'

# Ratpenat Units
Rb = 0.5 #kpc
c = 2.99792458e10 # cm/s
cdenst = 0.1 # cm^-3
mp = 1.6727e-24 # g
Gamma = 5.0/3.0
kb = 1.38e-16 # erg/K

RatpenatLength = Rb*3.085e21 #cm
RatpenatVelocity = c
RatpenatMass = cdenst * mp * RatpenatLength**3
RatpenatDensity = RatpenatMass / RatpenatLength / RatpenatLength / RatpenatLength
RatpenatPressure = RatpenatMass / RatpenatLength / RatpenatLength / RatpenatLength * RatpenatVelocity**2
RatpenatInternal = RatpenatVelocity**2

# Properties
gridsize = 250*Rb #kpc
numberofcells = 1000 #cells
resolution = gridsize/numberofcells #to kpc

# Data
mp = 1.6727e-24 # g
me = 9.10938356e-28 # g
xe = np.array(ratpenat['denstye'][::10, ::10, ::10])/np.array(ratpenat['densty'][::10, ::10, ::10])
meff = xe/me + (1-xe)/mp
density = np.array(ratpenat['densty'][::10, ::10, ::10])*RatpenatDensity
densitye = np.array(ratpenat['denstye'][::10, ::10, ::10])*RatpenatDensity
pressure = np.array(ratpenat['pres'][::10, ::10, ::10])*RatpenatPressure
internal = pressure/density/(Gamma - 1.0)
temperature = pressure/density/meff/kb
velx = np.array(ratpenat['velx'][::10, ::10, ::10])*RatpenatVelocity
vely = np.array(ratpenat['vely'][::10, ::10, ::10])*RatpenatVelocity
velz = np.array(ratpenat['velz'][::10, ::10, ::10])*RatpenatVelocity
vel = np.sqrt(velx**2 + vely**2 + velz**2) # to cm/s
kinetic = 0.5 * vel**2 # cgs
timej = ratpenat['timej'][0]*RatpenatLength/RatpenatVelocity/(365*24*3600)/(1e6)
print('Ratpenat time:',timej)
time = str(int(timej))

# RATPENAT Grid
X, Y = np.shape(density)[0],np.shape(density)[1]
x, z = np.shape(density)[0],np.shape(density)[2]
# X = np.linspace(-X/2,X/2,X+1) #to kpc
# Y = np.linspace(0,Y,Y+1) #to kpc
# x = np.linspace(-x/2,x/2,x+1) #to kpc
# z = np.linspace(-z/2,z/2,z+1) #to kpc
X = np.linspace(0,X,X+1) #to kpc
Y = np.linspace(0,Y,Y+1) #to kpc
x = np.linspace(0,x,x+1) #to kpc
z = np.linspace(0,z,z+1) #to kpc
X,Y = np.meshgrid(Y,X)
x,z = np.meshgrid(z,x)

############################################################# Maps

# Density 
fig, ax = plt.subplots(1,2,figsize=np.array([12.0, 5.0]), dpi=200)
im0 = ax[0].pcolormesh(X,Y,density[:,:,50],norm=LogNorm(),cmap='jet')
cbar0 = plt.colorbar(im0, ax=ax[0])
ax[0].set_xlabel('y [cells]',size=12)
ax[0].set_ylabel('x [cells]',size=12)
ax[0].set_aspect('equal')

im1 = ax[1].pcolormesh(x,z,density[:,300,:],norm=LogNorm(),cmap='jet')
cbar1 = plt.colorbar(im1, ax=ax[1])
cbar1.set_label('density[cgs]',size=12,weight='bold')
ax[1].set_xlabel('x [cells]',size=12)
ax[1].set_ylabel('z [kpc]',size=12)
ax[1].set_aspect('equal')

plt.suptitle('D02/'+time+'Myr')
fig.savefig('dens_'+time+'Myr.png',bbox_inches='tight')

# Pressure
fig, ax = plt.subplots(1,2,figsize=np.array([12.0, 5.0]), dpi=200)
im0 = ax[0].pcolormesh(X,Y,pressure[:,:,50],norm=LogNorm(),cmap='jet')
cbar0 = plt.colorbar(im0, ax=ax[0])
ax[0].set_xlabel('y [cells]',size=12)
ax[0].set_ylabel('x [cells]',size=12)
ax[0].set_aspect('equal')

im1 = ax[1].pcolormesh(x,z,pressure[:,300,:],norm=LogNorm(),cmap='jet')
cbar1 = plt.colorbar(im1, ax=ax[1])
cbar1.set_label('pressure[cgs]',size=12,weight='bold')
ax[1].set_xlabel('z [kpc]',size=12)
ax[1].set_ylabel('x [cells]',size=12)
ax[1].set_aspect('equal')

plt.suptitle('D02/'+time+'Myr')
fig.savefig('pres_'+time+'Myr.png',bbox_inches='tight')

# Internal energy
fig, ax = plt.subplots(1,2,figsize=np.array([12.0, 5.0]), dpi=200)
im0 = ax[0].pcolormesh(X,Y,internal[:,:,50],norm=LogNorm(),cmap='jet')
cbar0 = plt.colorbar(im0, ax=ax[0])
ax[0].set_xlabel('y [cells]',size=12)
ax[0].set_ylabel('x [cells]',size=12)
ax[0].set_aspect('equal')

im1 = ax[1].pcolormesh(x,z,internal[:,300,:],norm=LogNorm(),cmap='jet')
cbar1 = plt.colorbar(im1, ax=ax[1])
cbar1.set_label('specific internal energy[cgs]',size=12,weight='bold')
ax[1].set_xlabel('z [kpc]',size=12)
ax[1].set_ylabel('x [cells]',size=12)
ax[1].set_aspect('equal')

plt.suptitle('D02/'+time+'Myr')
fig.savefig('int_'+time+'Myr.png',bbox_inches='tight')

# Velocity
fig, ax = plt.subplots(1,2,figsize=np.array([12.0, 5.0]), dpi=200)
im0 = ax[0].pcolormesh(X,Y,vely[:,:,50]/c,cmap='jet')
cbar0 = plt.colorbar(im0, ax=ax[0])
ax[0].set_xlabel('y [cells]',size=12)
ax[0].set_ylabel('x [cells]',size=12)
ax[0].set_aspect('equal')

im1 = ax[1].pcolormesh(x,z,vely[:,300,:]/c,cmap='jet')
cbar1 = plt.colorbar(im1, ax=ax[1])
cbar1.set_label('vely [c]',size=12,weight='bold')
ax[1].set_xlabel('z [kpc]',size=12)
ax[1].set_ylabel('x [cells]',size=12)
ax[1].set_aspect('equal')

plt.suptitle('D02/'+time+'Myr')
fig.savefig('vely_'+time+'Myr.png',bbox_inches='tight')