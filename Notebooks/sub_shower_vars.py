# coding: utf-8
atm      = EGS5Simulations.LayeredAtmosphere('Parameters/atmprof6.dat')
ztop     = atm.topOfAtmosphere()
z0       = 0
nlayer   = 100
bfield   = None
nmedia   = 1
emax     = 10000000
zn       = 0.0
w0       = cos(zn/180.0*pi)
res      = 0.02
fov      = 5.0
scopes = [ ]
for i in range(0,5):
    scopes.append([i*60/sqrt(2),i*60/sqrt(2),12])
eprimary = 100000.0
eprune   = 100.0
nsim     = 200
nsubsim  = 20
mfpfrac  = 0.5
d        = 0.0
