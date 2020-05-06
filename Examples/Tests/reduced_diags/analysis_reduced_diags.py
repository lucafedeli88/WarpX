#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the reduced diagnostics.
# The setup is a uniform plasma with electrons, protons and photons.
# Particle energy, field energy and maximum field values will be
# outputed using the reduced diagnostics.
# And they will be compared with the data in the plotfiles.

# Tolerance: 1.0e-8 for particle energy, 1.0e-3 for field energy,
# 1.0e-9 for the maximum electric field and 1.0e-18 for the
# maximum magnetic field.
# The difference of the field energy is relatively large,
# because fields data in plotfiles are cell-centered,
# but fields data in reduced diagnostics are staggered.

# Possible running time: ~ 2 s

import sys
import yt
import numpy as np
import scipy.constants as scc

fn = sys.argv[1]

ds = yt.load(fn)
ad = ds.all_data()

# PART1: get results from plotfiles

EPyt = 0.0
# electron
px = ad['electrons','particle_momentum_x'].to_ndarray()
py = ad['electrons','particle_momentum_y'].to_ndarray()
pz = ad['electrons','particle_momentum_z'].to_ndarray()
w  = ad['electrons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_e**2*scc.c**4)-scc.m_e*scc.c**2)*w )
# proton
px = ad['protons','particle_momentum_x'].to_ndarray()
py = ad['protons','particle_momentum_y'].to_ndarray()
pz = ad['protons','particle_momentum_z'].to_ndarray()
w  = ad['protons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_p**2*scc.c**4)-scc.m_p*scc.c**2)*w )

# photon
px = ad['photons','particle_momentum_x'].to_ndarray()
py = ad['photons','particle_momentum_y'].to_ndarray()
pz = ad['photons','particle_momentum_z'].to_ndarray()
w  = ad['photons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt(px**2+py**2+pz**2)*scc.c)*w )

# Use raw field plotfiles to compare with maximum field reduced diag
Ex_raw = ad['raw','Ex_aux'].to_ndarray()[:,0]
Ey_raw = ad['raw','Ey_aux'].to_ndarray()[:,0]
Ez_raw = ad['raw','Ez_aux'].to_ndarray()[:,0]
Bx_raw = ad['raw','Bx_aux'].to_ndarray()[:,0]
By_raw = ad['raw','By_aux'].to_ndarray()[:,0]
Bz_raw = ad['raw','Bz_aux'].to_ndarray()[:,0]
max_Ex = np.amax(np.abs(Ex_raw))
max_Ey = np.amax(np.abs(Ey_raw))
max_Ez = np.amax(np.abs(Ez_raw))
max_E = np.sqrt(np.amax(Ex_raw**2+Ey_raw**2+Ez_raw**2))
max_Bx = np.amax(np.abs(Bx_raw))
max_By = np.amax(np.abs(By_raw))
max_Bz = np.amax(np.abs(Bz_raw))
max_B = np.sqrt(np.amax(Bx_raw**2+By_raw**2+Bz_raw**2))

ad = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex = ad['Ex'].to_ndarray()
Ey = ad['Ey'].to_ndarray()
Ez = ad['Ez'].to_ndarray()
Bx = ad['Bx'].to_ndarray()
By = ad['By'].to_ndarray()
Bz = ad['Bz'].to_ndarray()
Es = np.sum(Ex**2)+np.sum(Ey**2)+np.sum(Ez**2)
Bs = np.sum(Bx**2)+np.sum(By**2)+np.sum(Bz**2)
N  = np.array( ds.domain_width / ds.domain_dimensions )
dV = N[0]*N[1]*N[2]
EFyt = 0.5*Es*scc.epsilon_0*dV + 0.5*Bs/scc.mu_0*dV

# PART2: get results from reduced diagnostics

EFdata = np.genfromtxt("./diags/reducedfiles/EF.txt")
EPdata = np.genfromtxt("./diags/reducedfiles/EP.txt")
MFdata = np.genfromtxt("./diags/reducedfiles/MF.txt")
EF = EFdata[1][2]
EP = EPdata[1][2]
max_Exdata = MFdata[1][2]
max_Eydata = MFdata[1][3]
max_Ezdata = MFdata[1][4]
max_Edata  = MFdata[1][5]
max_Bxdata = MFdata[1][6]
max_Bydata = MFdata[1][7]
max_Bzdata = MFdata[1][8]
max_Bdata  = MFdata[1][9]

# PART3: print and assert

max_diffEmax = max(abs(max_Exdata-max_Ex),abs(max_Eydata-max_Ey),
                   abs(max_Ezdata-max_Ez),abs(max_Edata-max_E))
max_diffBmax = max(abs(max_Bxdata-max_Bx),abs(max_Bydata-max_By),
                   abs(max_Bzdata-max_Bz),abs(max_Bdata-max_B))

print('difference of field energy:', abs(EFyt-EF))
print('tolerance of field energy:', 1.0e-3)
print('difference of particle energy:', abs(EPyt-EP))
print('tolerance of particle energy:', 1.0e-8)
print('maximum difference of maximum electric field:', max_diffEmax)
print('tolerance of maximum electric field difference:', 1.0e-9)
print('maximum difference of maximum magnetic field:', max_diffBmax)
print('tolerance of maximum electric field difference:', 1.0e-18)

assert(abs(EFyt-EF) < 1.0e-3)
assert(abs(EPyt-EP) < 1.0e-8)
assert(max_diffEmax < 1.0e-9)
assert(max_diffBmax < 1.0e-18)