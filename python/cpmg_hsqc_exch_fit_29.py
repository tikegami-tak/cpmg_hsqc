# Simulates CPMG data, which include off-resonance effects and pulse-width miscalibration

'''
Modify some parameters listed from line 27 to line 60.

Type as follows for each execution.

reset -f
run cpmg_hsqc_exch_fit_29.py

Alternatively, from a Linux shell, not from Jupyterlab
date; python3 cpmg_hsqc_exch_fit_29.py >& log.txt; date
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import sys
from mpl_toolkits.mplot3d import Axes3D
import os
import shutil

pi = np.pi

class Ham:

	r2i = 0		# I spin transverse relaxation rate (/sec)
	r1i = 0		# I spin longitudinal relaxation rate (/sec)
	r2s = 3		# S spin transverse relaxation rate (/sec)
	r1s = 0		# S spin longitudinal relaxation rate (/sec)

	kab = 200		# exchange rate from A to B
	kba = 200		# exchange rate from B to A

	jis = 140		# J(IS) coupling constant

	dwib = 0		# I's chemical shift difference of state B from that of state A (Hz)
	dwsb = 50		# S's chemical shift difference of state B from that of state A (Hz)

	phi = 0			# I spin pulse phase (1, 2, 3, 4)
	v1i = 0			# I spin pulse field strength (Hz)
	phs = 0			# S spin pulse phase (1, 2, 3, 4)
	v1s = 0			# S spin pulse field strength (Hz)
	offi = 0		# I spin offset (Hz)
	offs = 0		# S spin offset (Hz)

	tex = 40*10**-3	# Total relaxation delay = 4 * each relax delay
	n90 = 20*10**-6	# 90 degree pulse width for spin S
	h90 = 10*10**-6	# 90 degree pulse width for spin I
	pmis = 0.9			# mis-calibration of the pulse width for spin S
	hpmis = 1.0		# mis-calibration of the pulse width for spin I
	dmis = 1.0			# mis-setting of 1/(2J) delay

	ofsrng = 4000		# The offset range of +- ofsrng, required for bloch[1]
	gmax = 40			# How many offset points lie within +- ofsrng, required for bloch[1]
	lctm = 18			# How many Vcpmg points, required for bloch[3]
	minhcw = 15000		# Min power of 1H for CW decoupling

	sn = 4				# Tcp is divided into 4 parts. Input 4 for all the experiments (hpi = 5, 6, 7, 8) considered here.
	hpi = 8
	# 4, simple 2*cpmg(2SyIz) - Pelement - 2*cpmg(Sx), the conventional method (Anti-TROSY)
	# 5, ST-CW-CPMG {y y x -x}
	# 6, in-phase 2*cpmg(-Sy) - 180x - 2*cpmg(Sy) with 1H CW decoupling
	# 7, simple 2*cpmg(2SyIz) - Pelement - 2*cpmg(Sx), the conventional method
	# 8, AFTAC including the 60-300-60 pulses

	def __init__(self, **kwargs):

		for kw in kwargs:
			setattr(self, kw, kwargs[kw])

		# Hamiltonian for state A

		wxi = 2*pi*self.v1i*np.cos(self.phi)
		wyi = 2*pi*self.v1i*np.sin(self.phi)
		wzi = 2*pi*self.offi
		tmp1 = np.array([[1, 0], [0, 0]])
		tmpi = np.array([[-self.r2i, -wzi, wyi], [wzi, -self.r2i, -wxi], [-wyi, wxi, -self.r1i]])

		wxs = 2*pi*self.v1s*np.cos(self.phs)
		wys = 2*pi*self.v1s*np.sin(self.phs)
		wzs = 2*pi*self.offs
		tmp2 = np.array([[0, 0], [0, 1]])
		tmps = np.array([[-self.r2s, -wzs, wys], [wzs, -self.r2s, -wxs], [-wys, wxs, -self.r1s]])

		tmpis = np.kron(tmp1, tmpi) + np.kron(tmp2, tmps)

		tmp3 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		tmpisz = np.kron(tmp3, tmpi) + np.kron(tmps,tmp3)

		tmpi = np.array([[0,0,0,0,0,0,0,-pi*self.jis,0], [0,0,0,0,0,0,pi*self.jis,0,0], [0,0,0,0,0,0,0,0,0], \
		[0,0,0,0,0,-pi*self.jis,0,0,0], [0,0,pi*self.jis,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]])

		tmps = np.array([[0,0,0,0,0,0,0,pi*self.jis,0], [0,0,0,0,0,0,-pi*self.jis,0,0], [0,0,0,0,0,0,0,0,0], \
		[0,0,0,0,0,pi*self.jis,0,0,0], [0,0,-pi*self.jis,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]])

		tmp1 = np.hstack([tmpis, tmpi])
		tmp2 = np.hstack([tmps.T, tmpisz])
		hamt = np.vstack([tmp1, tmp2])

		# Hamiltonian for state B

		wxi = 2*pi*self.v1i*np.cos(self.phi)
		wyi = 2*pi*self.v1i*np.sin(self.phi)
		bwzi = 2*pi*(self.offi+self.dwib)
		tmp1 = np.array([[1, 0], [0, 0]])
		tmpi = np.array([[-self.r2i, -bwzi, wyi], [bwzi, -self.r2i, -wxi], [-wyi, wxi, -self.r1i]])

		wxs = 2*pi*self.v1s*np.cos(self.phs)
		wys = 2*pi*self.v1s*np.sin(self.phs)
		bwzs = 2*pi*(self.offs+self.dwsb)
		tmp2 = np.array([[0, 0], [0, 1]])
		tmps = np.array([[-self.r2s, -bwzs, wys], [bwzs, -self.r2s, -wxs], [-wys, wxs, -self.r1s]])

		tmpis = np.kron(tmp1, tmpi) + np.kron(tmp2, tmps)

		tmp3 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		tmpisz = np.kron(tmp3, tmpi) + np.kron(tmps,tmp3)

		tmpi = np.array([[0,0,0,0,0,0,0,-pi*self.jis,0], [0,0,0,0,0,0,pi*self.jis,0,0], [0,0,0,0,0,0,0,0,0], \
			[0,0,0,0,0,-pi*self.jis,0,0,0], [0,0,pi*self.jis,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]])

		tmps = np.array([[0,0,0,0,0,0,0,pi*self.jis,0], [0,0,0,0,0,0,-pi*self.jis,0,0], [0,0,0,0,0,0,0,0,0], \
			[0,0,0,0,0,pi*self.jis,0,0,0], [0,0,-pi*self.jis,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]])

		tmp1 = np.hstack([tmpis, tmpi])
		tmp2 = np.hstack([tmps.T, tmpisz])
		hamtb = np.vstack([tmp1, tmp2])

		# Hamiltonian for state A again

		hama = np.pad(hamt, (1, 0), 'constant')		# to become a (16, 16) matrix

		# auto- and cross-correlated cross-relaxation, ignored
		hama[3, 6] = 0
		hama[6, 3] = 0
		hama[7, 11] = 0
		hama[11, 7] = 0
		hama[8, 10] = 0
		hama[10, 8] = 0

		# approximated relaxation
		hama[7, 7] = -self.r2s			# MQ
		hama[8, 8] = -self.r2s			# MQ
		hama[9, 9] = -self.r2s			# anti-phase, SxyIz
		hama[10, 10] = -self.r2s		# MQ
		hama[11, 11] = -self.r2s		# MQ
		hama[12, 12] = -self.r2s		# anti-phase, SxyIz
		hama[13, 13] = -self.r2i		# anti-phase, SzIxy
		hama[14, 14] = -self.r2i		# anti-phase, SzIxy
		hama[15, 15] = 0				# two-spin order

		# Hamiltonian for state B

		# auto- and cross-correlated cross-relaxation, ignored
		hamtb[2, 5] = 0
		hamtb[5, 2] = 0
		hamtb[6, 10] = 0
		hamtb[10, 6] = 0
		hamtb[7, 9] = 0
		hamtb[9, 7] = 0

		# approximated relaxation
		hamtb[6, 6] = -self.r2s		# MQ
		hamtb[7, 7] = -self.r2s		# MQ
		hamtb[8, 8] = -self.r2s		# anti-phase, SxyIz
		hamtb[9, 9] = -self.r2s		# MQ
		hamtb[10, 10] = -self.r2s	# MQ
		hamtb[11, 11] = -self.r2s	# anti-phase, SxyIz
		hamtb[12, 12] = -self.r2i	# anti-phase, SzIxy
		hamtb[13, 13] = -self.r2i	# anti-phase, SzIxy
		hamtb[14, 14] = 0			# two-spin order

		tmp1 = np.hstack([hama, np.zeros((16, 15))])
		tmp2 = np.hstack([np.zeros((15, 16)), hamtb])
		self.ham = np.vstack([tmp1, tmp2])

		# exchange involving auto- and cross-correlated cross-relaxation, ignored
		self.ham[3, 0] = 0		# 2 * self.r1i * self.kba/(self.kab + self.kba)
		self.ham[6, 0] = 0		# 2 * self.r1s * self.kba/(self.kab + self.kba)
		self.ham[15, 0] = 0
		self.ham[18, 0] = 0		# 2 * self.r1i * self.kab/(self.kab + self.kba)
		self.ham[21, 0] = 0		# 2 * self.r1s * self.kab/(self.kab + self.kba)
		self.ham[30, 0] = 0

		for i in range(1, 16):
			self.ham[i, i] -= self.kab

		for i in range(16, 31):
			self.ham[i, i] -= self.kba

		for i in range(1, 16):
			j = i + 15
			self.ham[i, j] += self.kba
			self.ham[j, i] += self.kab

class Vec:

	def __init__(self, **kwargs):

		self.ix = 0
		self.iy = 0
		self.iz = 0
		self.sx = 0
		self.sy = 0
		self.sz = 0
		self.ixsx = 0
		self.iysx = 0
		self.izsx = 0
		self.ixsy = 0
		self.iysy = 0
		self.izsy = 0	# -kba/(kab + kba)	# the initial state
		self.ixsz = 0
		self.iysz = 0
		self.izsz = 0

		self.bix = 0
		self.biy = 0
		self.biz = 0
		self.bsx = 0
		self.bsy = 0
		self.bsz = 0
		self.bixsx = 0
		self.biysx = 0
		self.bizsx = 0
		self.bixsy = 0
		self.biysy = 0
		self.bizsy = 0	# -kab/(kab + kba)	# the initial state
		self.bixsz = 0
		self.biysz = 0
		self.bizsz = 0

		for kw in kwargs:
			setattr(self, kw, kwargs[kw])

		self.vec = np.array([1/2, self.ix, self.iy, self.iz, self.sx, self.sy, self.sz, \
			self.ixsx, self.iysx, self.izsx, self.ixsy, self.iysy, self.izsy, self.ixsz, self.iysz, self.izsz, \
			self.bix, self.biy, self.biz, self.bsx, self.bsy, self.bsz, \
			self.bixsx, self.biysx, self.bizsx, self.bixsy, self.biysy, self.bizsy, self.bixsz, self.biysz, self.bizsz])

# ST-CW-CPMG y y x -x, Jiang et al. (2015) J. Magn. Reson. 257, 1.
def rcpulse_STCW_yyxx(de, lct, ofs, dw):

	k = 1
	if lct == 0:
		pcw = Ham.tex/4.0
	else:
		pcw = Ham.tex/lct/4.0		# the length for one 13C echo = the length of a 1H 360deg CW pulse
	while (pcw/k) > (1/Ham.minhcw):
		k += 1
	v1i = k/pcw
	# The above code can adjust the 1H cw power such that Vcw = 2n * Vcpmg
	# Activating the following line can fix the 1H cw power.
	# v1i = hamr.minhcw

	hamr = Ham(v1i = v1i, offs = ofs, dwsb = dw)
	hamrde = expm(hamr.ham * de)

	try:
		z1 = -hamr.kba/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z1 = -1
	try:
		z2 = -hamr.kab/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z2 = 0

	vec1 = Vec(sy = z1, bsy = z2)
	vec = np.copy(vec1.vec)

	# 1st CPMG
	ham1 = Ham(v1i = v1i, v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + 1 * pi/2, offs = ofs, dwsb = dw)
	ham1n180 = expm(ham1.ham * hamr.n90 *2)

	ham2 = Ham(v1i = v1i, v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + 0 * pi/2, offs = ofs, dwsb = dw)
	ham2n180 = expm(ham2.ham * hamr.n90 *2)

	ham3 = Ham(v1i = v1i, v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + 2 * pi/2, offs = ofs, dwsb = dw)
	ham3n180 = expm(ham3.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)		# y
		np.dot(hamrde, vec, vec)

		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)		# y
		np.dot(hamrde, vec, vec)

		np.dot(hamrde, vec, vec)
		np.dot(ham2n180, vec, vec)		# x
		np.dot(hamrde, vec, vec)

		np.dot(hamrde, vec, vec)
		np.dot(ham3n180, vec, vec)		# -x
		np.dot(hamrde, vec, vec)

	# {y, y}m (m=0 or 1) is omitted to match the number of 180deg pulses to those of the other sequences.

	return vec

# The in-phase CW 1H-decoupling pulse sequence
def rcpulse_inphase_CW(de, lct, p1, c1, ofs, dw):

	k = 1
	if lct == 0:
		pcw = Ham.tex/4.0
	else:
		pcw = Ham.tex/lct/4.0		# the length for one 13C echo = the length of a 1H 360deg CW pulse
	while (pcw/k) > (1/Ham.minhcw):
		k += 1
	v1i = k/pcw
	# The above code can adjust the 1H cw power such that Vcw = 2n * Vcpmg
	# Activating the following line can fix the 1H cw power.
	# v1i = hamr.minhcw

	hamr = Ham(v1i = v1i, offs = ofs, dwsb = dw)
	hamrde = expm(hamr.ham * de)

	try:
		z1 = -hamr.kba/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z1 = -1
	try:
		z2 = -hamr.kab/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z2 = 0

	vec1 = Vec(sy = z1, bsy = z2)
	vec = np.copy(vec1.vec)

	# 1st CPMG
	ham1 = Ham(v1i = v1i, v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c1 * pi/2, offs = ofs, dwsb = dw)
	ham1n180 = expm(ham1.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# 2nd CPMG

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# the central 180deg S pulse

	hamrde2 = expm(hamr.ham * 0.025/1000)
	ham9 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p1*pi/2 + pi, offs = ofs, dwsb = dw)
	ham9n60 = expm(ham9.ham * hamr.n90*60/90.0)
	ham10 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p1*pi/2, offs = ofs, dwsb = dw)
	ham10n300 = expm(ham10.ham * hamr.n90*300/90.0)

	np.dot(hamrde2, vec, vec)		# 25 us
	np.dot(ham9n60, vec, vec)		# 15N 60deg
	np.dot(ham10n300, vec, vec)		# 15N 300deg
	np.dot(ham9n60, vec, vec)		# 15N 60deg
	np.dot(hamrde2, vec, vec)		# 25 us

	# 3rd CPMG

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# 4th CPMG

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)
		np.dot(hamrde, vec, vec)

	return vec

# The conventional RC pulse sequence
def rcpulse_conventional(de, lct, p1, c1, c2, c3, c4, ofs, dw):

	hamr = Ham(offs = ofs, dwsb = dw)
	hamrde = expm(hamr.ham * de)

	try:
		z1 = -hamr.kba/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z1 = -1
	try:
		z2 = -hamr.kab/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z2 = 0

	vec1 = Vec(izsy = z1, bizsy = z2)

	# print('ID before = ', id(vec1.vec))
	vec = np.copy(vec1.vec)
	# print('ID after = ', id(vec))

	# 1st CPMG
	ham1 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c1 * pi/2, offs = ofs, dwsb = dw)
	ham1n180 = expm(ham1.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# 2nd CPMG
	ham2 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c2 * pi/2, offs = ofs, dwsb = dw)
	ham2n180 = expm(ham2.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham2n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# S: 60-300-60, I: 90-180-90
	if (hamr.n90 * 150/90 - hamr.h90 *2.0) <= 0:
		sys.exit()

	hamrj4 = expm(hamr.ham * hamr.dmis * (1.0/4/abs(hamr.jis)-hamr.n90 * 2.0/pi*2.2))
	ham3 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p1*pi/2 + pi, offs = ofs, dwsb = dw)
	ham3n60 = expm(ham3.ham * hamr.n90*60/90.0)
	ham4 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p1*pi/2, offs = ofs, dwsb = dw)
	ham4n150 = expm(ham4.ham * (hamr.n90*150/90.0-hamr.h90*2))
	ham5 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, v1i = 1/4/hamr.h90*hamr.hpmis, phs = hamr.phs + p1*pi/2, phi = 0, offs = ofs, dwsb = dw)
	ham5n150h90 = expm(ham5.ham * hamr.h90)
	ham6 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, v1i = 1/4/hamr.h90*hamr.hpmis, phs = hamr.phs + p1*pi/2, phi = pi/2, offs = ofs, dwsb = dw)
	ham6n300h180 = expm(ham6.ham * hamr.h90*2)

	np.dot(hamrj4, vec, vec)		# 1/(4j)
	np.dot(ham3n60, vec, vec)		# 15N 60deg
	np.dot(ham4n150, vec, vec)		# 15N 150
	np.dot(ham5n150h90, vec, vec)	# 1H 90
	np.dot(ham6n300h180, vec, vec)	# 15N 300 + 1H 180
	np.dot(ham5n150h90, vec, vec)	# 1H 90
	np.dot(ham4n150, vec, vec)		# 15N 150
	np.dot(ham3n60, vec, vec)		# 15N 60deg
	np.dot(hamrj4, vec, vec)		# 1/(4j)

	# 3rd CPMG
	ham7 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c3 * pi/2, offs = ofs, dwsb = dw)
	ham7n180 = expm(ham7.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham7n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# 4th CPMG
	ham8 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c4 * pi/2, offs = ofs, dwsb = dw)
	ham8n180 = expm(ham8.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham8n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# print('ID final = ', id(vec))
	return vec

# The new AFTAC pulse sequence
def rcpulse_AFTAC(de, lct, p1, p2, p3, c1, c2, c3, c4, ofs, dw):

	hamr = Ham(offs = ofs, dwsb = dw)
	hamrde = expm(hamr.ham * de)

	try:
		z1 = -hamr.kba/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z1 = -1
	try:
		z2 = -hamr.kab/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		z2 = 0

	vec1 = Vec(izsy = z1, bizsy = z2)

	vec = np.copy(vec1.vec)

	# 1st CPMG
	ham1 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c1 * pi/2, offs = ofs, dwsb = dw)
	ham1n180 = expm(ham1.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham1n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# S: 60-300-60, I: 90-180-90
	if (hamr.n90 * 150/90 - hamr.h90 *2.0) <= 0:
		sys.exit()

	hamrj4 = expm(hamr.ham * hamr.dmis * (1.0/4/abs(hamr.jis)-hamr.n90 * 2.0/pi*2.2))
	ham3 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p1*pi/2 + pi, offs = ofs, dwsb = dw)
	ham3n60 = expm(ham3.ham * hamr.n90*60/90.0)
	ham4 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p1*pi/2, offs = ofs, dwsb = dw)
	ham4n150 = expm(ham4.ham * (hamr.n90*150/90.0-hamr.h90*2))
	ham5 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, v1i = 1/4/hamr.h90*hamr.hpmis, phs = hamr.phs + p1*pi/2, phi = 0, offs = ofs, dwsb = dw)
	ham5n150h90 = expm(ham5.ham * hamr.h90)
	ham6 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, v1i = 1/4/hamr.h90*hamr.hpmis, phs = hamr.phs + p1*pi/2, phi = pi/2, offs = ofs, dwsb = dw)
	ham6n300h180 = expm(ham6.ham * hamr.h90*2)

	np.dot(hamrj4, vec, vec)		# 1/(4j)
	np.dot(ham3n60, vec, vec)		# 15N 60deg
	np.dot(ham4n150, vec, vec)		# 15N 150
	np.dot(ham5n150h90, vec, vec)	# 1H 90
	np.dot(ham6n300h180, vec, vec)	# 15N 300 + 1H 180
	np.dot(ham5n150h90, vec, vec)	# 1H 90
	np.dot(ham4n150, vec, vec)		# 15N 150
	np.dot(ham3n60, vec, vec)		# 15N 60deg
	np.dot(hamrj4, vec, vec)		# 1/(4j)

	# 2nd CPMG
	ham2 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c2 * pi/2, offs = ofs, dwsb = dw)
	ham2n180 = expm(ham2.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham2n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# the central 180deg S pulse

	hamrde2 = expm(hamr.ham * hamr.dmis * 1.22/1000)
	ham9 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p2*pi/2 + pi, offs = ofs, dwsb = dw)
	ham9n60 = expm(ham9.ham * hamr.n90*60/90.0)
	ham10 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p2*pi/2, offs = ofs, dwsb = dw)
	ham10n300 = expm(ham10.ham * hamr.n90*300/90.0)

	np.dot(hamrde2, vec, vec)		# 1.22 ms
	np.dot(ham9n60, vec, vec)		# 15N 60deg
	np.dot(ham10n300, vec, vec)		# 15N 300deg
	np.dot(ham9n60, vec, vec)		# 15N 60deg
	np.dot(hamrde2, vec, vec)		# 1.22 ms

	# 3rd CPMG
	ham7 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c3 * pi/2, offs = ofs, dwsb = dw)
	ham7n180 = expm(ham7.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham7n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# S: 60-300-60, I: 90-180-90
	if (hamr.n90 * 150/90 - hamr.h90 *2.0) <= 0:
		sys.exit()

	hamrj4 = expm(hamr.ham * hamr.dmis * (1.0/4/abs(hamr.jis)-hamr.n90 * 2.0/pi*2.2))
	ham11 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p3*pi/2 + pi, offs = ofs, dwsb = dw)
	ham11n60 = expm(ham11.ham * hamr.n90*60/90.0)
	ham12 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + p3*pi/2, offs = ofs, dwsb = dw)
	ham12n150 = expm(ham12.ham * (hamr.n90*150/90.0-hamr.h90*2))
	ham13 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, v1i = 1/4/hamr.h90*hamr.hpmis, phs = hamr.phs + p3*pi/2, phi = 0, offs = ofs, dwsb = dw)
	ham13n150h90 = expm(ham13.ham * hamr.h90)
	ham14 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, v1i = 1/4/hamr.h90*hamr.hpmis, phs = hamr.phs + p3*pi/2, phi = pi/2, offs = ofs, dwsb = dw)
	ham14n300h180 = expm(ham14.ham * hamr.h90*2)

	np.dot(hamrj4, vec, vec)		# 1/(4j)
	np.dot(ham11n60, vec, vec)		# 15N 60deg
	np.dot(ham12n150, vec, vec)		# 15N 150
	np.dot(ham13n150h90, vec, vec)	# 1H 90
	np.dot(ham14n300h180, vec, vec)	# 15N 300 + 1H 180
	np.dot(ham13n150h90, vec, vec)	# 1H 90
	np.dot(ham12n150, vec, vec)		# 15N 150
	np.dot(ham11n60, vec, vec)		# 15N 60deg
	np.dot(hamrj4, vec, vec)		# 1/(4j)

	# 4th CPMG
	ham8 = Ham(v1s = 1/4/hamr.n90*hamr.pmis, phs = hamr.phs + c4 * pi/2, offs = ofs, dwsb = dw)
	ham8n180 = expm(ham8.ham * hamr.n90 *2)

	for k in range(lct):
		np.dot(hamrde, vec, vec)
		np.dot(ham8n180, vec, vec)
		np.dot(hamrde, vec, vec)

	# print('ID final = ', id(vec))
	return vec

# Rex as a function of the Vcpmg
def rex_vs_vcpmg(hamr=Ham(), grp=1, flw=1, noise=0):

	magx = np.zeros((hamr.lctm+1, 2))

	np.random.seed()
	rec1 = 0
	rec2 = 0
	rec3 = 0
	rec4 = 0

	for lct in range(hamr.lctm+1):
		if lct == 0:
			d17 = 0.0
		else:
			d17 = hamr.tex/lct/8.0 - hamr.n90
			magx[lct][0] = 2*lct/hamr.tex

		if hamr.hpi == 6:
			vec1 = rcpulse_inphase_CW(d17, lct, 0, 1, hamr.offs, hamr.dwsb)		# 1H-CW in-phase cpmg
			rec1 = vec1[5]		# Sy

		elif hamr.hpi == 7:
			vec1 = rcpulse_conventional(d17, lct, 0, 1, 1, 0, 0, hamr.offs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
			rec1 = vec1[4]

		elif hamr.hpi == 4:
			vec1 = rcpulse_conventional(d17, lct, 0, 1, 1, 0, 0, hamr.offs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
			# rec1 = (vec1[4] - vec1[9])*0.5		# TROSY component
			rec1 = (vec1[4] + vec1[9])*0.5			# anti-TROSY component

		elif hamr.hpi == 3:		# The phase of the central 180deg pulse is shifted by 90deg.
			vec1 = rcpulse_conventional(d17, lct, 1, 1, 1, 0, 0, hamr.offs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
			# rec1 = (vec1[4] - vec1[9])*0.5		# TROSY component
			rec1 = (vec1[4] + vec1[9])*0.5			# anti-TROSY component

		elif hamr.hpi == 8:
			vec1 = rcpulse_AFTAC(d17, lct, 0, 0, 2, 1, 0, 2, 1, hamr.offs, hamr.dwsb)		# AFTAC
			rec1 = vec1[12]

		elif hamr.hpi == 5:
			vec1 = rcpulse_STCW_yyxx(d17, lct, hamr.offs, hamr.dwsb)		# ST-CW-CPMG
			rec1 = vec1[5]		# Sy

		if hamr.hpi == 6:
			vec2 = rcpulse_inphase_CW(d17, lct, 2, 1, hamr.offs, hamr.dwsb)		# 1H-CW in-phase cpmg
			rec2 = vec2[5]		# Sy

		elif hamr.hpi == 7:
			vec2 = rcpulse_conventional(d17, lct, 2, 1, 1, 0, 0, hamr.offs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
			rec2 = vec2[4]

		elif hamr.hpi == 4:
			vec2 = rcpulse_conventional(d17, lct, 2, 1, 1, 0, 0, hamr.offs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
			# rec2 = (vec2[4] - vec2[9])*0.5		# TROSY component
			rec2 = (vec2[4] + vec2[9])*0.5			# anti-TROSY component

		elif hamr.hpi == 3:		# The phase of the central 180deg pulse is shifted by 90deg.
			vec2 = rcpulse_conventional(d17, lct, 3, 1, 1, 0, 0, hamr.offs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
			# rec2 = (vec2[4] - vec2[9])*0.5		# TROSY component
			rec2 = (vec2[4] + vec2[9])*0.5			# anti-TROSY component

		elif hamr.hpi == 8:
			vec2 = rcpulse_AFTAC(d17, lct, 0, 2, 2, 1, 0, 2, 1, hamr.offs, hamr.dwsb)		# AFTAC
			rec2 = vec2[12]

		if hamr.hpi == 6:
			ns = 2
		elif hamr.hpi == 7 or hamr.hpi == 4 or hamr.hpi == 3:
			ns = 2
		elif hamr.hpi == 8:
			ns = 2
		elif hamr.hpi == 5:
			ns = 1
		else:
			ns = 1

		yrec = (rec1 + rec2 + rec3 + rec4)/ns

		if lct == 0:
			iyrec = yrec
		else:
			a = yrec/iyrec
			if a <= 0:
				a = 1.0e-20
				print("Warning: logarithm becomes infinity in rex_vs_vcpmg during synthesizing the RD curve: the data should be discarded.")
			magx[lct][1] = -np.log(a) / hamr.tex + noise * np.random.randn()

	magx0 = np.delete(magx, 0, axis=0)	# line 0 is deleted.

	if grp == 1:
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		# plt.rcParams["font.family"] = "Times New Roman"
		ax.scatter(magx0[:, 0], magx0[:, 1], s=16, color="blue", label="measured")
		ynew = interp1d(magx0[:, 0], magx0[:, 1], kind='cubic')
		xnew = np.linspace(min(magx0[:, 0]), max(magx0[:, 0]), num=100)
		ax.plot(xnew, ynew(xnew), color="magenta", label="fit", linestyle='dotted', linewidth = 2.0)
		ax.set_xlabel("CPMG frequency (/s)", fontsize=12)
		ax.set_ylabel("Rex (/s)", fontsize=12)
		ax.set_title('Rex as a function of Vcpmg')
		# ax.set_ylim([-1.0, 14.0])
		# plt.show()
		plt.savefig("testtest.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()
		# print(x)

	if flw == 1:
		with open("magpy.txt", "w") as fp:
			for i in range(len(magx0[:, 0])):
				print("{:.3f}   {:.4f}".format(magx0[i, 0], magx0[i, 1]), file=fp)

	return magx0

# Peak intensity as a function of the offset
def peak_vs_offset():

	hamr = Ham()
	magx  = np.zeros((hamr.gmax, 2))

	try:
		izsy = -hamr.kba/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		izsy = -1
	try:
		bizsy = -hamr.kab/(hamr.kab + hamr.kba)
	except ZeroDivisionError:
		bizsy = 0

	for h in range(4):
		if h == 0:
			lct = 0
		elif h == 1:
			lct = 1
		elif h == 2:
			lct = 2
		else:
			lct = 40

		if lct == 0:
			d17 = 0
		else:
			d17 = hamr.tex/lct/8.0 - hamr.n90

		for g in range(hamr.gmax):
			ofs = (g+1-hamr.gmax/2.0) * 2.0 * hamr.ofsrng / hamr.gmax
			rec1 = 0
			rec2 = 0
			rec3 = 0
			rec4 = 0

			if hamr.hpi == 6:
				vec1 = rcpulse_inphase_CW(d17, lct, 0, 1, ofs, hamr.dwsb)		# 1H-CW in-phase cpmg
				rec1 = vec1[5]		# Sy

			elif hamr.hpi == 7:
				vec1 = rcpulse_conventional(d17, lct, 0, 1, 1, 0, 0, ofs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
				rec1 = vec1[4]

			elif hamr.hpi == 4:
				vec1 = rcpulse_conventional(d17, lct, 0, 1, 1, 0, 0, ofs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
				rec1 = (vec1[4] + vec1[9])*0.5		# anti-TROSY component

			elif hamr.hpi == 3:
				vec1 = rcpulse_conventional(d17, lct, 1, 1, 1, 0, 0, ofs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
				rec1 = (vec1[4] + vec1[9])*0.5		# anti-TROSY component

			elif hamr.hpi == 8:
				vec1 = rcpulse_AFTAC(d17, lct, 0, 0, 2, 1, 0, 2, 1, ofs, hamr.dwsb)		# AFTAC, 60-300-60 pulse Best
				rec1 = vec1[12]			# state A
				# vec1 = rcpulse_AFTAC(d17, lct, 0, 0, 2, 3, 0, 2, 3, ofs, hamr.dwsb)		# 60-300-60 Varian
				# rec1 = vec1[27]		# +15, state B

			elif hamr.hpi == 5:
				vec1 = rcpulse_STCW_yyxx(d17, lct, ofs, hamr.dwsb)		# ST-CW-CPMG
				rec1 = vec1[5]		# Sy

			if hamr.hpi == 6:
				vec2 = rcpulse_inphase_CW(d17, lct, 2, 1, ofs, hamr.dwsb)		# 1H-CW in-phase cpmg
				rec2 = vec2[5]		# Sy

			elif hamr.hpi == 7:
				vec2 = rcpulse_conventional(d17, lct, 2, 1, 1, 0, 0, ofs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
				rec2 = vec2[4]

			elif hamr.hpi == 4:
				vec2 = rcpulse_conventional(d17, lct, 2, 1, 1, 0, 0, ofs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
				rec2 = (vec2[4] + vec2[9])*0.5		# anti-TROSY component

			elif hamr.hpi == 3:
				vec2 = rcpulse_conventional(d17, lct, 3, 1, 1, 0, 0, ofs, hamr.dwsb)	# 2*cpmg(2SyIz) - P - 2*cpmg(Sx)
				rec2 = (vec2[4] + vec2[9])*0.5		# anti-TROSY component

			elif hamr.hpi == 8:
				vec2 = rcpulse_AFTAC(d17, lct, 0, 2, 2, 1, 0, 2, 1, ofs, hamr.dwsb)		# AFTAC, 60-300-60 pulse Best
				rec2 = vec2[12]			# state A
				# vec2 = rcpulse_AFTAC(d17, lct, 0, 2, 2, 3, 0, 2, 3, ofs, hamr.dwsb)		# 60-300-60 Varian
				# rec2 = vec2[27]		# +15, state B

			if hamr.hpi == 6:
				ns = 2
			elif hamr.hpi == 7 or hamr.hpi == 4 or hamr.hpi == 3:
				ns = 2
			elif hamr.hpi == 8:
				ns = 2
			elif hamr.hpi == 5:
				ns = 1
			else:
				ns = 1

			magx[g][0] = ofs
			magx[g][1] = (rec1 + rec2 + rec3 + rec4)/ns
			# magx[g][1] = (rec1 + rec2 + rec3 + rec4)/ns/izsy		# state A

		if h == 0:
			magx0y = magx[:, 1].copy()
			magx0x = magx[:, 0].copy()
		elif h == 1:
			magx1y = magx[:, 1] / magx0y
		elif h == 2:
			magx2y = magx[:, 1] / magx0y
		else:
			magx3y = magx[:, 1] / magx0y

	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	# plt.rcParams["font.family"] = "Times New Roman"
	ax.scatter(magx0x, magx1y, s=6, color="magenta", label="n=1")
	ax.scatter(magx0x, magx2y, s=6, color="blue", label="n=2")
	ax.scatter(magx0x, magx3y, s=6, color="green", label="n=40")
	ax.set_xlabel("offset (Hz)", fontsize=12)
	ax.set_ylabel("peak intensity ratio", fontsize=12)
	ax.set_title('Peak intensity ratio as a function of offset')
	ax.set_ylim([-0.1, 1.1])
	plt.legend()
	plt.show()
	# plt.savefig("testtest.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()
	# print(x)
	return magx0x, magx1y, magx2y, magx3y

# Fitting to the direct CPMG data, instead of using the Carver-Richards Equation

def readin_dat(fnam):

	drex = []
	drey = []

	f = open(fnam, 'r')
	datalist = f.readlines()

	for dat in datalist:
		sdat = dat.split()

		if (len(sdat) == 0):	# only space, only tab, etc.
			continue
		elif ("#" in sdat[0]) or ("!" in sdat[0]):
			continue
		else:
			drex.append(float(sdat[0]))
			drey.append(float(sdat[1]))
	f.close()

	return drex, drey

def one_cpmg(lnu, lpa, lkex, lrx, ldw):

	ljc = 0			# J coupling (Hz)
	lpb = 1.0 - lpa
	dp = np.array([[2*pi*ljc*0.5j-lrx-lkex*lpb, lkex*lpa], [lkex*lpb, 2*pi*(ldw+0.5*ljc)*1j-lrx-lkex*lpa]])

	edp = expm(dp * 1.0/4.0/lnu)
	# edn = np.conjugate(edp)					# = expm(dn * dlt)

	edpn = np.dot(edp, np.conjugate(edp))		# = np.dot(edp, edn)
	ednp = np.conjugate(edpn)					# = np.dot(edn, edp)

	apn = edpn.copy()			# the matrix must be copied
	anp = ednp.copy()

	lnc = int(2.0*tcp*lnu + 0.5)

	for i in range(2, lnc+1):
		aa = np.dot(ednp, apn)
		apn = np.dot(edpn, anp)
		anp = aa.view()			# the address is copied

	apn = np.dot(apn, [lpa, lpb])
	anp = np.dot(anp, [lpa, lpb])

	pint = (0.5*(apn + anp)).real
	r = -np.log(pint[0]/lpa)/tcp
	return(r)				# the imaginary part is canceled.

# {ldlt - 180deg - ldlt}n1 - (180deg)p1 - {ldlt - 180deg - ldlt}n2 - (180deg)p2 -

def two_cpmg(n1, n2, n3, n4, p1, p2, p3, tcp, ldlt, lpa, lkex, lrx, ldw, jdlt, gdlt):

	ljc = 0			# J coupling (Hz)
	lpb = 1.0 - lpa
	dp = np.array([[2*pi*ljc*0.5j-lrx-lkex*lpb, lkex*lpa], [lkex*lpb, 2*pi*(ldw+0.5*ljc)*1j-lrx-lkex*lpa]])

	edp = expm(dp * ldlt)
	# edn = np.conjugate(edp)					# = expm(dn * dlt)
	jedp = expm(dp * jdlt)
	gedp = expm(dp * gdlt)

	edpn = np.dot(edp, np.conjugate(edp))		# = np.dot(edp, edn)
	ednp = np.conjugate(edpn)					# = np.dot(edn, edp)

	jedpn = np.dot(jedp, np.conjugate(jedp))
	jednp = np.conjugate(jedpn)
	gedpn = np.dot(gedp, np.conjugate(gedp))
	gednp = np.conjugate(gedpn)

	for k in range(2):		# k = 0; reference, k = 1; the actual intensity

		if (k == 0):
			ln1 = 0
			ln2 = 0
			ln3 = 0
			ln4 = 0
		else:
			ln1 = n1
			ln2 = n2
			ln3 = n3
			ln4 = n4

		aa = np.array([[1, 0], [0, 1]])			# the matrix must be copied
		bb = np.array([[1, 0], [0, 1]])

		i = 0

		if (i % 2):
			aa = np.dot(gednp, aa)
			bb = np.dot(gedpn, bb)
		else:
			aa = np.dot(gedpn, aa)
			bb = np.dot(gednp, bb)
		i += 1

		# cpmg sequence 1

		for j in range(ln1):
			if (i % 2):
				aa = np.dot(ednp, aa)
				bb = np.dot(edpn, bb)
			else:
				aa = np.dot(edpn, aa)
				bb = np.dot(ednp, bb)
			i += 1

		if (p1 == 1):	# 1/(2J) for AFTAC
			if (i % 2):
				aa = np.dot(jednp, aa)
				bb = np.dot(jedpn, bb)
			else:
				aa = np.dot(jedpn, aa)
				bb = np.dot(jednp, bb)
			i += 1
		elif (p1 == 2):		# A normal 180 deg pulse
			if (i % 2):
				aa = np.dot(gednp, aa)
				bb = np.dot(gedpn, bb)
			else:
				aa = np.dot(gedpn, aa)
				bb = np.dot(gednp, bb)
			i += 1

		# cpmg sequence 2

		for j in range(ln2):
			if (i % 2):
				aa = np.dot(ednp, aa)
				bb = np.dot(edpn, bb)
			else:
				aa = np.dot(edpn, aa)
				bb = np.dot(ednp, bb)
			i += 1

		if (p2 == 1):		# 1/(2J) for the conventional RC version
			if (i % 2):
				aa = np.dot(jednp, aa)
				bb = np.dot(jedpn, bb)
			else:
				aa = np.dot(jedpn, aa)
				bb = np.dot(jednp, bb)
			i += 1
		elif (p2 == 2):		# A normal 180 deg pulse with a pair of gradients for the new version, AFTAC
			if (i % 2):
				aa = np.dot(gednp, aa)
				bb = np.dot(gedpn, bb)
			else:
				aa = np.dot(gedpn, aa)
				bb = np.dot(gednp, bb)
			i += 1

		# cpmg sequence 3

		for j in range(ln3):
			if (i % 2):
				aa = np.dot(ednp, aa)
				bb = np.dot(edpn, bb)
			else:
				aa = np.dot(edpn, aa)
				bb = np.dot(ednp, bb)
			i += 1

		if (p3 == 1):	# 1/(2J) for AFTAC
			if (i % 2):
				aa = np.dot(jednp, aa)
				bb = np.dot(jedpn, bb)
			else:
				aa = np.dot(jedpn, aa)
				bb = np.dot(jednp, bb)
			i += 1
		elif (p3 == 2):		# A normal 180 deg pulse
			if (i % 2):
				aa = np.dot(gednp, aa)
				bb = np.dot(gedpn, bb)
			else:
				aa = np.dot(gedpn, aa)
				bb = np.dot(gednp, bb)
			i += 1

		# cpmg sequence 4

		for j in range(ln4):
			if (i % 2):
				aa = np.dot(ednp, aa)
				bb = np.dot(edpn, bb)
			else:
				aa = np.dot(edpn, aa)
				bb = np.dot(ednp, bb)
			i += 1

		if (i % 2):
			aa = np.dot(gednp, aa)
			bb = np.dot(gedpn, bb)
		else:
			aa = np.dot(gedpn, aa)
			bb = np.dot(gednp, bb)
		i += 1

		aa = np.dot(aa, [lpa, lpb])
		bb = np.dot(bb, [lpa, lpb])

		if (k == 0):
			pint0 = (0.5*(aa + bb)).real
		else:
			pint1 = (0.5*(aa + bb)).real

	#if pint1[0] <= 0:
	#	print(f" pint1:{pint1[0]:f}")
	#if pint0[0] <= 0:
	#	print(f" pint0:{pint0[0]:f}")
	#	print(lpa, lkex, lrx, ldw)

	a = pint1[0]/pint0[0]
	if a <= 0:
		print("Warning: logarithm becomes infinity in two-cpmg during fitting")
		return 1000
	else:
		return(-np.log(a)/tcp)				# the imaginary part is canceled.

def cpmg_plot(n, fnam, maxnc, pa, kex, rx, dw, tcp, sn, p1, p2, p3, jdlt, gdlt, drex, drey, grp=1, flw=0):

	x = list(range(maxnc))
	y = list(range(maxnc))

	for lnc in range(1, maxnc+1):    # 1, 2, 3, ...., maxnc
		# (ldlt - Pi - ldlt)*lnc*sn = tcp
		ldlt = tcp/lnc/2.0/sn
		lnu = 1.0/4.0/ldlt
		# y[lnc-1] = one_cpmg(lnu, pa, kex, rx, dw)
		y[lnc-1] = two_cpmg(n1=lnc, n2=lnc, n3=lnc, n4=lnc, p1=p1, p2=p2, p3=p3, tcp=tcp, ldlt=ldlt, lpa=pa, lkex=kex, lrx=rx, ldw=dw, jdlt=jdlt, gdlt=gdlt)
		x[lnc-1] = lnu

	if flw == 1:
		with open(fnam, "w") as fp: # change tabs in the following each line in FOR...
			for lnc in range(1, maxnc+1):    # 1, 2, 3, ...., maxnc
				print("{:.3f}   {:.4f}".format(x[lnc-1], y[lnc-1]), file=fp)

	if grp == 1:
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)

		# plt.rcParams["font.family"] = "Times New Roman"

		if n == 2:
			ax.scatter(drex, drey, s=18, color="blue", label="measured")

		ynew = interp1d(x, y, kind='cubic')
		xnew = np.linspace(min(x), max(x), num=100)
		ax.plot(xnew, ynew(xnew), color="magenta", label="fit", linestyle='dotted', linewidth = 2.0)
		ax.set_xlabel("CPMG frequency (/s)", fontsize=12)
		ax.set_ylabel("Rex (/s)", fontsize=12)
		plt.show()
		# plt.savefig("testtest.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()
		# print(x)

def calc_resd(prm, npx, npy, **kwargs):

	for kw in kwargs:
		if kw == 'tcp':
			tcp = kwargs[kw]
		elif kw == 'sn':
			sn = kwargs[kw]
		elif kw == 'p1':
			p1 = kwargs[kw]
		elif kw == 'p2':
			p2 = kwargs[kw]
		elif kw == 'p3':
			p3 = kwargs[kw]
		elif kw == 'jdlt':
			jdlt = kwargs[kw]
		elif kw == 'gdlt':
			gdlt = kwargs[kw]
		else:
			print('wrong1: ' + kw)
			sys.exit()

	resd = np.zeros_like(npx)

	for i in range(len(npx)):
		# resd[i] = npy[i] - one_cpmg(npx[i], prm[0], prm[1], prm[2], prm[3])
		ldlt = 1.0/4.0/npx[i]
		lnc = int(2.0*tcp*npx[i]/sn + 0.5)
		resd[i] = npy[i] - two_cpmg(n1=lnc, n2=lnc, n3=lnc, n4=lnc, p1=p1, p2=p2, p3=p3, tcp=tcp, ldlt=ldlt, lpa=prm[0], lkex=prm[1], lrx=prm[2], ldw=prm[3], jdlt=jdlt, gdlt=gdlt)
	return resd

def calc_resd_2B0(prm, npx, npy, **kwargs):

	for kw in kwargs:
		if kw == 'b0ratio':			# e.g. 800/500
			b0ratio = kwargs[kw]
		elif kw == 'tcp1':
			tcp1 = kwargs[kw]
		elif kw == 'tcp2':
			tcp2 = kwargs[kw]
		elif kw == 'sn':
			sn = kwargs[kw]
		elif kw == 'p1':
			p1 = kwargs[kw]
		elif kw == 'p2':
			p2 = kwargs[kw]
		elif kw == 'p3':
			p3 = kwargs[kw]
		elif kw == 'jdlt':
			jdlt = kwargs[kw]
		elif kw == 'gdlt':
			gdlt = kwargs[kw]
		elif kw == 'num1':			# The number of the elements in npx1 (e.g. 500 MHz)
			num1 = kwargs[kw]
		else:
			print('wrong2: ' + kw)
			sys.exit()

	resd = np.zeros_like(npx)

	for i in range(num1):
		ldlt = 1.0/4.0/npx[i]
		lnc = int(2.0*tcp1*npx[i]/sn + 0.5)
		resd[i] = npy[i] - two_cpmg(n1=lnc, n2=lnc, n3=lnc, n4=lnc, p1=p1, p2=p2, p3=p3, tcp=tcp1, ldlt=ldlt, lpa=prm[0], lkex=prm[1], lrx=prm[2], ldw=prm[3], jdlt=jdlt, gdlt=gdlt)

	for i in range(num1, len(npx)):
		ldlt = 1.0/4.0/npx[i]
		lnc = int(2.0*tcp2*npx[i]/sn + 0.5)
		resd[i] = npy[i] - two_cpmg(n1=lnc, n2=lnc, n3=lnc, n4=lnc, p1=p1, p2=p2, p3=p3, tcp=tcp2, ldlt=ldlt, lpa=prm[0], lkex=prm[1], lrx=prm[4], ldw=prm[3]*b0ratio, jdlt=jdlt, gdlt=gdlt)

	return resd

def calc_cost(prm, npx, npy, pnum, **kwargs):

	if pnum == 4:
		resd = calc_resd(prm, npx, npy, **kwargs)
	elif pnum == 5:
		resd = calc_resd_2B0(prm, npx, npy, **kwargs)
	else:
		print('calc_cost wrong' + pnum)
		sys.exit()
	return ((resd * resd).sum())/(len(npx) - pnum)		# reduced-X2 with 4 or 5 parameters
	return np.sqrt(((resd * resd).sum())/len(npx))		# rmsd instead of X2

def fit_cpmg(tcp, sn, p1, p2, p3, jdlt, gdlt, drex, drey, hamr=Ham()):

	maxpa = 1.0
	minpa = 0
	maxkex = 3000.0
	minkex = 0
	maxrx = 100.0
	minrx = -5.0
	maxdw = 300.0
	mindw = 0

	bd = ([minpa, minkex, minrx, mindw], [maxpa, maxkex, maxrx, maxdw])	# -np.inf
	x = np.array(drex)
	y = np.array(drey)
	kwargs={'tcp':tcp, 'sn':sn, 'p1':p1, 'p2':p2, 'p3':p3, 'jdlt':jdlt, 'gdlt':gdlt}

	# the initial values, these are to be optimized.
	np.random.seed()
	minrmsd = np.inf

	for pm in range(5):

		lkex = max(minkex+0.1, min(hamr.kab + hamr.kba + 10.0*np.random.randn(), maxkex-0.1))	# Hz
		lrx = max(minrx+0.1, min(hamr.r2s + 1.0*np.random.randn(), maxrx-0.1))					# intrinsic auto-relaxation (R0), Hz
		ldw = max(mindw+0.1, min(hamr.dwsb + 5.0*np.random.randn(), maxdw-0.1))				# Delta-omega, Hz (not rad/s)
		try:
			lpa = max(minpa+0.01, min(hamr.kba/(hamr.kab + hamr.kba) + 0.1*np.random.randn(), maxpa-0.01))
		except ZeroDivisionError:
			lpa = 0.999

		p0 = [lpa, lkex, lrx, ldw]

		res = least_squares(calc_resd, p0, loss='soft_l1', bounds=bd, args=(x, y), kwargs=kwargs)
		rmsd = calc_cost([res.x[0], res.x[1], res.x[2], res.x[3]], x, y, 4, **kwargs)
		print(f"Fit  pa {res.x[0]:.2f}    kex {res.x[1]:.2f} Hz   rx0 {res.x[2]:.2f} Hz    dw {res.x[3]:.2f} Hz    hpi {hamr.hpi:d}")
		print(f" Cost {res.cost:.10f}")
		print(f" Rmsd {rmsd:.10f}")
		if rmsd < minrmsd or pm == 0:
			pa, kex, rx, dw = res.x[0], res.x[1], res.x[2], res.x[3]
			minrmsd = rmsd

	print(f"*** Final fit  pa {pa:.2f}    kex {kex:.2f} Hz   rx0 {rx:.2f} Hz    dw {dw:.2f} Hz\n")
	return pa, kex, rx, dw, minrmsd

def fit_cpmg_2B0(tcp1, tcp2, num1, b0ratio, sn, p1, p2, p3, jdlt, gdlt, drex, drey):

	maxpa = 1.0
	minpa = 0
	maxkex = 3000.0
	minkex = 0
	maxrx = 100.0
	minrx = -5.0
	maxdw = 300.0
	mindw = 0

	bd = ([minpa, minkex, minrx, mindw, minrx], [maxpa, maxkex, maxrx, maxdw, maxrx])	# -np.inf
	x = np.array(drex)
	y = np.array(drey)
	kwargs={'tcp1':tcp1, 'tcp2':tcp2, 'b0ratio':b0ratio, 'sn':sn, 'p1':p1, 'p2':p2, 'p3':p3, 'jdlt':jdlt, 'gdlt':gdlt, 'num1':num1}

	# the initial values, these are to be optimized.
	np.random.seed()
	minrmsd = np.inf

	for pm in range(5):

		lkex = max(minkex+0.1, min(Ham.kab + Ham.kba + 10.0*np.random.randn(), maxkex-0.1))	# Hz
		lrx1 = max(minrx+0.1, min(Ham.r2s + 1.0*np.random.randn(), maxrx-0.1))					# intrinsic auto-relaxation (R0), Hz
		lrx2 = max(minrx+0.1, min(Ham.r2s + 1.0*np.random.randn(), maxrx-0.1))
		ldw = max(mindw+0.1, min(Ham.dwsb + 5.0*np.random.randn(), maxdw/b0ratio-0.1))				# Delta-omega, Hz (not rad/s)
		try:
			lpa = max(minpa+0.01, min(Ham.kba/(Ham.kab + Ham.kba) + 0.1*np.random.randn(), maxpa-0.01))
		except ZeroDivisionError:
			lpa = 0.999

		p0 = [lpa, lkex, lrx1, ldw, lrx2]

		res = least_squares(calc_resd_2B0, p0, loss='soft_l1', bounds=bd, args=(x, y), kwargs=kwargs)
		rmsd = calc_cost([res.x[0], res.x[1], res.x[2], res.x[3], res.x[4]], x, y, 5, **kwargs)
		print(f"Fit  pa {res.x[0]:.2f}    kex {res.x[1]:.2f} Hz   rx01 {res.x[2]:.2f} Hz  rx02 {res.x[4]:.2f} Hz   dw {res.x[3]:.2f} Hz    hpi {Ham.hpi:d}")
		print(f" Cost {res.cost:.10f}")
		print(f" Rmsd {rmsd:.10f}")
		if rmsd < minrmsd or pm == 0:
			pa, kex, rx1, dw, rx2 = res.x[0], res.x[1], res.x[2], res.x[3], res.x[4]
			minrmsd = rmsd

	print(f"*** Final fit  pa {pa:.2f}    kex {kex:.2f} Hz   rx01 {rx1:.2f} Hz  rx02 {rx2:.2f} Hz   dw {dw:.2f} Hz\n")
	return pa, kex, rx1, dw, rx2, minrmsd

def start_cpmg_sim_vcpmg():

	np.set_printoptions(precision=8, suppress=True)

	magx = rex_vs_vcpmg(grp=1, flw=1, noise=0)	# output to a graph (grp) and/or a file (flw)
	# print(magx)

def start_cpmg_sim_offset():

	np.set_printoptions(precision=4, suppress=True)

	magx0x, magx1y, magx2y, magx3y = peak_vs_offset()
	# print(magx1y)

def start_cpmg_fit(finnam='1x.inp', foutnam='1xout.txt'):

	if Ham.hpi == 6:	# in-phase 2*cpmg(Sx) - 180 - 2*cpmg(Sx) with 1H CW decoupling
		p1 = 0			# nothing: 0
		p2 = 2			# 180 deg
		p3 = 0			# nothing: 0
		gdlt = 25/1000000 # + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif Ham.hpi == 7 or Ham.hpi == 4 or Ham.hpi == 3:	# simple 2*cpmg(2SyIz) - Pelement - 2*cpmg(Sx), the conventional method
		p1 = 0			# nothing: 0
		p2 = 1			# 1/(2J), P-element: 1
		p3 = 0			# nothing: 0
		gdlt = 1.22/1000 # + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif Ham.hpi == 8:	# for the new version, AFTAC, 1-2-1
		p1 = 1			# 1/(2J), P-element: 1
		p2 = 2			# 180 deg with a pair of gradients: 2
		p3 = 1			# 1/(2J), P-element: 1
		gdlt = 1.22/1000 # + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif Ham.hpi == 5:	# ST-CW-CPMG
		p1 = 0			# nothing: 0
		p2 = 0			# nothing: 0
		p3 = 0			# nothing: 0
		gdlt = 0

	jdlt = 1.0/4/Ham.jis # - hamr.n90*2.0/pi*2.2 + hamr.n90*(60+300+60)/90/4	# 1/(4J)

	print("File name: " + finnam)

	drex, drey = readin_dat(finnam)
	maxnc = int(2.0*Ham.tex*max(drex)/Ham.sn + 0.5)
	# The max number of the spin-echo (ldlt - Pi - ldlt)lnc, L6
	pa, kex, rx, dw, rmsd = fit_cpmg(Ham.tex, Ham.sn, p1, p2, p3, jdlt, gdlt, drex, drey)
	# pa=0.30;    kex=497.34;   rx=3.35;    dw=85.49
	# pa=0.68;    kex=496.15;   rx=2.95;    dw=84.36
	# pa=0.50;    kex=496.15;   rx=2.95;    dw=84.36

	cpmg_plot(2, foutnam, maxnc, pa, kex, rx, dw, Ham.tex, Ham.sn, p1, p2, p3, jdlt, gdlt, drex, drey, grp=1, flw=0)	# cpmg_plot(1, 'sim_cpmg_2.txt', 24)	# Plot the cpmg curve, 2: scatter plot

def start_cpmg_sim_fit(grp=1, noise=0):

	magx = rex_vs_vcpmg(grp=0, flw=0, noise=noise)	# output to a graph (grp) and/or a file (flw)

	if Ham.hpi == 6:		# in-phase 2*cpmg(Sx) - 180 - 2*cpmg(Sx) with 1H CW decoupling
		p1 = 0				# nothing: 0
		p2 = 2				# 180 deg
		p3 = 0				# nothing: 0
		gdlt = 25/1000000	# + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif Ham.hpi == 7 or Ham.hpi == 4 or Ham.hpi == 3:		# simple 2*cpmg(2SyIz) - Pelement - 2*cpmg(Sx), the conventional method
		p1 = 0				# nothing: 0
		p2 = 1				# 1/(2J), P-element: 1
		p3 = 0				# nothing: 0
		gdlt = 1.22/1000	# + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif Ham.hpi == 8:		# for the new version, AFTAC, 1-2-1
		p1 = 1				# 1/(2J), P-element: 1
		p2 = 2				# 180 deg with a pair of gradients: 2
		p3 = 1				# 1/(2J), P-element: 1
		gdlt = 1.22/1000	# + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif Ham.hpi == 5:		# ST-CW-CPMG {y y x -x}
		p1 = 0				# nothing: 0
		p2 = 0				# nothing: 0
		p3 = 0				# nothing: 0
		gdlt = 0

	jdlt = 1.0/4/Ham.jis	# - hamr.n90*2.0/pi*2.2 # + hamr.n90*(60+300+60)/90/4		# 1/(4J)

	drex, drey = magx[:, 0], magx[:, 1]

	try:
		pa = Ham.kba/(Ham.kab + Ham.kba)
	except ZeroDivisionError:
		pa = 1
	kex = Ham.kab+Ham.kba
	rx = Ham.r2s
	dw = Ham.dwsb

	print(f"Sim  pa {pa:.2f}    kex {kex:.2f} Hz   rx0 {rx:.2f} Hz    dw {dw:.2f} Hz   offs {Ham.offs:.2f} Hz   pmis {Ham.pmis:.2f}   hpi {Ham.hpi:d}")

	maxnc = int(2.0*Ham.tex*max(drex)/Ham.sn + 0.5)
	# The max number of the spin-echo (ldlt - Pi - ldlt)lnc, L6
	pa, kex, rx, dw, rmsd = fit_cpmg(Ham.tex, Ham.sn, p1, p2, p3, jdlt, gdlt, drex, drey)

	if grp == 1:
		cpmg_plot(2, 'mag_out.txt', maxnc, pa, kex, rx, dw, Ham.tex, Ham.sn, p1, p2, p3, jdlt, gdlt, drex, drey, grp=1, flw=0)
	else:
		return pa, kex, rx, dw, rmsd, drex, drey

def start_cpmg_sim_fit_2B0(b01, b02, grp=1, noise=0):

	b0ratio = b02/b01
	hamr1 = Ham()
	hamr2 = Ham(dwsb = hamr1.dwsb * b0ratio, offs = hamr1.offs * b0ratio)

	magx1 = rex_vs_vcpmg(hamr=hamr1, grp=0, flw=0, noise=noise)	# output to a graph (grp) and/or a file (flw)
	magx2 = rex_vs_vcpmg(hamr=hamr2, grp=0, flw=0, noise=noise)
	magx = np.concatenate([magx1, magx2])

	if hamr1.hpi == 6:		# in-phase 2*cpmg(Sx) - 180 - 2*cpmg(Sx) with 1H CW decoupling
		p1 = 0				# nothing: 0
		p2 = 2				# 180 deg
		p3 = 0				# nothing: 0
		gdlt = 25/1000000	# + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif hamr1.hpi == 7 or hamr1.hpi == 4 or hamr1.hpi == 3:	# simple 2*cpmg(2SyIz) - Pelement - 2*cpmg(Sx), the conventional method
		p1 = 0				# nothing: 0
		p2 = 1				# 1/(2J), P-element: 1
		p3 = 0				# nothing: 0
		gdlt = 1.22/1000	# + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif hamr1.hpi == 8:	# for the new version, AFTAC, 1-2-1
		p1 = 1				# 1/(2J), P-element: 1
		p2 = 2				# 180 deg with a pair of gradients: 2
		p3 = 1				# 1/(2J), P-element: 1
		gdlt = 1.22/1000	# + hamr.n90*(60+300+60)/90/4		# each delay sandwitching the central 13C 180deg pulse (sec), 2.2 for 60-300-60 pulses

	elif hamr1.hpi == 5:	# ST-CW-CPMG {y y x -x}
		p1 = 0				# nothing: 0
		p2 = 0				# nothing: 0
		p3 = 0				# nothing: 0
		gdlt = 0

	jdlt = 1.0/4/hamr1.jis	# - hamr.n90*2.0/pi*2.2 # + hamr.n90*(60+300+60)/90/4		# 1/(4J)

	drex, drey = magx[:, 0], magx[:, 1]

	try:
		pa = hamr1.kba/(hamr1.kab + hamr1.kba)
	except ZeroDivisionError:
		pa = 1
	kex = hamr1.kab+hamr1.kba
	rx1 = hamr1.r2s
	dw = hamr1.dwsb

	print(f"Sim  pa {pa:.2f}    kex {kex:.2f} Hz   rx0 {rx1:.2f} Hz    dw {dw:.2f} Hz   offs {hamr1.offs:.2f} Hz   pmis {hamr1.pmis:.2f}   hpi {hamr1.hpi:d}")

	# ??? maxnc = int(2.0*hamr1.tex*max(drex)/hamr1.sn + 0.5)
	# The max number of the spin-echo (ldlt - Pi - ldlt)lnc, L6

	pa, kex, rx1, dw, rx2, rmsd = fit_cpmg_2B0(hamr1.tex, hamr2.tex, hamr1.lctm, b0ratio, hamr1.sn, p1, p2, p3, jdlt, gdlt, drex, drey)

	if grp == 1:
		cpmg_plot(2, 'mag_out1.txt', hamr1.lctm, pa, kex, rx1, dw, hamr1.tex, hamr1.sn, p1, p2, p3, jdlt, gdlt, drex[:hamr1.lctm], drey[:hamr1.lctm], grp=1, flw=0)
		cpmg_plot(2, 'mag_out2.txt', hamr2.lctm, pa, kex, rx2, dw*b0ratio, hamr2.tex, hamr2.sn, p1, p2, p3, jdlt, gdlt, drex[hamr1.lctm:], drey[hamr1.lctm:], grp=1, flw=0)
	else:
		return pa, kex, rx1, dw, rx2, rmsd, drex, drey

def start_cpmg_sim_param_cal_1B0():

	gmax = Ham.gmax
	ofsrng = Ham.ofsrng
	ini_offs = Ham.offs

	with open("cpmg_vs_offs_v.txt", "w") as fpv, open("cpmg_vs_offs_p.txt", "w") as fpp:
		print("# pm  pa  kex  rx  dw  rmsd", file=fpp)

		for g in range(gmax):
			pm = (g+1-gmax/2.0) * 2.0 * ofsrng / gmax
			Ham.offs = pm
			pa, kex, rx, dw, rmsd, drex, drey = start_cpmg_sim_fit(grp=0, noise=0)
			print("pm = {:.2f}".format(pm), file=fpv)
			for i in range(len(drex)):
				print("{:.2f}  {:.3f}  {:.3f}".format(pm, drex[i], drey[i]), file=fpv)
			print("mp", file=fpv)
			print("{:.2f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}".format(pm, pa, kex, rx, dw, rmsd), file=fpp)

	Ham.offs = ini_offs

def cpmg_rmsd_noexchange_plot_2d():

	f2 = open("cpmg_vs_offs_v2.txt", 'r')
	datalist = f2.readlines()

	drep2 = []
	drev2 = []
	rmsd = []
	aved = []
	npm2 = 0

	for dat in datalist:
		sdat = dat.split()

		if len(sdat) == 0:	# only space, only tab, etc.
			continue
		elif ("#" in sdat[0]) or ("!" in sdat[0]):
			continue
		else:
			pass

		if sdat[0] == 'pm':
			drep2.append(float(sdat[2]))
			drey2 = []
			npm2 += 1
		elif (npm2 == 1) and (sdat[0] != 'mp'):
			pass		# drev2.append(float(sdat[1]))
		else:
			pass

		if sdat[0] == 'mp':
			rmsd.append(np.sqrt((((drey2-np.average(drey2)) * (drey2-np.average(drey2))).sum())/len(drey2)))
			aved.append(np.average(drey2))
		elif sdat[0] != 'pm':
			drey2.append(float(sdat[2]))
		else:
			pass
	f2.close()

	fig = plt.figure(dpi=120)

	ax1 = fig.add_subplot(3, 2, 1)
	ax1.scatter(drep2, aved, s=4, color="cyan", label="measured")
	ax1.set_ylabel("average (/s)", fontsize=12, fontstyle='italic')
	ax1.set_xlabel("offset (Hz)", fontsize=12, fontstyle='italic')
	# ax1.set_ylim([-0.1, 8.0])

	ax2 = fig.add_subplot(3, 2, 2)
	ax2.scatter(drep2, rmsd, s=4, color="magenta", label="measured")
	ax2.set_ylabel("rmsd (/s)", fontsize=12, fontstyle='italic')
	ax2.set_xlabel("offset (Hz)", fontsize=12, fontstyle='italic')
	# ax2.set_ylim([-0.1, 8.0])

	# plt.show()
	plt.tight_layout()
	plt.savefig("test2d.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()
	# print(x)

def start_cpmg_sim_param_plot_2d():

	f = open("cpmg_vs_offs_p.txt", 'r')
	datalist = f.readlines()
	f.close()

	drex = []
	drep = []
	drek = []
	drer = []
	drew = []
	dred = []

	for dat in datalist:
		sdat = dat.split()

		if (len(sdat) == 0):	# only space, only tab, etc.
			continue
		elif ("#" in sdat[0]) or ("!" in sdat[0]):
			continue
		else:
			drex.append(float(sdat[0]))
			drep.append(float(sdat[1]))
			drek.append(float(sdat[2]))
			drer.append(float(sdat[3]))
			drew.append(float(sdat[4]))
			dred.append(float(sdat[5]))

	fig = plt.figure(dpi=120)
	ax1 = fig.add_subplot(3, 2, 1)
	ax2 = fig.add_subplot(3, 2, 2)
	ax3 = fig.add_subplot(3, 2, 3)
	ax4 = fig.add_subplot(3, 2, 4)
	ax5 = fig.add_subplot(3, 2, 5)

	# plt.rcParams["font.family"] = "Times New Roman"

	ax1.scatter(drex, drep, s=4, color="blue", label="measured")
	ynew = interp1d(drex, drep, kind='cubic')		# more than 4 points are required.
	xnew = np.linspace(min(drex), max(drex), num=100)
	# ax1.plot(xnew, ynew(xnew), color="grey", label="fit", linestyle='dotted', linewidth = 2.0)
	ax1.set_ylabel("p$_{a}$", fontsize=12, fontstyle='italic')
	ax1.set_ylim([-0.1, 1.1])

	ax2.scatter(drex, drek, s=4, color="red", label="measured")
	ynew = interp1d(drex, drek, kind='cubic')
	xnew = np.linspace(min(drex), max(drex), num=100)
	# ax2.plot(xnew, ynew(xnew), color="grey", label="fit", linestyle='dotted', linewidth = 2.0)
	ax2.yaxis.set_label_position("right")
	ax2.yaxis.tick_right()
	ax2.set_ylabel("k$_{ex}$ (/s)", fontsize=12, fontstyle='italic')
	ax2.set_ylim([0, 1500])

	ax3.scatter(drex, drer, s=4, color="cyan", label="measured")
	ynew = interp1d(drex, drer, kind='cubic')
	xnew = np.linspace(min(drex), max(drex), num=100)
	# ax3.plot(xnew, ynew(xnew), color="grey", label="fit", linestyle='dotted', linewidth = 2.0)
	ax3.set_xlabel("offset (Hz)", fontsize=12, fontstyle='italic')
	ax3.set_ylabel("R$_{0}$ (/s)", fontsize=12, fontstyle='italic')
	ax3.set_ylim([-0.1, 10.1])

	ax4.scatter(drex, drew, s=4, color="magenta", label="measured")
	ynew = interp1d(drex, drew, kind='cubic')
	xnew = np.linspace(min(drex), max(drex), num=100)
	# ax4.plot(xnew, ynew(xnew), color="grey", label="fit", linestyle='dotted', linewidth = 2.0)
	ax4.set_xlabel("offset (Hz)", fontsize=12, fontstyle='italic')
	ax4.yaxis.set_label_position("right")
	ax4.yaxis.tick_right()
	ax4.set_ylabel("Δω (Hz)", fontsize=12, fontstyle='italic')
	ax4.set_ylim([0, 200])

	ax5.scatter(drex, dred, s=4, color="green", label="measured")
	ynew = interp1d(drex, dred, kind='cubic')
	xnew = np.linspace(min(drex), max(drex), num=100)
	# ax5.plot(xnew, ynew(xnew), color="grey", label="fit", linestyle='dotted', linewidth = 2.0)
	ax5.set_xlabel("offset (Hz)", fontsize=12, fontstyle='italic')
	ax5.set_ylabel("rmsd (/s)", fontsize=12, fontstyle='italic')
	# ax5.set_ylim([-0.1, 10.1])

	# plt.show()
	plt.tight_layout()
	plt.savefig("test2d.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()
	# print(x)

def start_cpmg_sim_param_plot_3d():

	f = open("cpmg_vs_offs_v.txt", 'r')
	datalist = f.readlines()

	drep = []
	drev = []
	drer = []
	npm = 0

	for dat in datalist:
		sdat = dat.split()

		if len(sdat) == 0:	# only space, only tab, etc.
			continue
		elif ("#" in sdat[0]) or ("!" in sdat[0]):
			continue
		else:
			pass

		if sdat[0] == 'pm':
			drep.append(float(sdat[2]))
			drey = []
			npm += 1
		elif (npm == 1) and (sdat[0] != 'mp'):
			drev.append(float(sdat[1]))
		else:
			pass

		if sdat[0] == 'mp':
			drer.append(drey)
		elif sdat[0] != 'pm':
			drey.append(float(sdat[2]))
		else:
			pass

	f.close()

	x, y = np.meshgrid(drev, drep)
	z = np.array(drer)

	fig = plt.figure(dpi=120)
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	# fig.subplots_adjust(right=2)
	# plt.rcParams["font.family"] = "Times New Roman"
	# plt.rcParams['axes.labelpad'] = 16
	# plt.subplots_adjust(right=1.3)
	# plt.rcParams['figure.subplot.right'] = 0.48
	ax.plot_wireframe(x, y, z, rstride=1, cstride=0, color='blue', linewidth=0.7)
	# ax.plot_surface(x, y, z)
	# ax.set_box_aspect((1.2, 1, 1))
	# fontsize=10
	ax.set_xlabel("$ν_{cpmg}$ (/s)", fontstyle='italic')
	ax.set_zlabel("R$_{2}^{eff}$ (/s)", fontstyle='italic')
	ax.set_ylabel("offset (Hz)", fontstyle='italic')
	ax.set_title('R$_{2}^{eff}$ as a function of ν$_{cpmg}$ and offset', fontstyle='italic')
	ax.set_zlim([-3, 80])
	# plt.show()
	plt.tight_layout()
	plt.savefig("test3d.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()

def start_cpmg_sim_param_plot_3d_2B0():

	f1 = open("cpmg_vs_offs_v1.txt", 'r')
	datalist = f1.readlines()

	drep1 = []
	drev1 = []
	drer1 = []
	npm1 = 0

	for dat in datalist:
		sdat = dat.split()

		if len(sdat) == 0:	# only space, only tab, etc.
			continue
		elif ("#" in sdat[0]) or ("!" in sdat[0]):
			continue
		else:
			pass

		if sdat[0] == 'pm':
			drep1.append(float(sdat[2]))
			drey1 = []
			npm1 += 1
		elif (npm1 == 1) and (sdat[0] != 'mp'):
			drev1.append(float(sdat[1]))
		else:
			pass

		if sdat[0] == 'mp':
			drer1.append(drey1)
		elif sdat[0] != 'pm':
			drey1.append(float(sdat[2]))
		else:
			pass
	f1.close()

	f2 = open("cpmg_vs_offs_v2.txt", 'r')
	datalist = f2.readlines()

	drep2 = []
	drev2 = []
	drer2 = []
	npm2 = 0

	for dat in datalist:
		sdat = dat.split()

		if len(sdat) == 0:	# only space, only tab, etc.
			continue
		elif ("#" in sdat[0]) or ("!" in sdat[0]):
			continue
		else:
			pass

		if sdat[0] == 'pm':
			drep2.append(float(sdat[2]))
			drey2 = []
			npm2 += 1
		elif (npm2 == 1) and (sdat[0] != 'mp'):
			drev2.append(float(sdat[1]))
		else:
			pass

		if sdat[0] == 'mp':
			drer2.append(drey2)
		elif sdat[0] != 'pm':
			drey2.append(float(sdat[2]))
		else:
			pass
	f2.close()

	x1, y1 = np.meshgrid(drev1, drep1)
	z1 = np.array(drer1)

	x2, y2 = np.meshgrid(drev2, drep2)
	z2 = np.array(drer2)

	fig = plt.figure(dpi=120)

	ax1 = fig.add_subplot(1, 2, 1, projection='3d')
	ax2 = fig.add_subplot(1, 2, 2, projection='3d')

	# plt.rcParams["font.family"] = "Times New Roman"

	ax1.plot_wireframe(x1, y1, z1, rstride=1, cstride=0, color='blue', linewidth=0.7)
	# ax1.plot_surface(x1, y1, z1)
	ax2.plot_wireframe(x2, y2, z2, rstride=1, cstride=0, color='magenta', linewidth=0.7)
	# ax2.plot_surface(x2, y2, z2)

	ax1.set_xlabel("$ν_{cpmg}$ (/s)", fontsize=10, fontstyle='italic')
	ax1.set_zlabel("R$_{2}^{eff}$ (/s)", fontsize=10, fontstyle='italic')
	ax1.set_ylabel("offset (Hz)", fontsize=10, fontstyle='italic')
	ax1.set_title('R$_{2}^{eff}$ as a function of ν$_{cpmg}$ and offset', fontstyle='italic')
	# ax1.set_zlim([-3, 20])

	ax2.set_xlabel("$ν_{cpmg}$ (/s)", fontsize=10, fontstyle='italic')
	ax2.set_zlabel("R$_{2}^{eff}$ (/s)", fontsize=10, fontstyle='italic')
	ax2.set_ylabel("offset (Hz)", fontsize=10, fontstyle='italic')
	# ax2.set_title('R$_{2}^{eff}$ as a function of ν$_{cpmg}$ and offset', fontstyle='italic')
	# ax2.set_zlim([-3, 35])

	# plt.show()
	plt.tight_layout()
	plt.savefig("test3d.pdf", dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()

def start_cpmg_plot_2d(finnam='mag.dat', foutnam="test.pdf"):

	drex, drey = readin_dat(finnam)

	rmsd = np.sqrt((((drey-np.average(drey)) * (drey-np.average(drey))).sum())/len(drey))
	print(f"AVE: {np.average(drey):.3f}  RMSD: {rmsd:.3f}")

	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	# plt.rcParams["font.family"] = "Times New Roman"
	ax.scatter(drex, drey, s=18, color="blue", label="measured")
	ynew = interp1d(drex, drey, kind='cubic')
	xnew = np.linspace(min(drex), max(drex), num=100)
	ax.plot(xnew, ynew(xnew), color="magenta", label="fit", linestyle='dotted', linewidth = 2.0)
	ax.set_xlabel("$ν_{cpmg}$ (/s)", fontsize=12, fontstyle='italic')
	ax.set_ylabel("R$_{2}^{eff}$ (/s)", fontsize=12, fontstyle='italic')
	ax.set_title('R$_{2}^{eff}$ as a function of ν$_{cpmg}$', fontstyle='italic')
	ax.set_ylim([-1.0, 8.0])
	# plt.show()
	plt.savefig(foutnam, dpi=1000, bbox_inches="tight", transparent=True)		# must delete plt.show()
	# print(x)

def start_cpmg_sim_param_cal_2B0(b01, b02):	# b01: 500, b02: 800

	b0ratio = b02/b01
	ini_offs = Ham.offs

	with open("cpmg_vs_offs_v1.txt", "w") as fpv1, open("cpmg_vs_offs_v2.txt", "w") as fpv2, open("cpmg_vs_offs_p.txt", "w") as fpp:
		print("# pm  pa  kex  rx1  dw  rx2  rmsd", file=fpp)

		for g in range(Ham.gmax):
			pm = (g+1-Ham.gmax/2.0) * 2.0 * Ham.ofsrng / Ham.gmax
			Ham.offs = pm
			pa, kex, rx1, dw, rx2, rmsd, drex, drey = start_cpmg_sim_fit_2B0(b01, b02, grp=0, noise=0)

			print("pm = {:.2f}".format(pm), file=fpv1)
			for i in range(0, Ham.lctm):
				print("{:.2f}  {:.3f}  {:.3f}".format(pm, drex[i], drey[i]), file=fpv1)
			print("mp", file=fpv1)

			print("pm = {:.2f}".format(pm*b0ratio), file=fpv2)
			for i in range(Ham.lctm, Ham.lctm *2):
				print("{:.2f}  {:.3f}  {:.3f}".format(pm*b0ratio, drex[i], drey[i]), file=fpv2)
			print("mp", file=fpv2)

			print("{:.2f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}".format(pm, pa, kex, rx1, dw, rx2, rmsd), file=fpp)

	Ham.offs = ini_offs

def start_cpmg_sim_param_cal_fullauto_B():

	ini_hpi = Ham.hpi

	for hpi in [3, 4, 5, 6, 7, 8]:

		Ham.hpi = hpi

		dir = "hpi_" + str(hpi) + "_1B"			# with exchange, 1 B0
		os.makedirs(dir, exist_ok=True)
		os.chdir(dir)
		start_cpmg_sim_param_cal_1B0()
		start_cpmg_sim_param_plot_2d()
		start_cpmg_sim_param_plot_3d()
		shutil.move('test2d.pdf', dir + '_test2d.pdf')
		shutil.move('test3d.pdf', dir + '_test3d.pdf')
		plt.clf()
		plt.close("all")
		os.chdir('..')

		dir = "hpi_" + str(hpi) + "_2B"			# with exchange, 2 B0
		os.makedirs(dir, exist_ok=True)
		os.chdir(dir)
		start_cpmg_sim_param_cal_2B0(500, 800)
		start_cpmg_sim_param_plot_2d()
		start_cpmg_sim_param_plot_3d_2B0()
		shutil.move('test2d.pdf', dir + '_test2d.pdf')
		shutil.move('test3d.pdf', dir + '_test3d.pdf')
		plt.clf()
		plt.close("all")
		os.chdir('..')

		ini_kab = Ham.kab
		ini_kba = Ham.kba
		ini_r2s = Ham.r2s

		Ham.kab = 0
		Ham.kba = 0
		Ham.r2s = 0

		# dir = "hpi_" + str(hpi) + "_1B" + "_0"		# no exchange, 1 B0
		# os.makedirs(dir, exist_ok=True)
		# os.chdir(dir)
		# start_cpmg_sim_param_cal_1B0()
		# start_cpmg_sim_param_plot_2d()
		# start_cpmg_sim_param_plot_3d()
		# shutil.move('test2d.pdf', dir + '_test2d.pdf')
		# shutil.move('test3d.pdf', dir + '_test3d.pdf')
		# plt.clf()
		# plt.close("all")
		# os.chdir('..')

		dir = "hpi_" + str(hpi) + "_2B" + "_0"		# no exchange, 2 B0
		os.makedirs(dir, exist_ok=True)
		os.chdir(dir)
		start_cpmg_sim_param_cal_2B0(500, 800)
		cpmg_rmsd_noexchange_plot_2d()			# Rmsd from the average will be calculated only for v2
		start_cpmg_sim_param_plot_3d_2B0()		# Only v2 will be used.
		shutil.move('test2d.pdf', dir + '_test2d.pdf')
		shutil.move('test3d.pdf', dir + '_test3d.pdf')
		plt.clf()
		plt.close("all")
		os.chdir('..')

		Ham.kab = ini_kab
		Ham.kba = ini_kba
		Ham.r2s = ini_r2s

	Ham.hpi = ini_hpi

'''
reset -f
run cpmg_hsqc_exch_fit_29.py
'''

#start_cpmg_fit(finnam='1x.inp', foutnam='1xout.txt')	# fitting to an actual experimental RD curve

#start_cpmg_sim_vcpmg()									# simulation of an RD curve against Vcpmg
#start_cpmg_sim_offset()								# simulation of peak intensity ratios against offset

#start_cpmg_sim_fit()									# simulation and fitting on a single B0 static magnetic field
#start_cpmg_sim_fit_2B0(500, 800)						# simulation and fitting on two B0 static magnetic fields, 500: primary, 800: secondary

#start_cpmg_sim_param_cal_1B0()							# simulation of a series of RD curves against offset on a single B0 static magnetic field
#start_cpmg_sim_param_plot_2d()							# 2D plot of pa, Kex, R0, dW
#start_cpmg_sim_param_plot_3d()							# 3D plot of RD curves against offset on a single B0 static magnetic field

#start_cpmg_sim_param_cal_2B0(500, 800)					# simulation of a series of RD curves against offset on two B0 static magnetic fields
#start_cpmg_sim_param_plot_2d()							# 2D plot of pa, Kex, R0, dW
#start_cpmg_sim_param_plot_3d_2B0()						# 3D plot of RD curves against offset on two B0 static magnetic fields

#start_cpmg_sim_param_cal_fullauto_B()					# This is only for automatically creating a series of data figures to be submitted to a paper

#start_cpmg_plot_2d(finnam='mag.dat', foutnam="test.pdf")	# 2D plot of an actual RD curve