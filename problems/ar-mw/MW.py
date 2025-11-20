import numpy as np
from scipy.optimize import least_squares
from mdcraft.solver.thermostat import Langevin, Compound as CompoundThermostat

def linear(coeffs, xtrain, ytrain):
	A, B = coeffs
	return np.sum((A*xtrain + B - ytrain)**2)

class MovingWindow:

	def __init__(self, buffer, bxmin, bxmax, xin, xout, xfront, uin, uout = 0.0, bpos = 0.0):
		self.buffer  = buffer # atoms buffer to insert
		self.xin     = xin    # position of inflow atoms
		self.xout    = xout   # position of outflow atoms
		self.xfront  = xout   # current position of the shock front
		self.xfront0 = xfront # desired position of the shock front
		self.vshock  = 0.0    # shock wave speed
		self.uin     = uin    # inflow atoms velocity
		self.uout    = uout   # outflow atoms velocity
		# shock velocity estimation history
		self.history = {"xfront": np.zeros(10), "time": np.zeros(10)} 
		self.hindex  = 0
		self.hready  = False  # flag which shows that history is full
		# shift coordinates to 0
		self.buffer["r"]["x"] -= bxmin
		# initial position in buffer
		self.bpos = bpos
		# period of the buffer
		self.period = bxmax - bxmin

	def thermostats(self, Tin, Tout, dt, \
		betaInTight = 5,  betaInSoft = 0.5, betaOut = 5, \
		widthInTight = 5, widthInSoft = 10, widthOut = 5, \
		frozenWidth = 2):

		self.frozenWidth = frozenWidth
		thermostats = CompoundThermostat()
		self.thermostatInTight = Langevin(
		    beta        = betaInTight, # ps^-1 
		    temperature = Tin, # K 
		    time_step   = dt,  # ps
		    heat_x      = 1,
		    heat_y      = 1,
		    heat_z      = 1,
		    Ux          = self.uin,
		    Uy          = 0,
		    Uz          = 0,
		    xmin        = self.xin - frozenWidth - widthInTight,
		    xmax        = self.xin
		)
		self.thermostatInSoft = Langevin(
		    beta        = betaInSoft, # ps^-1 
		    temperature = Tin, # K 
		    time_step   = dt,  # ps
		    heat_x      = 1,
		    heat_y      = 1,
		    heat_z      = 1,
		    Ux          = self.uin,
		    Uy          = 0,
		    Uz          = 0,
		    xmin        = self.xin - frozenWidth - widthInTight - widthInSoft ,
		    xmax        = self.xin - frozenWidth - widthInTight 
		)
		self.thermostatOut = Langevin(
		    beta        = betaOut, # ps^-1 
		    temperature = Tout, # K 
		    time_step   = dt,   # ps
		    heat_x      = 1,
		    heat_y      = 1,
		    heat_z      = 1,
		    Ux          = self.uout,
		    Uy          = 0,
		    Uz          = 0,
		    xmin        = self.xout,
		    xmax        = self.xout + widthOut 
		)

		thermostats.add(self.thermostatInTight)
		thermostats.add(self.thermostatInSoft)
		thermostats.add(self.thermostatOut)

		return thermostats

	def frozen_save(self, atoms):
		self.ifrozen = np.where(atoms.data["r"]["x"] > self.xin - self.frozenWidth)
		self.frozen = np.array(atoms.data[self.ifrozen])

	def frozen_step(self, atoms, dt):
		atoms.data["r"][self.ifrozen] = self.frozen["r"][:]
		atoms.data["v"][self.ifrozen] = self.frozen["v"][:]
		atoms.data["r"]["x"][self.ifrozen] += self.uin*dt

	def update(self, atoms, xfront, dt):
		# remove atoms
		iremove = np.where(atoms.data["r"]["x"] < self.xout)
		atoms.erase(iremove)
		# insert atoms
		Lpass = abs(self.uin * dt)
		iinsert = []
		xinsert = []
		if self.bpos + Lpass < self.period:
			iinsert = np.where(
				(self.buffer["r"]["x"] >= self.bpos)*
				(self.buffer["r"]["x"] < self.bpos + Lpass)
			)
			xinsert = self.buffer["r"]["x"][iinsert] - self.bpos + self.xin - Lpass
			iinsert = iinsert[0]
			self.bpos += Lpass
		else:
			iinsert1 = np.where(self.buffer["r"]["x"] >= self.bpos)
			iinsert2 = np.where(self.buffer["r"]["x"] < Lpass - (self.period - self.bpos))
			xinsert1 = self.buffer["r"]["x"][iinsert1] - self.bpos + self.xin - Lpass
			xinsert2 = self.buffer["r"]["x"][iinsert2] + self.xin + (self.period - self.bpos) - Lpass
			iinsert = np.concatenate((iinsert1[0],iinsert2[0]))
			xinsert = np.concatenate((xinsert1,xinsert2))
			self.bpos = Lpass - (self.period - self.bpos)

		oldsize = atoms.data.size
		atoms.resize(oldsize + iinsert.size)
		atoms.data[oldsize:] = self.buffer[iinsert]
		atoms.data["r"]["x"][oldsize:] = xinsert
		atoms.data["v"]["x"][oldsize:] += self.uin

		# shift velocities according to shock front propagation
		omega = (xfront - self.xout)/(self.xfront0 - self.xout)
		self.history["xfront"][self.hindex] = self.history["xfront"][self.hindex - 1] + \
		                                (xfront - self.xfront) + abs(self.uin)*dt
		self.history["time"][self.hindex] = self.history["time"][self.hindex - 1] + dt
		self.hindex += 1
		if self.hindex > 9:
			self.hready = True
			self.hindex = 0
		
		self.xfront = xfront

		vshock = abs(self.uin)
		if self.hready:
			x = self.history["xfront"] - np.amin(self.history["xfront"])
			t = self.history["time"] - np.amin(self.history["time"])
			result = least_squares(linear, [1.0, 0.0], args=(t, x))
			vshock = result.x[0]

		dv = 0.005
		if omega > 1.0:
			dv = abs(self.uin) - vshock
		else:
			dv = min(abs(self.uout), dv)

		self.uin += dv
		self.uout = min(self.uout + dv, 0.0)
		self.thermostatInTight.set_average_velocity(self.uin)
		self.thermostatInSoft.set_average_velocity(self.uin)
		self.thermostatOut.set_average_velocity(self.uout)
		self.vshock = vshock

		atoms.data["v"]["x"] += dv