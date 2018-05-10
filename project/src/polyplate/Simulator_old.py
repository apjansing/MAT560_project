from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
from math import sin, cos
from numpy import array, dot


class Simulator_old():
	"""docstring for Simulator"""
	def __init__(self, mesh, U, eigenVectors, eigenValues,
	 pExpression = Expression('sin(x[0])', degree=1), 
	 vExpression = Expression('cos(x[0])', degree=1)):
		# Test for PETSc and SLEPc
		if not has_linear_algebra_backend("PETSc"):
			print("DOLFIN has not been configured with PETSc. Exiting.")
			exit()

		if not has_slepc():
			print("DOLFIN has not been configured with SLEPc. Exiting.")
			exit()
		self.mesh = mesh
		self.U = U
		self.fSpace = FunctionSpace(self.mesh, self.U)
		self.eigenVectors = eigenVectors
		self.eigenValues = eigenValues
		self.pExpression = pExpression
		self.vExpression = vExpression
		self.p = None
		self.v = None


	def interp(self, expression):
		return fem.interpolation.interpolate(expression, self.fSpace)

	def interpolatePositionPortion(self):
		self.p = self.interp(self.pExpression)
		return self.p
	
	def interpolateVelocityPortion(self):
		self.v = self.interp(self.vExpression)
		return self.v

	def vecAtTime(self, t):
		if self.v == None:
			self.interpolateVelocityPortion()
		if self.p == None:
			self.interpolatePositionPortion()
		V = list(array(self.v.vector()))
		V += V
		P = list(array(self.p.vector()))
		P += P
		S = 0
		for i in range(len(V)):
			S += cos(self.eigenValues[i]*t)*dot(dot(V, self.eigenVectors[i])/dot(self.eigenVectors[i],self.eigenVectors[i]), self.eigenVectors[i]) + sin(self.eigenValues[i]*t)*dot(dot(P, self.eigenVectors[i])/dot(self.eigenVectors[i],self.eigenVectors[i]), self.eigenVectors[i])
		return S

	def getOmegaAtT(self, t = 0):
		# omega = sum[cos(wt)(ep/ee)e] + sum[sin(wt)(ev/ee)e]
		vec = self.vecAtTime(t)# + self.getPositionPortion(t)
		f = Function(self.fSpace)
		print(f.vector()[:])
		print(len(vec))
		f.vector()[:] = vec
		return f