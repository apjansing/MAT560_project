from __future__ import print_function
from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
import pylab
import numpy as np
# from re import sub

class EigenSolver:
	"""docstring for EigenSolver"""
	def __init__(self, U, V, mesh):
		self.U = U
		self.V = V
		self.mesh = mesh

	def getEigenVectorValue(self, saveData = False, saveDataDir = None, saveDateName = "output"):
		# Test for PETSc and SLEPc
		if not has_linear_algebra_backend("PETSc"):
			print("DOLFIN has not been configured with PETSc. Exiting.")
			exit()

		if not has_slepc():
			print("DOLFIN has not been configured with SLEPc. Exiting.")
			exit()

		if(self.mesh == None):
			print('Mesh not present.')
			exit()

		if(self.U == None):
			self.U = FiniteElement('CG', triangle, 1)
		if(self.V == None):
			self.V = FiniteElement('CG', triangle, 1)
		parameters['reorder_dofs_serial'] = False
		W = FunctionSpace(self.mesh, self.V * self.U)

		class DirichletBoundary(SubDomain):
			def inside(self, x, on_boundary):
				return on_boundary

		# Define boundary condition
		u0 = Constant(0.0)
		boundaryCondition1 = DirichletBC(W.sub(0), u0, DirichletBoundary())
		boundaryCondition2 = DirichletBC(W.sub(1), u0, DirichletBoundary())

		# Define the bilinear form
		(u, v) = TrialFunction(W)
		(f1, f2) = TestFunction(W)
		a = -(dot(grad(u), grad(f2)) + dot(grad(v), grad(f1)))*dx
		L = (u*f1 + v*f2)*dx

		# Create the matrices
		A = PETScMatrix()
		b = PETScMatrix()
		assemble(a, tensor = A)
		assemble(L, tensor = b)
		boundaryCondition1.apply(A)
		boundaryCondition2.apply(A)

		# Create eigensolver
		eigensolver = SLEPcEigenSolver(A, b)

		# Compute all eigenvalues of A x = \lambda x
		print("Computing eigenvalues. This can take a minute.")
		eigensolver.solve()

		print("found: ", eigensolver.get_number_converged())

		eigenVectors = []
		eigenValues = []

		# If we're saving the data, make the directories to save the data to.
		if saveData:
			if saveDataDir == None:
				saveDataDir = 'eigen'
			mkdirs(saveDataDir)

		u = Function(W)
		for i in reversed(range(eigensolver.get_number_converged())):
			r, c, rx, cx = eigensolver.get_eigenpair(i)
			eigenVectors += [ rx.vec().getArray()[:len(rx)/2] ]
			eigenValues += [ r ]
			if saveData:
				u.vector()[:] = rx
				plot(u.sub(0))
				# save images
				pylab.savefig('%s/images/%s%d_%04d.png' % (saveDateName, saveDateName, m, i), 
					bbox_inches='tight')

		eigenVectors = np.array(eigenVectors)
		eigenValues = np.array(eigenValues)

		if saveData:
			np.save('%s/%sEigenVectors.npy' % (saveDateName, saveDateName), eigenVectors)
			np.save('%s/%sEigenValues.npy' % (saveDateName, saveDateName), eigenValues)
		return eigenVectors, eigenValues


	def mkdirs(s):
		try:
			mkdir('%s' % (s))
		except Exception as e:
			pass
		try:
			mkdir('%s/images' % (s))
		except Exception as e:
			pass
		try:
			mkdir('%s/vectors' % (s))
		except Exception as e:
			pass
