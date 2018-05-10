from polyplate import *
from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
from sys import argv
import numpy as np
import pylab

def main(xmlMeshFilename):
        parameters['reorder_dofs_serial'] = False
        mesh = Mesh(xmlMeshFilename)
        U = FiniteElement('CG', triangle, 1)
        V = FiniteElement('CG', triangle, 1)
        W = FunctionSpace(mesh, V)

        solver = EigenSolver(U, V, mesh)
        eigenVectors, eigenValues = solver.getEigenVectorValue()

        p = Expression('x[0] - x[1] * x[0]', degree=2)
        v = Expression('sin((0.25 - pow(x[0], 2) - pow(x[1], 2)) * 30)', degree=2) 
        simulator = Simulator(mesh, U, eigenVectors, eigenValues, p, v, 300)
        max_val = simulator.upperBound() * 0.8
        for t in range(200):
                print("Frame %3d / %d" % (t + 1, 200))
                plot(simulator.evaluate(t), range_min=-max_val, range_max=max_val)
                pylab.savefig('../solution/frame%04i.png' % (t + 1), bbox_inches='tight')



if __name__ == '__main__':
        if len(argv) < 2:
                main("gmsh/demo/meshes/plate/plate040.xml")
        else:
                main(argv[1])

