#!/usr/bin/env python

import newton
import unittest
import numpy as N
import functions as F
import math

class TestNewton(unittest.TestCase):
    def testLinear(self):
        f = lambda x : 3.0 * x + 6.0
        solver = newton.Newton(f, tol=1.e-15, maxiter=10, r = 5)
        x = solver.solve(2.0)
        self.assertEqual(x, -2.0)

    def testLinear2var(self):
    	f = lambda x: N.matrix([[x[0,0]+x[1,0]+2],[-x[0,0]+x[1,0]+4]])
    	solver = newton.Newton(f, tol=1.e-15, maxiter=20)
    	x0 = N.matrix([[0],[-3.5]])
    	x = solver.solve(x0)
    	N.testing.assert_array_almost_equal(x,N.matrix([[1.0],[-3.0]]))

    def testLinear3var(self):
    	f= lambda x: N.matrix([[x[0,0]+x[1,0]+5*x[2,0]],[-x[0,0]-x[1,0]+3*x[2,0]+8],[2*x[0,0]-x[1,0]-x[2,0]-2]])
    	solver = newton.Newton(f, tol=1.e-15, maxiter=20)
    	x0 = N.matrix([[1],[4],[-1.5]])
    	x = solver.solve(x0)
    	N.testing.assert_array_almost_equal(x,N.matrix([[2.0],[3.0],[-1.0]]))

    def testQuadratic(self):
    	#tests x^2 -4x -12
    	f = F.Polynomial([1,-4,-12])
    	solver = newton.Newton(f, tol=1.e-15, maxiter=200)
    	x = solver.solve(-2.3)
    	self.assertAlmostEqual(x, -2.0)

    def testPeriodic(self):
    	f = lambda x: N.matrix([[math.sin(x[0,0])],[math.tan(x[1,0])]])
    	solver = newton.Newton(f, tol=1.e-15,maxiter=20)
    	x0 = N.matrix([[1],[1]])
    	x = solver.solve(x0)
    	N.testing.assert_array_almost_equal(x, N.matrix([[0],[0]]))

    def testStep(self):
    	f = F.Polynomial([1,4,4])
    	apsolver = newton.Newton(f, tol=1.e-15, maxiter=1)
    	x0 = -4
    	x=apsolver.step(x0)
    	self.assertTrue(x>x0)
    	Df = lambda x: 2*x + 4
    	ansolver = newton.Newton(f, tol=1.e-15, maxiter=1, Df=Df)
    	xan = ansolver.step(x0)
    	self.assertTrue(x != xan) #So we know a different step function is getting called
    	self.assertTrue(abs(x-xan)<0.0000005) #test that the difference is very small btw the two results


    def testSolve(self):
    	f = F.Polynomial([1,4,4])
    	apsolver = newton.Newton(f, tol=1.e-15, maxiter=100)
    	x0 = -3
    	x=apsolver.solve(x0)
    	self.assertTrue(x>x0)
    	Df = lambda x: 2*x + 4
    	ansolver = newton.Newton(f, tol=1.e-15, maxiter=100, Df=Df)
    	xan = ansolver.solve(x0)
    	self.assertTrue(x != xan) #So we know a different step function is getting called
    	self.assertTrue(abs(x-xan)<0.0000005) #test that the difference is very small btw the two results
    	self.assertAlmostEqual(-2.0,x)
    	self.assertAlmostEqual(-2.0,xan)

    def testNotConverge(self):
    	f = F.Polynomial([1,4,4])
    	solver = newton.Newton(f, tol=1.e-15,maxiter=2)
    	x0 = 1.8
    	self.assertRaises(Exception, solver.solve, x0)

    def testExceedRadius(self):
    	f = F.Polynomial([1,4,4])
    	x0=-6
    	apsolver=newton.Newton(f, tol=1.e-15,maxiter=1000)
    	self.assertRaises(Exception, apsolver.solve, x0)
    	df = lambda x: 2*x + 4
    	ansolver=newton.Newton(f, tol=1.e-15,maxiter=1000, Df=df)
    	self.assertRaises(Exception, ansolver.solve, x0)
    






if __name__ == "__main__":
    unittest.main()
