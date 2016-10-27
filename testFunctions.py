#!/usr/bin/env python

import functions as F
import numpy as N
import unittest

class TestFunctions(unittest.TestCase):
    def testApproxJacobian1(self):
        slope = 3.0
        def f(x):
            return slope * x + 5.0
        x0 = 2.0
        dx = 1.e-3
        Df_x = F.ApproximateJacobian(f, x0, dx)
        self.assertEqual(Df_x.shape, (1,1))
        print Df_x
        print slope
        self.assertAlmostEqual(Df_x, slope)

    def testApproxJacobian2(self):
        A = N.matrix("1. 2.; 3. 4.")
        def f(x):
            return A * x
        x0 = N.matrix("5; 6")
        dx = 1.e-6
        Df_x = F.ApproximateJacobian(f, x0, dx)
        self.assertEqual(Df_x.shape, (2,2))
        N.testing.assert_array_almost_equal(Df_x, A)

    def testApproxJacobian3(self):
        A = N.matrix("1. 2. 3.; 4. 5. 6.; 7. 8. 9.")
        def f(x):
            return A * x
        x0 = N.matrix("10; 11; 12")
        dx = 1.e-6
        Df_x = F.ApproximateJacobian(f, x0, dx)
        self.assertEqual(Df_x.shape, (3,3))
        N.testing.assert_array_almost_equal(Df_x, A)

    def testPolynomial(self):
        # p(x) = x^2 + 2x + 3
        p = F.Polynomial([1, 2, 3])
        for x in N.linspace(-2,2,11):
            self.assertEqual(p(x), x**2 + 2*x + 3)

    def testAnalytic1(self):
        f = F.Polynomial([1,4,4])
        Df = F.Polynomial([0,2,4])
        x0=-3
        Df_x = F.AnalyticalJacobian(Df,x0)
        self.assertAlmostEqual(Df_x,-2.0)

    def testAnalytic2(self):
        f = lambda x: N.matrix([[6*x[0,0]+2*x[1,0]],[-9*x[0,0]-3*[1,0]]])
        Df = lambda x: N.matrix([[6+2*x[1,0]],[-3-9*x[0,0]]])
        x0 = N.matrix([[2],[-4]])
        Df_x = F.AnalyticalJacobian(Df,x0)
        N.testing.assert_array_almost_equal(Df_x,N.matrix([[-2.0],[-21.0]]))

    def testAnalytictoApprox(self):
        f = F.Polynomial([1,-7,12])
        dx = 1.e-6
        Df = F.Polynomial([0,2,-7])
        x0 = 5
        Df_xap = F.ApproximateJacobian(f,x0,dx)
        Df_xan = F.AnalyticalJacobian(Df,x0)
        self.assertTrue(abs(Df_xap-Df_xan)<0.000005) #Not close enough to pass assert almost equal   

if __name__ == '__main__':
    unittest.main()



