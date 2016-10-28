# apc_hw3

These files contain a solver for finding the roots of a function f(x) = 0 using
Newton's method.

functions.py: Contains two functions-- a) ApproximateJacobian, which returns an
	approximation of the Jacobian, and b) AnalyticalJacobian, which returns 
	an analytical solution of the Jacobian. It also contains the Polynomial
	class, which we can use to construct a polynomial object.

newton.py: Contains the Newton class, which is an object useful for finding
	roots of f(x) = 0 from Newton's method. This class has two member 
	functions: a)solve, which returns the root of f(x) = 0 after being given
	an initial guess x0, and b) step, which iterates through a single step
	of Newton's method.

testNewton.py: Contains 9 tests to make sure the Newton class works correctly.

testFunctions.py: Contains 7 tests to make sure the Jacobians are being computed
	correctly and that the Polynomial class works as expected.
