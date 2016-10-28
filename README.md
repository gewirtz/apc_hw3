# apc_hw3

These files contain a solver for finding the roots of a function f(x) = 0 using
Newton's method (a Newton-Raphson solver).

functions.py: Contains two functions-- a) ApproximateJacobian, which returns an
	approximation of the Jacobian, and b) AnalyticalJacobian, which returns 
	an analytical solution of the Jacobian. It also contains the Polynomial
	class, which we can use to construct a polynomial object.

newton.py: Contains the Newton class, which is an object useful for finding
	roots of f(x) = 0 from Newton's method. Arguments for this class are:
	-f: the function
	-tol: the tolerance, basically iterate until |f(x)| <  tol
	-maxiter: maximum number of iterations to do
	-dx: step size
	-Df: the analytic Jacobian, absence of which leads to use of the approximate Jacobian
	-r: The approximate solution has to lie within r of the initial guess x0, or else exception
	This class has two member functions: a)solve, which returns the root of f(x) = 0 after being given
	an initial guess x0, and b) step, which iterates through a single step
	of Newton's method.

testNewton.py: Contains 9 tests to make sure the Newton class works correctly:
	-testLinear: tests solver for linear polynomial
	-testLinear2var: tests solver for linear polynomial with 2 variables
	-testLinear 3var: tests solver for linear polynomial w/ 3 vars
	-testQuadratic: tests solver for quadratic polynomial
	-testPeriodic: tests solver for a periodic function
	-testStep: tests that a) solver is stepping in correct direction, b) if Df is provided that the analytic Jacobian function is getting called, and c) that the approximate solutions found by the analytic and approximate Jacobian are almost equal
	-testSolve: tests the same 3 things for testStep, just calls the solve method vs step.
	-testNotConverge: tests that an exception is raised when the method doesn't converge before maxiter
	-testExceedRadius: tests that an exception is raised when r is exceeded.

testFunctions.py: Contains 7 tests to make sure the Jacobians are being computed
	correctly and that the Polynomial class works as expected:
	-testApproxJacobian1: tests that the approximateJacobian function returns a matrix of the correct dimensions and the correct value for a 1D function
	-testApproxJacobian2: same tests as above for a 2D function
	-testApproxJacobian3: same tests as above for a function with 3 variables.
	-testPolynomial: tests that the polynomial class works as expected
	-testAnalytic1: tests that the AnalyticalJacobian function returns expected output for one variable
	-testAnalytic2: same tests as testAnalytic1 but for 2 variables
	-testAnalytictoApprox: For single variable polynomials, tests that the results from calling the analyticalJacobian and approximateJacobian functions directly are about equal.
