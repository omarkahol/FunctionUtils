package com.okahol.FunctionUtils;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.Future;

public class Calculus {

	private double tol; // Tolerance on approximations
	private int itmax; // Maximum number of iterations
	private static int count; // Static counter--do not touch--it is used to manage recursive algorithm
	private double alfa; // "gradient descent rate"
	private double h;// "the step value--used for derivatives

//Default constructor --> Recommended values
	public Calculus() {
		this.tol = 1E-15;
		this.itmax = 5000;
		this.alfa = 0.1;
		this.h = 1E-5;
	}

//constructor
	public Calculus(double tol, int itmax, double alfa, double h) {
		this.tol = tol;
		this.itmax = itmax;
		this.alfa = alfa;
		this.h = h;
	}

	/*
	 * RECURSIVE ALGORITHM TO FIND THE ZERO OF A FUNCTION USING BISECTION METHOD (NOT VERY PRECISE)
	 * [a,b]--> the interval in which the zero will be searched FunctionTool fun-->
	 * A lambda expression or an anonymous class representing a mathematical
	 * function PrintIter--> display the number of iterations OUTPUTS The zero of
	 * the function (if found)
	 * 
	 */
	public double fzerobisect(double a, double b, FunctionTool fun, boolean PrintIter) {
		if (fun.f(a) * fun.f(b) < 0) {
			count++;
			double xm = (a + b) / 2;
			if (fun.f(xm) == 0) {
				System.out.println("Exact solution found!!");
				return xm;
			}
			if (count - 1 < itmax && tol < Math.abs(fun.f(xm))) {
				if (fun.f(a) * fun.f(xm) < 0) {
					return fzerobisect(a, xm, fun, PrintIter);
				} else {
					return fzerobisect(xm, b, fun, PrintIter);
				}
			} else {
				if (count - 1 == itmax) {
					System.out.println("Itmax reached, result might be inaccurate");
					return xm;
				} else {
					int iter = count;
					if (PrintIter) {
						System.out.println("Solution found with " + iter + " iterations");
					}
					count = 0;
					return xm;
				}
			}
		} else {
			System.out
					.println("The function doesn't seem to have a zero as it doesn't verify the condition f(a)*f(b)<0");
			return 0;
		}
	}

	/*
	 * FIRST DERIVATIVE OF A FUNCTION INPUTS x0--> a point fun--> the function
	 * OUTPUTS the first derivative of the function in x0
	 * 
	 */
	public double df(double x0, FunctionTool fun) {
		double xp = (fun.f(x0 + h) - fun.f(x0 - h)) / (2 * h);
		return xp;
	}

	/*
	 * SECOND DERIVATIVE OF A FUNCTION -->Similar to the first derivative. It return
	 * the second derivative
	 * 
	 */
	public double df2(double x0, FunctionTool fun) {
		double xp2 = (fun.f(x0 + h) + fun.f(x0 - h) - 2 * fun.f(x0)) / Math.pow(h, 2);
		return xp2;
	}

	/*
	 * RECURSIVE IMPLEMENTATION OF THE 'GRADIENT DESCENT' ALGORITHM FOR FINDING
	 * MINIMA INPUTS x0--> the starting point, where the descent will start fun-->
	 * the function PrintIter--> print final number of iterations on screen OUTPUTS
	 * the minima
	 * 
	 */
	public double minf(double x0, FunctionTool fun, boolean PrintIter) {
		count++;
		double xp = df(x0, fun);
		if (count < itmax && tol < Math.abs(xp)) {
			return minf(x0 - alfa * xp, fun, PrintIter);
		} else {
			int iter = count;
			if (count == itmax) {
				System.out.println("Minima not found, Itmax reached");
			}
			if (PrintIter) {
				System.out.println("Solution found with " + iter + " iterations");
			}
			count = 0;
			return x0;
		}
	}

	/*
	 * FIND MAXIMA OF A FUNCTION --> See minf for more info
	 */
	public double maxf(double x0, FunctionTool fun, boolean PrintIter) {
		count++;
		double xp = df(x0, fun);
		if (count < this.itmax && tol < Math.abs(xp)) {
			return maxf(x0 + alfa * xp, fun, PrintIter);
		} else {
			int iter = count;
			if (count == itmax) {
				System.out.println("Maximum might be false, Itmax reached.");
			}
			if (PrintIter) {
				System.out.println("Solution found with " + iter + " iterations.");
			}
			count = 0;
			return x0;
		}
	}

	/*INTEGRAL OF A FUNCTION USING SIMPSON'S METHOD (PRECISE)
	 * INPUTS 
	 * 		[a,b]--> the interval of integration
	 * 		n --> number of subdivisions of the main interval [a,b]. Suggestion (5000 should work)
	 * 		fun--> the function
	 * OUTPUTS
	 * 		the integral value, if it exists
	 */
	public double quads(double a, double b, int n, FunctionTool fun) {

		double space = (b - a) / n;
		Vector x = linspace(a, b, space / 2);
		x = vfeval(x, fun);
		Vector weights = new Vector(2 * n + 1, "1");

		for (int k = 1; k <= 2 * n - 1; k++) {
			if (k % 2 == 0) {
				weights.setElement(k, 2);
			} else {
				weights.setElement(k, 4);

			}
		}
		double Integral = dot(x, weights);
		Integral = (Integral * space) / 6;
		return Integral;
	}
	/*
	 * INTEGRAL USING TRAPEZOIDAL RULE (NOT PRECISE)
	 * See quads for more info
	 */
	public double quadt(double a, double b, int n, FunctionTool fun) {
		double space = (b - a) / n;
		Vector fx = linspace(a, b, h);
		fx = vfeval(fx, fun);
		Vector weights = new Vector(n + 1, "1");
		weights.setElement(0, 0.5);
		weights.setElement(n, 0.5);
		weights.lambdaprod(h);
		double Integral = dot(fx, weights);
		return Integral;

	}

	/*RECURSIVE ALGORITHM TO FIND THE ZERO OF A FUNCTION USING NEWTON'S RULE (PRECISE, FAST AND STABLE)
	 * INPUTS
	 * 		x0-->starting point. It should be as close as possible to the real zero
	 * 		fun--> function
	 * 		PrintIter--> prints iterations
	 * 
	 */
	public double fzeronewton(double x0, FunctionTool fun, boolean PrintIter) {
		count++;
		double x1 = x0 - ((fun.f(x0)) / (df(x0, fun)));
		double scarto = Math.abs(x1 - x0);
		if (count - 1 < itmax && scarto > tol) {
			return fzeronewton(x1, fun, PrintIter);
		} else {
			if (count - 1 == itmax) {
				System.out.println("Itmax reached results might be inaccurate");
			}
			int iter = count - 1;
			if (PrintIter) {
				System.out.println("Solution found with " + iter + " iterations");
			}
			count = 0;
			return x1;
		}
	}
	/*
	 * ANOTHER RECURSIVE ALGORITHM FOR FINDING THE ZERO
	 * The inputs are similar to newton's. It is faster but less precise and unstable
	 */
	public double fZeroFast(double x0, FunctionTool fun, boolean PrintIter) {
		count++;
		double fx = fun.f(x0);
		double df = df(x0, fun);
		double df2 = df2(x0, fun);
		double sqrtarg = Math.abs(Math.pow(df, 2) - 2 * fx * df2);
		double x1 = x0 - ((df - Math.signum(Math.pow(df, 2) - 2 * fx * df2) * Math.sqrt(sqrtarg)) / (df2));
		double scarto = Math.abs(x1 - x0);
		if (count - 1 < itmax && scarto > tol) {
			return fZeroFast(x1, fun, PrintIter);
		} else {
			if (count - 1 == itmax) {
				System.out.println("Itmax reached results might be inaccurate");
			}
			int iter = count - 1;
			if (PrintIter) {
				System.out.println("Solution found with " + iter + " iterations");
			}
			count = 0;
			return x1;
		}
	}
	/*
	 * ANOTHER ALGORITHM FOR FINDING THE ZERO
	 * This one is slow but efficient and precise. The ck coefficient determines how fast it should be. Recommended value is one as
	 * an increase in velocity causes a decrease in stability and precision --> Anyways it should not be smaller or equal to zero
	 */
	public double fzero(double x0, double ck, FunctionTool fun, boolean PrintIter) {
		if (ck <= 0)
			throw new IllegalValueException("coefficient must be greater than zero");
		count++;
		double fx = fun.f(x0);
		double x1 = x0 + fx / ck;
		double scarto = Math.abs(x1 - x0);
		if (count - 1 < itmax && scarto > tol) {
			return fzero(x1, ck, fun, PrintIter);
		} else {
			if (count - 1 == itmax) {
				System.out.println("Itmax reached results might be inaccurate");
			}
			int iter = count - 1;
			if (PrintIter) {
				System.out.println("Solution found with " + iter + " iterations");
			}
			count = 0;
			return x1;
		}
	}

	/*INTEGRAL APPROXIMATION USING RICHARDSON RULE (VERY VERY VERY ... PRECISE)
	 * 
	 * The inputs and the output are similar to the other quad methods.
	 * Enhanced precision (if set to true) makes it even better
	 * Recommended value for n is 2000
	 */
	public double quadr(double a, double b, int n, boolean EnhancePrecision, FunctionTool fun) {
		if (!(n % 2 == 0)) {
			n = n + 1;
		}
		if (EnhancePrecision) {
			double Q1 = quads(a, b, (n / 2), fun);
			double Q2 = quads(a, b, n, fun);
			double QR = (16 * Q2 - Q1) / 15;
			return QR;
		} else {
			double Q1 = quadt(a, b, (n / 2), fun);
			double Q2 = quadt(a, b, n, fun);
			double QR = (4 * Q2 - Q1) / 3;
			return QR;
		}
	}
	/*
	 * Calculates Vector and Matrix Product
	 */
	public Vector vmprod(Matrix m, Vector v) {
		Vector vmp = new Vector(m.getRow());
		if (m.getCol() == v.length()) {
			for (int i = 0; i < m.getRow(); i++) {
				vmp.setElement(i, dot(v, m.getArow(i)));
			}
			return vmp;
		} else {
			return null;
		}
	}
	/*
	 * Gaussian elimination and LU decomposition. Returns ONLY the U matrix
	 */
	public Matrix LuGauss(Matrix m) {
		int n = m.getRow();
		Matrix L = new Matrix(n, "eye");
		Matrix U = m.clone();
		for (int k = 0; k < n - 1; k++) {
			for (int i = k + 1; i < n; i++) {
				L.setElement(U.getElements()[i][k] / U.getElements()[k][k], i, k);
				for (int j = k; j < n; j++) {
					U.setElement(U.getElements()[i][j] - (L.getElements()[i][k] * U.getElements()[k][j]), i, j);
				}

			}
		}
		return U;
	}
	/*
	 * Calculates the determinant of a Matrix
	 */
	public double det(Matrix m) {
		Matrix U = LuGauss(m);
		double prod = 1;
		for (int i = 0; i < U.getRow(); i++) {
			prod = prod * U.getElements()[i][i];
		}
		return prod;
	}
	/*
	 * maps a function to each element of a vector
	 */
	public Vector vfeval(Vector v, FunctionTool fun) {
		Vector vclone = v.clone();
		for (int i = 0; i < v.length(); i++) {
			vclone.setElement(i, fun.f(v.getElements()[i]));
		}
		return vclone;
	}
	/*
	 * Calculates the trace of a matrix
	 */
	public double trace(Matrix m) {
		double sum = 0;
		for (int i = 0; i < m.getRow(); i++) {
			sum += m.getElements()[i][i];
		}
		return sum;
	}
	/*
	 * Returns an equally-spaced vector
	 */
	public Vector linspace(double start, double end, double step) {
		double size = Math.floor((end - start) / step) + 1;
		Vector v = new Vector((int) size);
		for (int i = 0; i < (int) size; i++) {
			v.setElement(i, start + i * step);
		}
		return v;
	}

	
	/*
	 * maps the first derivative of a function to a vector
	 */
	public Vector vdfeval(Vector v, FunctionTool fun) {
		for (int i = 0; i < v.length(); i++) {
			v.setElement(i, df(v.getElements()[i], fun));
		}
		return v;
	}
	
	/*
	 * maps the second derivative of a function to a vector
	 */
	public Vector vdf2eval(Vector v, FunctionTool fun) {
		for (int i = 0; i < v.length(); i++) {
			v.setElement(i, df2(v.getElements()[i], fun));
		}
		return v;
	}
	
	/*
	 * Approximates the derivative of a function known only in a set of points
	 * The x vector represents the x coordinate of the points
	 * The y vector represents the y coordinate of the points
	 */
	public Vector vdfapprox(Vector x, Vector y) {
		Vector df = new Vector(x.length() - 2);
		for (int i = 1; i < x.length() - 1; i++) {
			df.setElement(i, y.getElements()[i + 1]
					- y.getElements()[i - 1] / (x.getElements()[i + 1] - x.getElements()[i - 1]));
		}
		return df;
	}
	
	/*
	 * Approximates the second derivative of a function known only in a set of points
	 * The x vector represents the x coordinate of the points
	 * The y vector represents the y coordinate of the points
	 */
	public Vector vdf2approx(Vector x, Vector y) {
		Vector df2 = new Vector(x.length() - 2);
		for (int i = 1; i < x.length() - 1; i++) {
			df2.setElement(i, (y.getElements()[i + 1] + y.getElements()[i - 1] - 2 * y.getElements()[i])
					/ (Math.pow(x.getElements()[i] - x.getElements()[i - 1], 2)));
		}
		return df2;
	}
	/*
	 * Returns the sum of all the elements of a vector
	 */
	public double tracev(Vector v) {
		double trace = 0;
		for (double el : v.getElements()) {
			trace += el;
		}
		return trace;
	}
	
	/*
	 * Element-wise product of two vectors
	 */
	public Vector vprodv(Vector vector1, Vector vector2) {
		Vector p = new Vector(vector1.length());
		for (int i = 0; i < vector1.length(); i++) {
			p.setElement(i, vector1.getElements()[i] * vector2.getElements()[i]);
		}
		return p;
	}
	/*
	 * dot product between two Vectors
	 */
	public double dot(Vector v1, Vector v2) {
		if (v1.length() == v2.length()) {
			double prod = 0;
			for (int i = 0; i < v1.length(); i++) {
				prod = prod + (v1.getElements()[i] * v2.getElements()[i]);
			}
			return prod;
		} else {
			System.out.println("Vector Dimensions must agree");
			return 0;
		}
	}
	/*
	 * Product between two matrixes 
	 */
	public Matrix mprod(Matrix m1, Matrix m2) {
		if (m1.getCol() == m2.getRow()) {
			Matrix prod = new Matrix(m1.getRow(), m2.getCol(), "0");
			for (int i = 0; i < prod.getCol(); i++) {
				for (int k = 0; k < m1.getRow(); k++) {
					prod.setElement(dot(m1.getArow(k), m2.getAcol(i)), k, i);
				}
			}
			return prod;
		} else {
			return null;
		}
	}
	
	
	/*
	 * multiplies lambda to every element of a vector
	 */
	public Vector vlambda(Vector v, double lambda) {
		Vector v1 = new Vector(v.length());
		v1 = sumv(v1, v);
		v1.lambdaprod(lambda);
		return v1;
	}
	
	/*
	 * multiplies lambda to every element of a matrix
	 */
	public Matrix mlambda(Matrix m, double lambda) {
		Matrix m1 = new Matrix(m.getRow(), m.getCol());
		m1 = msum(m1, m);
		m1.lambda(lambda);
		return m1;
	}
	/*
	 * sums two vectors
	 */
	public Vector sumv(Vector v1, Vector v2) {
		Vector vsum = new Vector(v1.length(), "0");
		vsum.sum(v1);
		vsum.sum(v2);
		return vsum;
	}

	/*
	 * sums two matrixes
	 */
	public Matrix msum(Matrix m1, Matrix m2) {
		Matrix sum = new Matrix(m1.getRow(), m1.getCol());
		sum.sum(m1);
		sum.sum(m2);
		return sum;
	}
	/*
	 * Works out the inverse of a matrix. At least, it tries to...
	 * works better if LuGaussFirst is set to True
	 */
	public Matrix inv(Matrix m, boolean LuGaussFirst) {
		if (m.getCol() == m.getRow() && Math.abs(det(m)) > 0) {
			Matrix theinv = new Matrix(m.getCol(), "0");
			Matrix work = new Matrix(m.getCol(), "0");
			work = msum(work, m);

			if (LuGaussFirst) {
				work = LuGauss(work);
			}

			for (int i = 0; i < m.getCol(); i++) {
				theinv.setAcol(solve(work, (new Matrix(work.getCol(), "eye")).getAcol(i),
						new Vector(work.getCol(), "0"), false), i);
			}
			return theinv;
		}
		return null;
	}

	/*
	 * maps the power function to each element
	 */
	public Vector vpow(Vector v, double pow) {
		Vector v1 = v.clone();
		v1.pow(pow);
		return v1;
	}

	/*
	 * Solves  linear system of equations using Gauss-Seidel algorithm
	 * IPUTS
	 * 		A-->Matrix 
	 * 		b-->Vector
	 * 		x0-->Starting vector ( just use new Vector(n,"n:0") )
	 * OUPUS
	 * 		a vector, the solution		
	 */
	public Vector solve(Matrix A, Vector b, Vector x0, boolean PrintIter) {
		int iter = 0;
		double diff = this.tol + 1;
		while (iter < this.itmax && diff > tol) {

			for (int i = 0; i < x0.length(); i++) {
				double sum = 0;
				for (int j = 0; j < x0.length(); j++) {
					if (!(j == i)) {
						sum -= A.getElements()[i][j] * x0.getElements()[j];
					}
				}
				x0.setElement(i, (Math.pow(A.getElements()[i][i], -1)) * (b.getElements()[i] + sum));

			}
			diff = norm((sumv(b, vlambda(vmprod(A, x0), -1))));
			iter += 1;
		}
		if (PrintIter) {
			System.out.println("Solution found with " + iter + " iterations");
		}
		return x0;
	}

	/*
	 * SIMILAR TO SOLVE BUT BETTER
	 * it implements the successive-over-relaxation (SOR) algorithm 
	 * W should be between 0 and 2
	 * Better use solve if you do not know the optimal value for w
	 */
	public Vector solvew(Matrix A, Vector b, Vector x0, double w, boolean PrintIter) {
		int iter = 0;
		Vector x1 = new Vector(x0.length());
		x1.setElements(x0.getElements());
		double diff = this.tol + 1;

		while ((iter < this.itmax) && (diff > this.tol)) {
			for (int i = 0; i < x0.length(); i++) {
				double somma = 0;
				for (int j = 0; j < x0.length(); j++) {
					if (j < i) {
						somma -= A.getElements()[i][j] * x1.getElements()[j];
					}
					if (j > i) {
						somma -= A.getElements()[i][j] * x0.getElements()[j];
					}
				}
				x1.setElement(i, w * (Math.pow(A.getElements()[i][i], -1)) * (b.getElements()[i] + somma)
						+ (1 - w) * x0.getElements()[i]);
			}
			diff = norm((sumv(b, vlambda(vmprod(A, x0), -1))));
			iter += 1;
			x0 = x1;

		}
		if (PrintIter) {
			System.out.println("Solution found with " + iter + " iterations");
		}
		return x1;
	}
	
	
	//THOSE METHODS ARE NOT FULLY TESTED YET
	
	
	public double norm(Vector v) {
		return dot(v, v);
	}

	public double limit(FunctionTool fun, double value) {
		double limit = fun.f(value - this.h);
		if (Math.abs(Math.floor(limit) - limit) < this.tol)
			return Math.floor(limit);
		else if (Math.abs(Math.ceil(limit) - limit) < this.tol)
			return Math.ceil(limit);
		else
			return limit;
	}

	public Matrix multiThreadInverse(Matrix m) throws InterruptedException, ExecutionException {
		Matrix inv = new Matrix(m.getCol(), "0");
		Matrix eye = new Matrix(m.getCol(), "eye");
		Future<Vector>[] storage = new Future[m.getCol()];
		ExecutorService ex = Executors.newCachedThreadPool();
		for (int i = 0; i < m.getCol(); i++) {
			storage[i] = ex.submit(new SolveThread(m, eye.getAcol(i), this));
		}
		ex.shutdown();
		ex.awaitTermination((long) this.itmax, TimeUnit.SECONDS);
		int counter = 0;
		for (Future<Vector> col : storage) {
			inv.setAcol(col.get(), counter);
			counter++;
		}
		return inv;
	}


}
