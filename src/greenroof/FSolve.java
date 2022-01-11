package greenroof;

import java.util.ArrayList;
import java.util.List;

//import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
//import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
//import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver;
//import org.apache.commons.math3.exception.DimensionMismatchException;

import com.mathlibrary.equation.NonLinearSys;
import com.mathlibrary.function.FunctionXs;
import com.mathlibrary.polynomial.Polynomial;

public class FSolve
{
	
	public static final double STUCK_VALUE = 1E-15;
	public static final double ERROR_RETURN = -9999.;
	public static final int ITERATIONS = 10000;
	double a = -1; /* left point */
	double b = 50; /* right point */
	boolean DEBUG=true;
	boolean SKIP=true;
	double TOL = 1.0E-6; /* tolerance */

	public static void main(String[] args)
	{
		FSolve f = new FSolve();
//		f.test5();
//		f.bisectionTesting();
		f.testSimple();

	}
	
	
//	public double[] balance_points(double x_min, double x_max, double y_min, double y_max, double delta, double eps)
//	{
//		eps=2E-32;
//
//	    double width = x_max - x_min;
//	    double height = y_max - y_min;
//	    double x_middle = (x_min + x_max)/2.;
//	    double y_middle = (y_min + y_max)/2.;
//
//	    double Fx = Math.min(f(C)*f(B), f(D)*f(A));
//	    double Fy = Math.min(f(A)*f(B), f(D)*f(C));
//	    double Gx = Math.min(g(C)*g(B), g(D)*g(A));
//	    double Gy = Math.min(g(A)*g(B), g(D)*g(C));
//
//	    double F = Math.min(Fx, Fy);
//	    double G = Math.min(Gx, Gy);
//
//	    double largest_dim = Math.max(width, height);
//	    double guaranteed_contain_zeros = Math.max(F, G) < 0 or (f(C) == 0 and g(C) == 0);
//
//	    if (guaranteed_contain_zeros and largest_dim <= eps)
//	    {
//	        return new double[]{x_middle, y_middle};
//	    }
//	    else if (guaranteed_contain_zeros or largest_dim > delta)
//	    {
//	        if (width >= height:)
//	        {
//	            return balance_points(x_min, x_middle, y_min, y_max, delta) + balance_points(x_middle, x_max, y_min, y_max, delta);
//	        }
//	        else
//	        {
//	            return balance_points(x_min, x_max, y_min, y_middle, delta) + balance_points(x_min, x_max, y_middle, y_max, delta);
//	        }
//	    }
//	    else
//	    {
//	        return new double[0];
//	    }
//	}
	
	public void testSimple()
	{
//		double a = 4;
//		double b = 5;
		
		double D_theta = -0.22;
		double psi = 1.11;
		double K = 1.5;
		double dt = 1.4;
		double F0 = 2.1142;
		
		double a = 0.25;
		double b = 10000 ;
		
		double value = bisectionalSolver(a, b, D_theta, psi, K, dt, F0);
		System.out.println(value);
	}
	
//	public void testGreenAmpt()
//	{
//		
//		
//		F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
//		
//	}
	
	
//	public void bisectionTesting()
//	{
//		double dz=15.0 ;
//		double ref_ta =23.5 ;
//		double UTb =  0.21266443039771887 ;
////		double[] mod_U_TaRef = new double[1];
//		double mod_U_TaRef =  0.1  ;
//		double Ri_rur = -370.55369919381894;
//		//TODO this should be 22.760359533246
//		
//		double c = bisectionalSolver(dz, ref_ta, UTb, mod_U_TaRef, Ri_rur);
//	}
	
	// adapted from https://x-engineer.org/bisection-method/
	public double bisectionalSolver2(double a, double b)
	{
		double localA = a;
		double localB = b;
		
		double NMAX = ITERATIONS; /* maximum number of iterations */
		double c = 0; /* estimated root */
		int index = 0; /* index */
		int stuckCount = 0;
		 
		c = (localA + localB)/2.0; /* calculate the midpoint */
		 /* Evaluate loop until the result is less than the tolerance
		 * maximum number of iterations is not yet reached*/
		 
		 if (calculateExpression(c) == 0)
		 {
		  /* If the first midpoint gives f(c) = 0, c is the root */
	//		 System.out.println("root is " + c);
		 }
		 else
		 {
			 while ((Math.abs(calculateExpression(c)) > TOL) && (index<=NMAX))
			 {
				 if (sign(calculateExpression(c)) == sign(calculateExpression(localA)))
				 {
					 /* f(c) has same sign as f(a) */
					 localA = c;
				 }
				 else
				 {
					 /* f(c) has same sign as f(b) */
					 localB = c;
				 } 
				 c = (localA+localB)/2.0; /* midpoint update */
				 double difference = Math.abs(localA-localB);
				 if (DEBUG)
				 {
					 System.out.println( (index+1) + " " + localA + " " + localB + " " + difference + " " + c );
				 }
				
				 //did it get stuck? Try different variations of end points
				 if (difference < STUCK_VALUE  )
				 {					 
					 if (stuckCount > 1)
					 {
						 localA = -1. - Math.random();
						 localB = 1. + Math.random();
					 }					 
					 if (stuckCount > 2)
					 {
						 localA = 0.0-Math.random();
						 localB =  0.0;
					 }
					 if (stuckCount > 3)
					 {				 
						 if (stuckCount % 2 == 0)
						 {
							 localA = -0.001;
							 localB =  0.0-Math.random()/1000.;
						 }
						 else
						 {
							 localA = -0.001;
							 localB =  Math.random();
						 }
					 }
					 if (stuckCount < 1)
					 {
						 c = c - 400;
					 } 
					 stuckCount ++;
				 }
	
				 index++; /* index increment */
			 }
		 } 
		 
		 if (index>=NMAX)
		 {
			 if (DEBUG)
			 {
				 System.out.println("Root not found " + c + " after " + index + " iterations");	
			 }
			 
			 return ERROR_RETURN;
		 }
		 
		 /* Display results */
		 if (DEBUG)
		 {
			 System.out.println("Root is " + c + " found after " + index + " iterations");		 
		 }
		
		 return c;
	}
	
	
	// adapted from https://x-engineer.org/bisection-method/
	public double bisectionalSolver(double a, double b, double D_theta, double psi, double K, double dt, double F0)
	{
		double localA = a;
		double localB = b;
		
		double NMAX = ITERATIONS; /* maximum number of iterations */
		double c = 0; /* estimated root */
		int index = 0; /* index */
//		int stuckCount = 0;
		 
		c = (localA + localB)/2.0; /* calculate the midpoint */
		 /* Evaluate loop until the result is less than the tolerance
		 * maximum number of iterations is not yet reached*/
		
		 double[] parametersC = new double[]{D_theta ,psi, K, dt, F0, c};
		 double[] parametersA = new double[]{D_theta ,psi, K, dt, F0, localA};
		 if (calculateExpression2(parametersC) == 0)
		 {
		  /* If the first midpoint gives f(c) = 0, c is the root */
	//		 System.out.println("root is " + c);
		 }
		 else
		 {
			 double calEx2PCInitial = Math.abs(calculateExpression2(parametersC));
			 while ( (calEx2PCInitial > TOL || Double.isNaN(calEx2PCInitial) ) && (index<=NMAX))
			 {
				 double cex2PC = sign(calculateExpression2(parametersC));
				 double cex2PA = sign(calculateExpression2(parametersA));
				 if (cex2PC == cex2PA)
				 {
					 /* f(c) has same sign as f(a) */
					 localA = c;
				 }
				 else
				 {
					 /* f(c) has same sign as f(b) */
					 localB = c;
				 } 
				 c = (localA+localB)/2.0; /* midpoint update */
				 double difference = Math.abs(localA-localB);
				 if (DEBUG)
				 {
					 System.out.println( (index+1) + " | " + cex2PC + " | " + cex2PA + " | " + localA + " | " + localB + " | " + difference + " | " + c );
				 }
				
				 //did it get stuck? Try different variations of end points
//				 if (difference < STUCK_VALUE  )
//				 {					 
//					 if (stuckCount > 1)
//					 {
//						 localA = -1. - Math.random();
//						 localB = 1. + Math.random();
//					 }					 
//					 if (stuckCount > 2)
//					 {
//						 localA = 0.0-Math.random();
//						 localB =  0.0;
//					 }
//					 if (stuckCount > 3)
//					 {				 
//						 if (stuckCount % 2 == 0)
//						 {
//							 localA = -0.001;
//							 localB =  0.0-Math.random()/1000.;
//						 }
//						 else
//						 {
//							 localA = -0.001;
//							 localB =  Math.random();
//						 }
//					 }
//					 if (stuckCount < 1)
//					 {
//						 c = c - 400;
//					 } 
//					 stuckCount ++;
//				 }
				 parametersA = new double[]{D_theta ,psi, K, dt, F0, localA};
				 parametersC = new double[]{D_theta ,psi, K, dt, F0, c};
	
				 index++; /* index increment */
				 calEx2PCInitial = Math.abs(calculateExpression2(parametersC));
			 }
		 } 
		 
		 if (index>=NMAX)
		 {
			 if (DEBUG)
			 {
				 System.out.println("Root not found " + c + " after " + index + " iterations");	
			 }
			 
			 return ERROR_RETURN;
		 }
		 
		 /* Display results */
		 if (DEBUG)
		 {
			 System.out.println("Root is " + c + " found after " + index + " iterations");		 
		 }
		
		 return c;
	}
	
	public double calculateExpression2(double[] parameters)
	{
		double D_theta = parameters[0];
		double psi = parameters[1];
		double K = parameters[2];
		double dt = parameters[3];
		double F0 = parameters[4];
		double F = parameters[5];
		
		
//		double returnValue = Math.pow(x, 3) + Math.pow(x, 2) - (3 * x) - 3;		
//		F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
		
		double F_psi_D_theta = F+psi*D_theta;
		double F0_psi_D_theta = F0+psi*D_theta;
		
//		double returnValue = F-F0-psi*D_theta*Math.log((F_psi_D_theta)/(F0_psi_D_theta))-K*(dt);
//		double returnValue = F-F0-psi*D_theta*Math.log((F_psi_D_theta)/(F0_psi_D_theta))-K*(dt);
		double returnValue = (F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt));
		
		return returnValue;
	}
	
	public double calculateExpression(double x)
	{
		double returnValue = Math.pow(x, 3) + Math.pow(x, 2) - (3 * x) - 3;		
//		F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
		return returnValue;
	}
	
	
//	public double calculateExpression()
//	{
//		double returnValue;
//		
//		returnValue = (F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt);
//		
////		F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
//		
//		return returnValue;
//	}
	
	public double bisectionTest()
	{
		
//		double dz=-1.0;
//		double ref_ta = 21.9;
//		double UTb = 2.63038178;
//		double mod_U_TaRef = 3.05879268;
////		int i=0;
//		double Ri_rur = 0.24555776;
////		double testValue = c;
//		
////		1 -1.0 21.9 2.63038178 3.05879268 0.24555776
////		21.7995777279967
////		double returnValue = calculateExpression(dz, ref_ta, UTb, mod_U_TaRef, Ri_rur, c);
		
		
		
		
		
		double dz=15.0 ;
		double ref_ta =23.5 ;
		double UTb =  0.21266443039771887 ;
//		double[] mod_U_TaRef = new double[1];
		double mod_U_TaRef =  0.1  ;
		double Ri_rur = -370.55369919381894;
		//TODO this should be 22.760359533246
		
		
		double a = -20; /* left point */
		double b = 50; /* right point */
		double TOL = 1.0E-12; /* tolerance */
		double NMAX = 1000; /* maximum number of iterations */
		double c = 0; /* estimated root */
		int i = 0; /* index */
		 
		c = (a + b)/2; /* calculate the midpoint */
		 
		 /* Display the table header and initial data */
		// printf("i\ta\t\tb\t\tc\t\tf(a)\t\tf(b)\t\tf(c)\n");
		// f_disp(i, a, b, c);
		 
		 /* Evaluate loop until the result is less than the tolerance
		 * maximum number of iterations is not yet reached*/
		 
		 if (calculateExpression(dz, ref_ta, UTb, mod_U_TaRef, Ri_rur, c) == 0)
		 {
		  /* If the first midpoint gives f(c) = 0, c is the root */
//		  printf("Root is: %f \n", c);
			 System.out.println("root is " + c);
		 }
		 else
		 {
			 while ((Math.abs(calculateExpression(dz, ref_ta, UTb, mod_U_TaRef, Ri_rur, c)) > TOL) && (i<=NMAX))
			 {
				 if (sign(calculateExpression(dz, ref_ta, UTb, mod_U_TaRef, Ri_rur, c)) == sign(calculateExpression(dz, ref_ta, UTb, mod_U_TaRef, Ri_rur, a)))
				 {
					 /* f(c) has same sign as f(a) */
					 a = c;
				 }
				 else
				 {
					 /* f(c) has same sign as f(b) */
					 b = c;
				 } 
				 c = (a+b)/2; /* midpoint update */
				 //				    f_disp(i+1, a, b, c); /* display current data */
				 System.out.println( (i+1) + " " + a + " " + b + " " + c);
				 i++; /* index increment */
			 }
		 }
		 
		 /* Display results */
//		 printf("\nRoot is c=%.7f found after %d iterations\n", c, i);
		 System.out.println("Root is " + c + " found after " + i + " iterations");
//		 printf("The value of the function f(c) is: %.10f\n", f(c));
//		 return 0;
		 
		 return c;
	}
	

	public double bisection()
	{
		double a = -2; /* left point */
		double b = 5; /* right point */
		double TOL = 1E-6; /* tolerance */
		double NMAX = 1000; /* maximum number of iterations */
		double c = 0; /* estimated root */
		int i = 0; /* index */
		 
		c = (a + b)/2; /* calculate the midpoint */
		 
		 /* Display the table header and initial data */
		// printf("i\ta\t\tb\t\tc\t\tf(a)\t\tf(b)\t\tf(c)\n");
		// f_disp(i, a, b, c);
		 
		 /* Evaluate loop until the result is less than the tolerance
		 * maximum number of iterations is not yet reached*/
		

		 
		 if (f(c) == 0)
		 {
		  /* If the first midpoint gives f(c) = 0, c is the root */
//		  printf("Root is: %f \n", c);
			 System.out.println("root is " + c);
		 }
		 else
		 {
			 while ((Math.abs(f(c)) > TOL) && (i<=NMAX))
			 {
				 if (sign(f(c)) == sign(f(a)))
				 {
					 /* f(c) has same sign as f(a) */
					 a = c;
				 }
				 else
				 {
					 /* f(c) has same sign as f(b) */
					 b = c;
				 } 
				 c = (a+b)/2; /* midpoint update */
				 //				    f_disp(i+1, a, b, c); /* display current data */
				 System.out.println( (i+1) + " " + a + " " + b + " " + c);
				 i++; /* index increment */
			 }
		 }
		 
		 /* Display results */
//		 printf("\nRoot is c=%.7f found after %d iterations\n", c, i);
		 System.out.println("Root is " + c + " found after " + i + " iterations");
//		 printf("The value of the function f(c) is: %.10f\n", f(c));
//		 return 0;
		 
		 return c;
	}
	
	public double calculateExpressionBi(double dz, double ref_ta, double UTb, double mod_U_TaRef, double Ri_rur, double Thi_tb)
	{
		
		double expressionValue = 9.806 * dz *(Thi_tb- ref_ta)*2.0/(Thi_tb + ref_ta )/ Math.pow((UTb-mod_U_TaRef),2.0) -  Ri_rur;
		
		return expressionValue;
	}


	public double sign(double x)
	{
		return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
	}
	
	public double f(double x)
	{
		double y;
		y = 10 - x*x;
		return y;
	}
	
	public void test5() 
	{
		double F00=0.5;
		double F0=F00;
		double psi=1.11;
		double K=1.5;
		double dt=1.4;
		double Se=1.2;
		double theta_e=1.1;
		double AWI=10.0;
		double F = 0;
		
		double Dtheta = (1.-Se)*theta_e; 
		double I=AWI/dt;
		double f = K*(Math.max(psi*Dtheta/F0,0)+1);
		double psiDTheta = psi*Dtheta;
		
		double a = F;
		double b= F0;
		double c=psiDTheta;
		double g=K;
		double j=dt;
		
		String f1 = "a-b-c*log((a+c)/(b+c))-g*j";

		List<String> functions = new ArrayList<String>();
		functions.add(f1);

		List<String> vars = new ArrayList<String>();
		vars.add("a");
		vars.add("b");
		vars.add("c");
		vars.add("g");
		vars.add("j");


		
		List<Double> xo = new ArrayList<Double>();

		
		xo.add(a);
		xo.add(b);
		xo.add(c);
		xo.add(g);
		xo.add(j);

		
		FunctionXs fXs1 = new FunctionXs(f1);

		
//		List<Double> test = new ArrayList<Double>();
//		test.add(14.1355467);
//		test.add(10.1303043);
//		test.add(43.9596517);
		
        try 
        {
			NonLinearSys nls = new NonLinearSys(functions,vars,xo);
			List<Double> result = nls.calc(0.1e-6, 1000);
			for (Double d : result) 
			{
				System.out.println(d);				
			}
			System.out.println("Test1:"+fXs1.getValue(result, vars));


			
		} catch (Exception e) 
        {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	public void test4()
	{
		double[] x ={1,-5,6};
		
	      try {
	    
			Polynomial.rootCalc(x);
			System.out.println("***************");
		} catch (Exception e)
	     {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public void test3() {
		String f1 = "-6x^2+26x-59";
//		String f2 = "4.67*e^(-3)*x^1.75+20-z";
//		String f3 = "3.72*e^(-2)*y^1.75+15-z";
		List<String> functions = new ArrayList<String>();
		functions.add(f1);
//		functions.add(f2);
//		functions.add(f3);
		List<String> vars = new ArrayList<String>();
		vars.add("x");
//		vars.add("y");
//		vars.add("z");
		List<Double> xo = new ArrayList<Double>();
		xo.add(16.0);
//		xo.add(7.0);
//		xo.add(50.0);
		
		FunctionXs fXs1 = new FunctionXs(f1);
//		FunctionXs fXs2 = new FunctionXs(f2);
//		FunctionXs fXs3 = new FunctionXs(f3);
		
//		List<Double> test = new ArrayList<Double>();
//		test.add(14.1355467);
//		test.add(10.1303043);
//		test.add(43.9596517);
		
        try 
        {
			NonLinearSys nls = new NonLinearSys(functions,vars,xo);
			List<Double> result = nls.calc(0.1e-6, 1000);
			for (Double d : result) 
			{
				System.out.println(d);				
			}
			System.out.println("Test1:"+fXs1.getValue(result, vars));
//			System.out.println("Test2:"+fXs2.getValue(result, vars));
//			System.out.println("Test3:"+fXs3.getValue(result, vars));

			
		} catch (Exception e) 
        {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	public void test() {
		String f1 = "2.35*e^(-3)*(x+y)^1.75-75+z";
		String f2 = "4.67*e^(-3)*x^1.75+20-z";
		String f3 = "3.72*e^(-2)*y^1.75+15-z";
		List<String> functions = new ArrayList<String>();
		functions.add(f1);
		functions.add(f2);
		functions.add(f3);
		List<String> vars = new ArrayList<String>();
		vars.add("x");
		vars.add("y");
		vars.add("z");
		List<Double> xo = new ArrayList<Double>();
		xo.add(16.0);
		xo.add(7.0);
		xo.add(50.0);
		
		FunctionXs fXs1 = new FunctionXs(f1);
		FunctionXs fXs2 = new FunctionXs(f2);
		FunctionXs fXs3 = new FunctionXs(f3);
		
		List<Double> test = new ArrayList<Double>();
		test.add(14.1355467);
		test.add(10.1303043);
		test.add(43.9596517);
		
        try 
        {
			NonLinearSys nls = new NonLinearSys(functions,vars,xo);
			List<Double> result = nls.calc(0.1e-6, 1000);
			for (Double d : result) {
				System.out.println(d);				
			}
			System.out.println("Test1:"+fXs1.getValue(result, vars));
			System.out.println("Test2:"+fXs2.getValue(result, vars));
			System.out.println("Test3:"+fXs3.getValue(result, vars));

			
		} catch (Exception e) 
        {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public void test2()
	{
		double dz=15.0 ;
		double ref_ta =16.6 ;
		double UTb =  1.8069344530637708 ;
		double mod_U_TaRef = 0.8496646334718477 ;
		double Ri_rur = -6.956094923672547;
//		Result : 15.8958756964009

		
		
		String f1 = "9.806*dz*(Thi_tb-ref_ta)*2.0/(Thi_tb+ref_ta)/(UTb-mod_U_TaRef)^2.0-Ri_rur";
		
		String f2 = "9.806*dz*(Thi_tb-ref_ta)*2.0/(Thi_tb+ref_ta)/(UTb-mod_U_TaRef)^2.0-Ri_rur";
		
		List<String> functions = new ArrayList<String>();
		functions.add(f1);
		List<String> vars = new ArrayList<String>();
		vars.add("dz");
		
		vars.add("ref_ta");
		vars.add("UTb");
		vars.add("mod_U_TaRef");
		vars.add("Ri_rur");
//		vars.add("Thi_tb");
		List<Double> xo = new ArrayList<Double>();
		xo.add(dz);
		
		xo.add(ref_ta);
		xo.add(UTb);
		xo.add(mod_U_TaRef);
		xo.add(Ri_rur);
//		xo.add(Thi_tb);
		
		FunctionXs fXs1 = new FunctionXs(f1);
		
		List<Double> test = new ArrayList<Double>();
		test.add(14.1355467);
		test.add(10.1303043);
		test.add(43.9596517);
		
        try 
        {
			NonLinearSys nls = new NonLinearSys(functions,vars,xo);
			List<Double> result = nls.calc(0.1e-6, 1000);
			for (Double d : result) 
			{
				System.out.println(d);				
			}
			System.out.println("Test1:"+fXs1.getValue(result, vars));

			
		} catch (Exception e) 
        {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public double calculateExpression(double dz, double ref_ta, double UTb, double mod_U_TaRef, double Ri_rur, double Thi_tb)
	{
		
		double expressionValue = 9.806 * dz *(Thi_tb- ref_ta)*2.0/(Thi_tb + ref_ta )/ Math.pow((UTb-mod_U_TaRef),2.0) -  Ri_rur;
		
		return expressionValue;
	}
	
	
//	public void test()
//	{
//		 NewtonRaphsonSolver test = new NewtonRaphsonSolver();
//	        UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() 
//	        {
//
//	            public double value(double x) 
//	            {
//	                return Math.cos(x);
//	            }
//	            
////	        	public double value(double dz, double ref_ta, double UTb, double mod_U_TaRef, double Ri_rur, double Thi_tb)
////	        	{
////	        		
////	        		double expressionValue = 9.806 * dz *(Thi_tb- ref_ta)*2.0/(Thi_tb + ref_ta )/ Math.pow((UTb-mod_U_TaRef),2.0) -  Ri_rur;
////	        		
////	        		return expressionValue;
////	        	}
//
//	            @Override
//	            public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException 
//	            {
//	                return t.cos();
//	            }
//
////				@Override
////				public double value(double arg0)
////				{
////					// TODO Auto-generated method stub
////					return 0;
////				}
//	        };
//
//	        for (int i = 1; i <= 500; i++) 
//	        {
//	            System.out.println(test.solve(1000, f, i, i+0.1));
//	        }
//	       
//	}
	

}
