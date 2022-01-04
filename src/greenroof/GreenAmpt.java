package greenroof;

import java.util.TreeMap;

public class GreenAmpt
{
	boolean DEBUG=true;
	double TOL = 1.0E-6; /* tolerance */
	public static final double ERROR_RETURN = -9999.;
	public static final int ITERATIONS = 10000;

	public TreeMap<String,Double> Green_Ampt_(double F0,double psi,double K,double dt,double Se,double theta_e,double AWI)
	{
	
	//INPUT PARAMETERS
	//F0         	Cumulative infiltration at the beginning of the time step
	//psi            Bubbling pressure (mm)
	//K              Saturated hydraulic conductivity (mm/h)
	//dt             Time step (h)
	//Se             Relative saturation
	//theta_e        Effective soil water content
	//AWI            Available water to infiltrate (mm)
	
	//OUTPUT PARAMETERS
	//f              Infiltration rate (mm/h)
	//F              Cumulative infiltration (mm)
	double f;
	double F;
	//Soil water content change at the wetting front
	double D_theta = (1.0-Se)*theta_e;  
	//Intensity of Available water to infiltrate
	double I=AWI/dt;
	
	//Initial infiltration rate
	f = K*(Math.max(psi*D_theta/F0,0)+1.0);
	F=0;
	
	//Assumption 1: When AWI=0, then infiltration will be null
	if (I==0)
	{
	    f=0;
	    F=F0;
	}
	  
	//Assumption 2: If I<f, no ponding occurs at the beginning of the time step    
	else if (f>I)
	{
	    //Possible values
	    double F_t = F0+I*dt;
	    double f_t = (K*(psi*D_theta/F_t+1.0));
	
	    
	    //Assumption 3: If I<f_t, no ponding occurs during the time step 
	    if (f_t>I)
	    {
	        F = F_t;
	        f = I; //Infiltration rate is equal to intensity of AWI
	    }
	    
	    //Assumption 4: If I>f_t, ponding is possible during the time step
	    else 
	    {
	        double F_p = K*psi*D_theta/(I-K);
	        double dt_p = (F_p-F0)/I;
			double a = 0.25;
			double b = 10000 ;
	        double value = bisectionalSolverAssm4(a, b, D_theta, psi, K, dt, F_p, dt_p);
	        F=value;
//	        double x0=F0+I*dt;
//	        F=(F-F_p-psi*D_theta*Math.log((F+psi*D_theta)/(F_p+psi*D_theta))-K*(dt-dt_p));
	        
//	        F = fsolve(@(F)(F-F_p-psi*D_theta*Math.log((F+psi*D_theta)/(F_p+psi*D_theta))-K*(dt-dt_p)),F0+I*dt); //F0+I*dt
	        f = K*(psi*D_theta/F+1.0);
	    }
	}
	//Assumption 5: If I>f, ponding is possible at the begining of time step
	else
	{
//		double x0=F0;
//		F=(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt));
		
		double a = 0.25;
		double b = 10000 ;
//		FSolve fsolve = new FSolve();
//		double value = fsolve.bisectionalSolver(a, b, D_theta, psi, K, dt, F0);		
		double value = bisectionalSolverAssm5(a, b, D_theta, psi, K, dt, F0);		
//	        F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
		F=value;
		f = (K*(psi*D_theta/F+1));
		
//	        f=f[GreenRoofCommon.ONE];
	}
	    
	
	// return [f[i],Ft[i]]
	TreeMap<String,Double> returnValues = new TreeMap<String,Double>();
	returnValues.put("f", f);
	returnValues.put("Ft", F);
	return returnValues;
	}
    
    
	// adapted from https://x-engineer.org/bisection-method/
	public double bisectionalSolverAssm5(double a, double b, double D_theta, double psi, double K, double dt, double F0)
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
		 if (calculateExpressionAssm5(parametersC) == 0)
		 {
		  /* If the first midpoint gives f(c) = 0, c is the root */
	//		 System.out.println("root is " + c);
		 }
		 else
		 {
			 double calEx2PCInitial = Math.abs(calculateExpressionAssm5(parametersC));
			 while ( (calEx2PCInitial > TOL || Double.isNaN(calEx2PCInitial) ) && (index<=NMAX))
			 {
				 double cex2PC = sign(calculateExpressionAssm5(parametersC));
				 double cex2PA = sign(calculateExpressionAssm5(parametersA));
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
				 calEx2PCInitial = Math.abs(calculateExpressionAssm5(parametersC));
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
	
	public double calculateExpressionAssm5(double[] parameters)
	{
		double D_theta = parameters[0];
		double psi = parameters[1];
		double K = parameters[2];
		double dt = parameters[3];
		double F0 = parameters[4];
		double F = parameters[5];
		
		
//		double returnValue = Math.pow(x, 3) + Math.pow(x, 2) - (3 * x) - 3;		
//		F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
		
//		double F_psi_D_theta = F+psi*D_theta;
//		double F0_psi_D_theta = F0+psi*D_theta;
		
//		double returnValue = F-F0-psi*D_theta*Math.log((F_psi_D_theta)/(F0_psi_D_theta))-K*(dt);
//		double returnValue = F-F0-psi*D_theta*Math.log((F_psi_D_theta)/(F0_psi_D_theta))-K*(dt);
		double returnValue = (F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt));
		
		return returnValue;
	}
	
	
	public double calculateExpressionAssm4(double[] parameters)
	{
		double D_theta = parameters[0];
		double psi = parameters[1];
		double K = parameters[2];
		double dt = parameters[3];
		double F_p = parameters[4];
		double dt_p = parameters[5];
		double F = parameters[6];
		
//		 F = fsolve(@(F)(F-F_p-psi*D_theta*Math.log((F+psi*D_theta)/(F_p+psi*D_theta))-K*(dt-dt_p)),F0+I*dt); //F0+I*dt
		double returnValue = (F-F_p-psi*D_theta*Math.log((F+psi*D_theta)/(F_p+psi*D_theta))-K*(dt-dt_p));
		
		return returnValue;
	}
	
	
	// adapted from https://x-engineer.org/bisection-method/
	public double bisectionalSolverAssm4(double a, double b, double D_theta, double psi, double K, double dt, double F_p, double dt_p)
	{
		double localA = a;
		double localB = b;
		
		double NMAX = ITERATIONS; /* maximum number of iterations */
		double c = 0; /* estimated root */
		int index = 0; /* index */
		 
		c = (localA + localB)/2.0; /* calculate the midpoint */
		 /* Evaluate loop until the result is less than the tolerance
		 * maximum number of iterations is not yet reached*/
		
		 double[] parametersC = new double[]{D_theta ,psi, K, dt, F_p, dt_p, c};
		 double[] parametersA = new double[]{D_theta ,psi, K, dt, F_p, dt_p, localA};
		 if (calculateExpressionAssm4(parametersC) == 0)
		 {
		  /* If the first midpoint gives f(c) = 0, c is the root */
	//		 System.out.println("root is " + c);
		 }
		 else
		 {
			 double calEx2PCInitial = Math.abs(calculateExpressionAssm4(parametersC));
			 while ( (calEx2PCInitial > TOL || Double.isNaN(calEx2PCInitial) ) && (index<=NMAX))
			 {
				 double cex2PC = sign(calculateExpressionAssm4(parametersC));
				 double cex2PA = sign(calculateExpressionAssm4(parametersA));
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
				 parametersA = new double[]{D_theta ,psi, K, dt, F_p, dt_p, localA};
				 parametersC = new double[]{D_theta ,psi, K, dt, F_p, dt_p, c};
	
				 index++; /* index increment */
				 calEx2PCInitial = Math.abs(calculateExpressionAssm4(parametersC));
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
	
	public double sign(double x)
	{
		return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
	}
}
