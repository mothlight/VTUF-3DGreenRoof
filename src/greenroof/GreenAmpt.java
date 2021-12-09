package greenroof;

public class GreenAmpt
{

	public void Green_Ampt_(double F0,double psi,double K,double dt,double Se,double theta_e,double AWI)
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
	double D_theta = (1-Se)*theta_e;  
	//Intensity of Available water to infiltrate
	double I=AWI/dt;
	
	//Initial infiltration rate
	f = K*(Math.max(psi*D_theta/F0,0)+1);
	
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
	    double f_t = (K*(psi*D_theta/F_t+1));
	
	    
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
	        F = fsolve(@(F)(F-F_p-psi*D_theta*Math.log((F+psi*D_theta)/(F_p+psi*D_theta))-K*(dt-dt_p)),F0+I*dt); //F0+I*dt
	        f = K*(psi*D_theta/F+1);
	    }
	}
	//Assumption 5: If I>f, ponding is possible at the beggining of time step
	else
	{
	        F = fsolve(@(F)(F-F0-psi*D_theta*Math.log((F+psi*D_theta)/(F0+psi*D_theta))-K*(dt)),F0);
	        f = real(K*(psi*D_theta/F+1));
	        F=(F);
	        f=f(1);
	}
	    
	
	
	return
	}
    
    
    
}
