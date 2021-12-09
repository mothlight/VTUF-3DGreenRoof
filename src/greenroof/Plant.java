package greenroof;

public class Plant
{

	public double rho=.06;     //sw reflectivity plant
	public double em=0.98;    //emissivity planta
	public double ks=0.14; // Extintion coefficient   
	public double ks_ir = 0.14; // IR extintion coefficient
	public double LAI=4;        //leaf area index A FEB Tabares  
	        
	public double k=0.5;      // thermal conductivity planta
	public double rsmin=500;      //[s/m] resistencia estomatica 
	public double height=0.55;   //  height of plant [m]
	public double z=2;       // [m] height of wind measurements
	        
	public double Zog = 0.001; // VARIA SEGUN SMOOOTH
	        
	//% Sailor requirements
	//fc=0.01;    // fractional vegetation coverage %
	                
    
	               
	       
	 public double fc()
	 {
	         double   s = 0.9 - 0.7*Math.exp(-0.75*LAI);
	         return s;
	 }
	        
	 public double tau_fsol()
	 {
	         double r = Math.exp(-ks*LAI);
	         return r;
	}
	       
    public double  tau_fir()
    {
       double r = Math.exp(-ks_ir*LAI);
       return r;
    }
    
    public double foliage_rho()
    {
       double r = (1-tau_fsol())*rho; 
       return r;
    }
	   
	
}
