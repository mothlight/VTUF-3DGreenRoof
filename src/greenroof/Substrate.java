package greenroof;

public class Substrate
{
    
	public double rho=0.1853;     //sw reflectivity substrate
	public double depth=0.10;     // depth substrate  (ex t)
	public double Lsubm=0.05;    // depth substrate middepth sensor
	public double em = 0.96;
	public double dens =1100; //kg/m^3
	public double Cp=1000; // J/kgK;
	public double k=0.13; //W/mK;
	public double VWCsat=0.8;     // VWC saturated
	public double VWCfc = 0.35;
	public double VWCwilting=0.02;    // VWC wilting point
	public double VWCresidual = 0.02;
	public double phi=0.85;           

    
  public double thermalConductivity(double VWC)
  {
           //k = (obj.k + VWC*0.591); 
           double k_ = k+VWC*0.5811;
           //k = (obj.k+VWC*0.7845);
           return k_;
  }
       
  public double density(double VWC)
  {
	  double rho = dens + VWC*1000;  
	  return rho;
  }
  
  public double rhoCp(double VWC)
  {
	  //r = obj.density(VWC);     
	  // Alexantri, Jones, 2006
	  double rcp = (1-VWCsat)*dens*Cp + VWC * 1000*4182;       
	  return rcp;
  }
       

}
