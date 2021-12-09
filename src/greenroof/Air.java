package greenroof;

public class Air
{
	public static final double Cp=1005.6;     //% specific heat air
	public static final double R=286;      //% gas constant of air [J/kg K]
	public static final double dens=40;  //% molar density of air
	public static final double k=0.02554;    //% aire [W/m K]
	public static final double gamma=0.06884;        
	public static final double Kv=0.4;
	public static final double Sch=0.63;      //% turbulent Schmidt number
	public static final double Pr=0.71;      //% turbulent Prandtl number
	
	
    public double density(double Pa, double Tair)
    {
            double d = Pa/(Air.R*Tair); 
            return d;
    }
        
    public double kinematicViscosity(double T)
    {
        //% Function obtained from fitting data 
        //% available in http://www.engineeringtoolbox.com/dry-air-properties-d_973.html           
        double mu = 5.409413e-11*T*T+7.490993e-8*T-1.163350e-5;
        return mu;
    }
}





