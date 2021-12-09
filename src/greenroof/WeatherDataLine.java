package greenroof;

public class WeatherDataLine
{

        //// measured variables
        public double Tair; // Air temperature
        public double RH; // Relative Humidity
        public double U; // Wind speed
        public double rainfall; //        
        public double R_sh; //Incomming solar Radiation
        //VWC // Volumetric water content in soil
        //Tin // Interior temperature (below the ceiling)
        public double Pa; // Atmospheric pressure
        public double irrigation;
        public double IR; // Infrared radiation
        //// derived variables
        public double Tsky; // Sky temperature
        

        

   public WeatherDataLine(double[] data)
   {
            
            //c = Constants;
            
            //// Measured
	   rainfall = data[Constants.ONE]; // Rainfall             
	   R_sh = data[Constants.TWO]; // Solar irradiance
	   Tair = data[Constants.THREE] + 273.15; // Exterior Air Temperature
	   RH = data[Constants.FOUR]/100.; // Exterior air RH 
	   U = data[Constants.FIVE]; // Wind speed 
	   Pa = 101300.; //Pa... assumed constant                        
	   irrigation = data[Constants.SIX]; //data(8)/10;

                       
            
        //// Derived
  
            // Assume clear sky
        double b=18.678;
        double c = 275.14;
        double d = 234.5;
        double g = Math.log(RH*Math.exp((b-(Tair-273.15)/d)*((Tair-273.15)/(c+(Tair-273.15)))));            
        double dewPointT = c*g/(b-g);
        double skyEmissivity = 0.787+0.764*Math.log((dewPointT+273.15)/274);
        IR = skyEmissivity*Constants.SB*Math.pow(Tair,4);
        Tsky = Math.pow((IR/Constants.SB),0.25);        
        
//        double[] newDataLine = new double[]{rainfall,R_sh,Tair,RH,U,Pa,irrigation,IR,Tsky};
//        return newDataLine;
        
   }
    

}
