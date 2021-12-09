package greenroof;

public class EPWWeatherDataLine
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
        public double IR;     // Infrared radiation           
        public double Tsky; // Sky temperature
        

        public void EPWWeatherDataLine_(double[] data, double Z, double SiteWindBLHeight, double SiteWindExp)
        {
            
             double WeatherFileWindModCoeff = 1.5863;
             
             // SKIPPED
             
             // --. N1, \field Year
             // --. N2, \field Month
             // --. N3, \field Day
             // --. N4, \field Hour
             // --. N5, \field Minute
             // --. A1, \field Data Source and Uncertainty Flags
             
            // 1. N6, \field Dry Bulb Temperature 
            Tair = data[Constants.ONE]+273.15;
            // 2. N7, \field Dew Point Temperature                   
            // 3. N8, \field Relative Humidity  
            RH = data[Constants.THREE]/100;
            
            Pa = data[Constants.FOUR];// 4. N9, \field Atmospheric Station Pressure
            // 5. N10, \field Extraterrestrial Horizontal Radiation
            // 6. N11, \field Extraterrestrial Direct Normal Radiation
            IR = data[Constants.SEVEN];// 7. N12, \field Horizontal Infrared Radiation Intensity            
            R_sh = data[Constants.EIGHT]; // 8. N13, \field Global Horizontal Radiation            
            // 9. N14, \field Direct Normal Radiation
            // 10. N15, \field Diffuse Horizontal Radiation
            // 11. N16, \field Global Horizontal Illuminance
            // 12. N17, \field Direct Normal Illuminance
            // 13. N18, \field Diffuse Horizontal Illuminance                
            // 14. N19, \field Zenith Luminance
            // 15. N20, \field Wind Direction
            U = data[Constants.SIXTEEN]*WeatherFileWindModCoeff*Math.pow((Z/SiteWindBLHeight),SiteWindExp);// 16. N21, \field Wind Speed
            // 17. N22, \field Total Sky Cover
            // 18. N23, \field Opaque Sky Cover (used if Horizontal IR Intensity missing)
            // 19. N24, \field Visibility
            // 20. N25, \field Ceiling Height
            // 21. N26, \field Present Weather Observation
            // 22. N27, \field Present Weather Codes
            // 23. N28, \field Precipitable Water
            // 24. N29, \field Aerosol Optical Depth
            // 25. N30, \field Snow Depth
            // 26. N31, \field Days Since Last Snowfall
            // 27. N32, \field Albedo
            rainfall = data[Constants.TWENTYEIGHT];// 28. N33, \field Liquid Precipitation Depth
            // 29. N34 \field Liquid Precipitation Quantity
            
            //// Not sure                      
            irrigation = 0; 
            
            Tsky = Math.pow((IR/Constants.SB),0.25);
        }
    
    

}
