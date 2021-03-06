package greenroof;

import java.util.TreeMap;

public class Interception
{
	//INTERCEPTION FROM VEGETATION OR STORAGE SURFACE

	public TreeMap<String,Double> Interception_(double P,double[] ET0,double P_cum_event_bf,int[] time_events,double k,double S,double dt)
	{
	//INPUT
	//P                              Precipitation(mm)
	//ET0                            Evapotranspiration (mm/h)
	//P_cum_event_bf                 Cumulative effective precipitation until the beginning of the time interval dt (mm)
	//time_events                    Indicates the amount of dt from the beginning of the events
	//k                              Leaf area index
	//S                              Maximum water stored (mm)
	//dt                             Modeling time interval (h)

	//OUTPUT
	//Peffect                        Precipitation minus interception (mm)
	//P_cum_event                    Cumulative effective precipitation at the end the time interval dt (mm)

	//OTHER VARIABLES
	//Evap                           Evaporation during the time interval dt (mm)
	//Evap_cum_event                 Cumulative evaporation until time [t t-1]
	//S_interval
	//I_interval
	//S_cum
	//S_cum_bf
	double Evap;
	double[] Evap_cum_event=new double[2];
	double S_interval;
	double I_interval;
	double S_cum;
	double S_cum_bf;
	double P_cum_event;
	double Peffect;

	//NOTE
	//time_events*dt is the time (in hour) from the beginning of the event
	//ET0 and time_events are vectors with two values, (1) is the value at time (t) 
	//and (2) is the value at time (t-1)

	    //Evaporation from the leaf surface or storage surface
	    if (time_events[Constants.ONE]==Constants.ONE) //At the beginning of the event
	    {
	        //Cumulative evaporation until time t
	        Evap_cum_event[Constants.ONE]=ET0[Constants.ONE]*time_events[Constants.ONE]*dt*k;
	        //Evaporation during the interval
	        Evap=Evap_cum_event[Constants.ONE];      
	    }
	    else    //During the event
	    {
	        //Cumulative evaporation until time t
	        Evap_cum_event[Constants.ONE]=ET0[Constants.ONE]*time_events[Constants.ONE]*dt*k;
	        //Cumulative evaporation until time t-1
	        Evap_cum_event[Constants.TWO]=ET0[Constants.TWO]*time_events[Constants.TWO]*dt*k; 
	        //Evaporation during the interval
	        Evap=Math.max(0,Evap_cum_event[Constants.ONE]-Evap_cum_event[Constants.TWO]);
	    }

	    //Cumulative precipitation
	    if (P==0) //If the rain is null
	    {
	        //Some of the cumulative precipitation evaporates
	        P_cum_event=Math.max(P_cum_event_bf-Evap,0);
	    }
	    else //If it rains accumulated water increases
	    {
	        P_cum_event=P+P_cum_event_bf;
	    }
	    
	    //Cumulative leaf storage or cumulative surface storage
	    //Cumulative storage at t-1 (beginning of dt)
	    S_cum_bf=S*(1-Math.exp(-P_cum_event_bf/S));  
	    //Cumultive storage at t (ending of dt)
	    S_cum=S*(1-Math.exp(-P_cum_event/S));          
	    
	    //Leaf storage or surface storage during dt
	    if (P==0) //If the rain is null, no water is added to storage
	    {
	       S_interval=0;
	    }
	    else //If it is rainy, water is added to storage
	    {
	       S_interval=S_cum-S_cum_bf;
	    }

	    //Intercepted water is equal to storage water plus evaporated water during dt
	    I_interval=Math.max(Evap+S_interval,0);
	    //Effective precipitation is equal to precipitation minus intercepted water
	    Peffect=Math.max(P-I_interval,0); 
	    
	    TreeMap<String,Double> returnValues = new TreeMap<String,Double>();
	    returnValues.put("Peffect", Peffect);
	    returnValues.put("Ptot_cum_event", P_cum_event);
	    return returnValues;
	    
 }
}
