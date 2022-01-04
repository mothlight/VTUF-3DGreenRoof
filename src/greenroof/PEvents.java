package greenroof;

public class PEvents
{

//EVENTS DURATION
//Efective precipitation events are separated by a certain number of hours

	public int[] p_events_(double[] Ptotal, int t_s, int dt)
	{
	
		//INPUT
		//Ptotal         Total precipitacion (mm)
		//t_s            Time interval between events (h)
		
		//VARIABLE INITIATION 
		//time_events indicates the amount of dt from the beginning of the events
		int[] time_events=new int[Ptotal.length];
		
		//The first value of the precipitation series starts the first event
		for (int i=Constants.ONE;i<Math.min(t_s/dt,Ptotal.length);i++)
		{
		    time_events[i]=i;
		}
		
		//When the time between two records of total precipitation is greater than
		//t_s, then that events ends and another begins
		int counter=t_s/dt;
		for (int i=counter+1-1;i<Ptotal.length;i++)
		{
		    if (Ptotal[i]>0)
		    {
		        int events0=0;
		        int k=Constants.ONE;
		        for (int j=Constants.ONE;j<counter;j++)
		        {
		            if (Ptotal[i-k]>0)
		            {
		                events0=events0+1;
		            }
		            k=k+1;
		        }
		        if (events0>0)
		        {
		            time_events[i]=time_events[i-1]+1;
		        }
		        else
		        {
		            time_events[i]=1;
		        }
		    }
		    else
		    {
		        time_events[i]=time_events[i-1]+1;
		    }
		}
		
		
		return time_events;
	
	
	
	}


}
