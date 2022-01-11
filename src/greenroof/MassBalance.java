package greenroof;

import java.util.TreeMap;

public class MassBalance
{
	//HYDROLOGICAL PROCESSES FOR PERVIOUS subarea

//	function [theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event]= 
		public TreeMap MassBalance_(int tstep, double[] P, double[] R, double[] runon, double[] E_T, int[] general_inf, double[] outflow_data, 
				double[] sub_data, double[] plant_data2, double[][] theta, double[][] DthetaDt, double[] f, double[][] pe, double[][] red, double[] F,
				double[] Ft, double[] Peffect, double[] AWI, double[] esc, double irr_vol, double Ptot_cum_event)
		{
	//INPUT
	//P                                      //Precipitation (mm)
	//runon                                  //Runon volume coming from other subareas (mm)
	//drainage                               //Drainage connection between this subarea and others (//)
	//ET0                                    //Evapotranspiration (mm/h)
	//R                                      //Irrigation (mm)
	//time_inf                               //Temporal information
	//input_sub                              //subarea parameters
	//input_lay                              //All layers parameters
	//Veg_data                             //Kcrop information for the vegetation involved


	//Temporal information
	double dt=general_inf[Constants.TWO]/60.;                         //Modeling time interval (h)
	int t_s=general_inf[Constants.THREE];                        //Time interval between events (h)

	//Parameters for estimating hydrological processes

	double A=outflow_data[Constants.THREE];                         //subarea area (m2)
	double D=outflow_data[Constants.FIVE];                         //Surface storage depth (mm)                                       
	                     
	int nL=general_inf[Constants.FOUR];                       //Number of layers in the subarea

	double vegcov=plant_data2[Constants.ONE]/100.;          //Vegetation coverage
	double S=plant_data2[Constants.THREE];                   //Maximum water stored (mm)
	double k=plant_data2[Constants.TWO];                   //Leaf area index


	//TODO, I don't think layer is used again
	//Layer's information
	//Information for each layer is ordered in 'layer' matrix. Each columnn correspond to one layer
//	double[][] layer=new double[10][nL];
//	for (int i=0;i<10;i++)
//	{
//		for (int j=0;j<nL;j++)
//		{
//			layer[i][j] = Double.NaN;
//		}
//	}
//	for (int j=Constants.ONE;j<nL;j++)
//	{
//		layer[:][j]=sub_data(1:10);
//	}
	//The parameters below are defined for all layers
//	double[] d=layer(1,:)*1000/nL;                      //Layer thickness (mm)
//	double[] theta_i=layer(2,:);                     //Initial soil water content  (m3/m3)
//	double[] theta_r=layer(3,:);                     //Residual soil water content (m3/m3)
//	double[] theta_s=layer(4,:);                     //Saturated soil water content (m3/m3)
//	double[] n=layer(5,:);                           //Curve shape parameter
//	double[] m=1-1./n;                               //Second curve shape parameter
//	double[] L=layer(6,:);                           //Empirical pore tortuosity
//	double[] Ks=layer(7,:);                          //Saturated hydraulic conductivity (mm/h)
//	double[] psi_b=layer(8,:);                       //Bubbling pressure (mm)
//	//WP=layer(9,:);                          //Wilting point (m3/m3)
//	double[] FC=layer(10,:);                         //Field capacity (m3/m3)
	
	double[] d=GreenRoofCommon.fillArray(nL, sub_data[Constants.ONE]*1000./nL);                      //Layer thickness (mm)
	double[] theta_i=GreenRoofCommon.fillArray(nL, sub_data[Constants.TWO]);                     //Initial soil water content  (m3/m3)
	double[] theta_r=GreenRoofCommon.fillArray(nL, sub_data[Constants.THREE]);                     //Residual soil water content (m3/m3)
	double[] theta_s=GreenRoofCommon.fillArray(nL, sub_data[Constants.FOUR]);                     //Saturated soil water content (m3/m3)
	double[] n=GreenRoofCommon.fillArray(nL, sub_data[Constants.FIVE]);                           //Curve shape parameter
	double[] m=GreenRoofCommon.fillArray(nL,1.-1./sub_data[Constants.FIVE]);                               //Second curve shape parameter
	double[] L=GreenRoofCommon.fillArray(nL, sub_data[Constants.SIX]);                           //Empirical pore tortuosity
	double[] Ks=GreenRoofCommon.fillArray(nL, sub_data[Constants.SEVEN]);                          //Saturated hydraulic conductivity (mm/h)
	double[] psi_b=GreenRoofCommon.fillArray(nL, sub_data[Constants.EIGHT]);                       //Bubbling pressure (mm)
	//WP=layer(9,:);                          //Wilting point (m3/m3)
	double[] FC=GreenRoofCommon.fillArray(nL, sub_data[Constants.TEN]);                         //Field capacity (m3/m3)
	
	
	double theta_e=theta_s[Constants.ONE]-theta_r[Constants.ONE];          //Effective soil water content (m3/m3)


	//OUTPUT
	//runoff                                 //Runoff volume that drains from this subarea (m3)
	//flp                                    //Summary of the firts layer processes
	    //Peffect                            //Precipitation plus runoff minus interception (mm)
	    //F                                  //Cumulative infiltration (mm)
	    //esc                                //effective surface runoff (mm)
	    //f                                  //Infiltration rate (mm/h)
	    //E_T                                //Evaporation or Evapotranspiration (mm/h)
	//theta                                  //Soil water content (m3/m3)
	//pe                                     //Percolation (mm/h)
	//red                                    //Redistribution (mm/h)
	//Q_sur                                  //Flow of surface runoff(m3/h)
	//Q_sub                                  //Flow of subsurface runoff (m3/h)
	//irr_vol                                //Irrigation volume (m3)

	//VARIABLES'S INITIATION
	double[] e=new double[P.length];                   //Variable assistant to calculate surface runoff (mm)
	double[] theta_p=new double[nL];                    //Variable assistant to calculate soil water content (m3/m3)
	double[] K_p=new double[nL];                        //Variable assistant to calculate unsaturated hydraulic conductivity (mm/h)
	double[] psi_p=new double[nL];                      //Variable assistant to calculate suction head at wetting front (mm)

	double[][] Aeq=new double[nL][2+2*nL];                   //Matix assistant to calculate mass balance

	for (int j=Constants.ONE;j<nL;j++)
	{
	    if (j==Constants.ONE) //For the first layer
	    {
	        Aeq[Constants.ONE][Constants.ONE]=dt/d[Constants.ONE];
	        Aeq[Constants.ONE][Constants.TWO]=-dt/d[Constants.ONE];
	        Aeq[Constants.ONE][Constants.THREE]=-dt/d[Constants.ONE];
	        Aeq[Constants.ONE][Constants.THREE+nL]=-1.;    
	    }
	    else //For the deeper layers
	    {
	        Aeq[j][j+2]=-dt/d[j];
	        Aeq[j][j+1]=dt/d[j];
	        Aeq[j][2+nL+j]=-1.;
	    }  
	}

	//TOTAL PRECIPITATION	
//	double[] Ptotal=P+runon/A*1000; //(mm)
	double[] Ptotal=GreenRoofCommon.add(GreenRoofCommon.multiply( (1.0/A*1000.), runon), P);
	
	//Ptotal includes rainfall plus runoff from upstream subarea
	//Ptotal is uniformly distribuited


	//INITIAL CONDITIONS
	for (int count=0;count<theta_i.length;count++)
	{
//		theta[Constants.ONE][:]=theta_i(:);                  //Initial soil water content (m3/m3)
		theta[Constants.ONE][count]=theta_i[count];
	}
	
	int i=tstep;

	//IRRIGATION
	double Ri=R[i]; 
	irr_vol=irr_vol+Ri*A/1000.; //Irrigation volume (m3)
	    
	//INTERCEPTION (mm) 
	    //Total precipitation minus interception
	        Peffect[i]=Ptotal[i]+Ri;
	        
	    //What are the new precipitation events?
	        //time_events indicates the amount of dt from the beginning of the events
//	        [time_events]=p_events(Peffect(1:i),t_s,dt);
	        PEvents pev = new PEvents();
	        int[] time_events=pev.p_events_(GreenRoofCommon.subsetArray(Peffect, 0, i),t_s,dt);
	        //Interception
	        double[] ET0_p;
	        int[] time_events_p;
	        if (i==Constants.ONE) //At the beginning
	        {
	            ET0_p=new double[]{E_T[i]};
	            time_events_p=new int[]{time_events[i]};
	            Ptot_cum_event=Peffect[i];
	        }
	        else
	        {
	            ET0_p=new double[]{E_T[i], E_T[i-1]};
	            time_events_p=new int[]{time_events[i],time_events[i-1]};    
	        } 
	        
//	        [Peffect[i],Ptot_cum_event]=
	        Interception interception = new Interception();
	        TreeMap<String,Double> interceptionReturn = 
	        		interception.Interception_(Peffect[i],ET0_p,Ptot_cum_event,time_events_p,k,S*vegcov,dt);
	        Peffect[i]= interceptionReturn.get("Peffect");
	        Ptot_cum_event= interceptionReturn.get("Ptot_cum_event");
	           
	       
	//AVAILABLE WATER TO INFILTRATE (mm)
	//If there is a reservoir, AWI includes the amount of water stored in the previous period
	    AWI[i]=AWI[i]+Peffect[i];
	    
	 double Se=0,K,t;  
	//PERCOLATION (mm/h) - For all layers    
	    for (int j=Constants.ONE;j<nL;j++)
	    {
	        //van Genuchten parameters
	        Se=(theta[i][j]-theta_r[j])/(theta_s[j]-theta_r[j]);     //Relative saturation(mm3/mm3)
	        K=(Ks[j]*Math.pow(Se,L[j]))* Math.pow((1.- Math.pow((1.-Math.pow(Se,(1./m[j]))),(m[j])) ),2);         //Unsaturated hydraulic conductivity (mm/h)
	        t=Math.max((theta[i][j]-FC[j])*d[j]/K,0);                     //Travel time through the layer(h)                        
	        //Percolation (mm/h)
	        if (theta[i][j]>FC[j])
	        {
	            pe[i][j]=((theta[i][j]-FC[j])*(1.-Math.exp(-dt/t))*d[j]/dt); //se multipica por d/dt para tener consistencia en las unidades
	        }
	        else
	        {
	            pe[i][j]=0;
	        }
	    }

	    
	    
	//INFILTRATION (mm/h) - for the first layer
	    //Updating new precipitation events
	    time_events=pev.p_events_(GreenRoofCommon.subsetArray(AWI, 0, i),t_s,dt);
	    double F00=0;
	    double F0;
	    //Initial flooding for Green Ampt calculation
	    if (i==Constants.ONE)
	    {
	         F0=F00;
	    }
	    //No flooding is assumed for each new event
	    else if (time_events[i]==1)
	    {
	         F0=0;
	    }
	    else
	    {
	    //Initial flooding is equal to cumulative infiltration from previous time interval
	         F0=Ft[i-1];
	    }
	    //Suction head at wetting front (mm)
	    double psi=psi_b[Constants.ONE]*Math.pow((Math.pow(((theta[i][Constants.ONE]-theta_r[Constants.ONE])/(theta_s[Constants.ONE]-theta_r[Constants.ONE])),(-1./m[Constants.ONE]))-1.),(1./n[Constants.ONE]));
	    //Infiltration (mm/h)
	    GreenAmpt greenAmpt = new GreenAmpt();
//	    [f[i],Ft[i]]=
	    TreeMap<String,Double> greenAmptReturnValues = greenAmpt.Green_Ampt_(F0,psi,Ks[Constants.ONE],dt,Se,theta_e,AWI[i]);
	    f[i]=greenAmptReturnValues.get("f");
	    Ft[i]=greenAmptReturnValues.get("Ft");
	    		
	    
	//MASS BALANCE - For all layers
	    double limit_sup=0; //Auxiliary variable
	    for (int j=Constants.ONE;j<nL;j++)
	    {
	        //Soil water content considering the above hydrological processes
	        if (j==Constants.ONE) //First layer
	        {
	            theta_p[j]=theta[i][j]+(f[i]-E_T[i]-pe[i][j])*dt/d[j];
	        }
	        else //Deeper layers
	        {
	            theta_p[j]=theta[i][j]+(pe[i][j-1]-pe[i][j])*dt/d[j];
	        }
	        //If the soil water content of the first layer exceeds theta_s, then
	        //the infiltration rate is decreased
	        if (theta_p[Constants.ONE]>theta_s[Constants.ONE])
	        {
	            Ft[i]=Ft[i]-(f[i]-((theta_s[Constants.ONE]-theta[i][Constants.ONE])*d[Constants.ONE]/dt+E_T[i]+pe[i][Constants.ONE]))*dt;
	            f[i]=((theta_s[Constants.ONE]-theta[i][Constants.ONE])*d[Constants.ONE]/dt+E_T[i]+pe[i][Constants.ONE]);
	        }
	        //Try again
	        if (j==Constants.ONE) //First layer
	        {
	            theta_p[j]=theta[i][j]+(f[i]-E_T[i]-pe[i][j])*dt/d[j];
	        }
	        else //Deeper layers
	        {
	            theta_p[j]=theta[i][j]+(pe[i][j-1]-pe[i][j])*dt/d[j];
	        }
	        //Analysis of the limits of soil water content: 
	        //Register if theta of some of the layers exceeds theta_s or is less than theta_r
	        if (theta_p[j]<theta_r[j] || theta_p[j]>theta_s[j])
	        {
	            limit_sup=1;
	        }
	    }
	    //If theta is maintained within the limits, its value is updated
	    if (limit_sup==0 )
	    {
//	        theta(i+1,:)=theta_p;
	        for (int count = 0;count<theta_p.length;count++)
	        {
	        	theta[i+1][count]=theta_p[count];
	        }
	    }
	    //If theta exceeds a limit, hydrological rates are adjusted as optimization problem     
	    //TODO not sure this gets called
//	    else
//	    {
//	    	//TODO fix this whole section
//	    	
//	        //Initial values
//	        x0=new double[]{f[i],E_T[i],pe(i,:),theta_p};   //x=[f E_T pe theta]
//	        //Limits
//	        li=[0,0,zeros(1,nL),theta_r*1.01]; //Lower limit
//	        ls=[f[i],E_T[i],pe(i,:),theta_s*0.99]; //Upper limit
//	        options=optimset('Algorithm','active-set');
//	        //The problem is solved by maximizing the sum of rates f+E_T+pe
//	        [x]=fmincon(@(x)(-sum(x(1:2+nL))),x0,[],[],Aeq,-theta(i,:),li,ls,[],options);
//	        //Rates are updated
//	        //If theta exceeds theta_s, the infiltration rate must be reduced
//	        if (x[Constants.ONE]<f[i])
//	        {
//	            Ft[i]=Ft[i]-(f[i]-x[Constants.ONE])*dt;
//	            f[i]=x[Constants.ONE];
//	        }
//	        //If theta is less than theta_r, the evapotranspiration or
//	        //percolation rate must be reduced
//	        //E_T[i]=x(2);
//	        pe(i,:)=x(3:2+nL);
//	        //The rate of change of soil water content is updated
//	        for (int j=Constants.ONE;j<nL;j++)
//	        {
//	            if (j==Constants.ONE) //First layer
//	            {
//	                theta[i+1][j]=theta[i][j]+(f[i]-E_T[i]-pe[i][j])*dt/d[j];
//	            }
//	            else //Deeper layers
//	            {
//	                theta[i+1][j]=theta[i][j]+(pe[i][j-1]-pe[i][j])*dt/d[j];
//	            }
//	        }
//	    }
	    
	    
	//CUMULATIVE INFILTRATION (mm) - For the first layer
	//Infiltration depth during the time step dt. It value is calculated based on
	//cumulative infiltration for each precipitation events
	    double dF;
	    if (i==Constants.ONE)
	    {
	        dF=Ft[i]-F00;
	    }
	    else if (time_events[i]==Constants.ONE) //At the beginning of the event
	    {
	        dF=Ft[i];
	    }
	    else //During the event
	    {
	        dF=Ft[i]-Ft[i-1];
	    }
	    //Cumulative infiltration considering the complete simulation
	    if (i!=Constants.ONE)
	    {
	        F[i]=(F[i-1]+dF);
	    }
	    else
	    {
	        F[i]=Ft[i];
	    }
	    
	    
	//SURFACE RUNOFF (mm) - For the first layer
	    //Infiltration depth during the time step dt. It value is calculated based on
	    //cumulative infiltration considering the complete simulation
	    double Fp;
	    if (i==Constants.ONE) //At the beginning of the simulation
	    {
	        Fp=F[i];
	    }
	    else //During the simulation
	    {
	        Fp=F[i]-F[i-1];
	    }
	    //The rainfall excess (mm) is calculated in each time as the difference
	    //between AWI and the infiltration depth
	    e[i]=Math.max(AWI[i]-Fp,0);
	    //If the subarea no considers a surface storage, the runoff is
	    //equal to rainfall excess
	    if (D==0)
	    {
	        esc[i]=e[i];
	    //If the subarea considers a surface storage, then surface runoff
	    //is generated only once the storage capacity is full
	    }
	    else if (D>0) 
	    {
	       if (e[i]>D) //If storage is full
	       {
	           esc[i]=e[i]-D;
	           AWI[i+1]=D;
	       }
	       else //If storage is not full, no runoff is generated
	       {
	           AWI[i+1]=e[i];
	           esc[i]=0;
	       }
	    }
	    
	//REDISTRIBUTION (mm/h) - for all layers
	    //Redistribution occurs when the infiltration is null and the
	    //subarea is composed of more than one layer
	    if (f[i]<1e-12 && nL>1)
	    {
//	        rr=ones(1,nL-1); //Redistribution factor
	        double[] rr = GreenRoofCommon.fillArray(nL-1, 1);
	        for (int j=Constants.ONE;j<nL;j++)
	        {
	            //van Genuchten parameters considering the theta values after
	            //considering the above hydrological processes
	            Se=(theta[i+1][j]-theta_r[j])/(theta_s[j]-theta_r[j]);	            
	            psi_p[j]=(psi_b[j]*Math.pow((Math.pow(Se,(Math.pow(-m[j],(-1))))-1),(Math.pow(n[j],(-1)))));	            
	            K_p[j]=(Ks[j]*Math.pow(Se,L[j])*Math.pow((1-Math.pow((1-Math.pow(Se,(n[j]/(n[j]-1)))),((n[j]-1)/n[j]))),2));
	        }
	        for (int j=Constants.TWO;j<nL-1;j++)
	        {
	            //If the soil is composed of 3 layers or more and the flow from the middle layers is
	            //established both upward and downward, the total redistribution
	            //flow obtained is split according to a factor 'rr'
	            if (psi_p[j-1]>psi_p[j] && psi_p[j+1]>psi_p[j])
	            {
	                rr[j-1]=Math.abs(psi_p[j]-psi_p[j-1])/(Math.abs(psi_p[j]-psi_p[j-1])+Math.abs(psi_p[j]-psi_p[j+1]));
	                rr[j]=Math.abs(psi_p[j]-psi_p[j+1])/(Math.abs(psi_p[j]-psi_p[j-1])+Math.abs(psi_p[j]-psi_p[j+1]));
	            }
	        }
	        for (int j=Constants.ONE;j<nL-1;j++)
	        {
	            if (psi_p[j]<=psi_p[j+1]) //The flow downward
	            {
	                red[i][j]=K_p[j]*((psi_p[j]-psi_p[j+1])/(0.5*(d[j]+d[j+1]))-1)*rr[j];//(negative)
	            }
	            else //The flow upward
	            {
	                red[i][j]=K_p[j+1]*((psi_p[j]-psi_p[j+1])/(0.5*(d[j]+d[j+1]))-1)*rr[j];//(positive)
	            }
	        }
	        limit_sup=0; //Auxiliary variable
	        for (int j=Constants.ONE;j<nL;j++)
	        {
	            if (j==Constants.ONE) //For the firt layer
	            {
	                theta_p[j]=theta[i+1][j]+red[i][j]*dt/d[j];
	            }
	            else if (j>Constants.ONE && j<nL) //For the middle layers
	            {
	                theta_p[j]=theta[i+1][j]+(red[i][j]-red[i][j-1])*dt/d[j];
	            }
	            else if (j==nL) //For the last layer
	            {
	                theta_p[j]=theta[i+1][j]-red[i][j-1]*dt/d[j];
	            }
	            //Analysis of the limits of soil water content: 
	            //Register if theta of some of the layers exceeds theta_s or is less than theta_r
	            if (theta_p[j]<theta_r[j] || theta_p[j]>theta_s[j])
	            {
	                limit_sup=1;
	            }
	        }
	        //If theta exceeds a limit, redistribution rate is adjusted as
	        //optimization problem  
	        //TODO, not sure if this is called
//	        if (limit_sup==1  )
//	        {
//	            dir_red=sign(red(i,1:nL-1)); //Redistribution flow orientation
//	            //Matrix Aer
//	            Aer=zeros(nL,2*nL-1);
//	            for (int j=Constants.ONE;j<nL;j++)
//	            {
//	                if (j==Constants.ONE) //For the first layer
//	                {
//	                    Aer[j][j]=1;
//	                    Aer[j][nL+j]=-dir_red[j]*dt/d[j]; 
//	                }
//	                else if (j>Constants.ONE && j<nL) //For the middle layers
//	                {
//	                    Aer[j][j]=1;
//	                    Aer[j][nL+j-1]=dir_red[j-1]*dt/d[j];
//	                    Aer[j][nL+j]=-dir_red[j]*dt/d[j];
//	                }
//	                else //For the last layer
//	                {
//	                    Aer[j][j]=1;
//	                    Aer[j][nL+j-1]=dir_red[j-1]*dt/d[j];
//	                }
//	            }
//	            //Initial values
//	            y0=[theta_p,abs(red(i,1:nL-1))];   //y=[theta(nL) red(nL-1)]
//	            //Limits
//	            li=[theta_r*1.01,zeros(1,nL-1)]; //Lower limit
//	            ls=[theta_s*0.99,abs(red(i,1:nL-1))]; //Upper limit
//	            options=optimset('Algorithm','active-set');
//	            //The problem is solved by maximizing the sum of rates red
//	            [y]=fmincon(@(y)(-sum(y(nL+1:end).*(1-dir_red)/2)),y0,[],[],Aer,theta(i+1,:),li,ls,[],options);
//	            //Rates are updated
//	            red(i,1:nL-1)=dir_red.*y(nL+1:end);               
//	        }    
	    }
	  
	//RATE OF CHANGE OF SOIL WATER CONTENT
	    for (int j=Constants.ONE;j<nL;j++)
	    {
	        if (j==Constants.ONE)
	        {
	            DthetaDt[i][j]=f[i]-E_T[i]-pe[i][j]+red[i][Constants.ONE];
	        }
	        else
	        {
	            DthetaDt[i][j]=pe[i][j-1]-pe[i][j]-red[i][j-1]+red[i][j];
	        }
	        theta[i+1][j]=theta[i][j]+DthetaDt[i][j]*dt/d[j];
	    }
	    
//	return items    theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event
	    TreeMap returnValues = new TreeMap();
	    returnValues.put("theta", theta);
	    returnValues.put("DthetaDt", DthetaDt);
	    returnValues.put("f", f);
	    returnValues.put("pe", pe);
	    returnValues.put("red", red);
	    returnValues.put("F", F);
	    returnValues.put("Ft", Ft);
	    returnValues.put("Peffect", Peffect);
	    returnValues.put("AWI", AWI);
	    returnValues.put("esc", esc);
	    returnValues.put("irr_vol", irr_vol);
	    returnValues.put("Ptot_cum_event", Ptot_cum_event);
	
	    
	    return returnValues;

		}
		
}
