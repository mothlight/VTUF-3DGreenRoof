package greenroof;

import java.util.ArrayList;

public class MainMix
{
	
    int meteo_dt;
    int sim_dt ;
    int event_dt;
    int n_layers;
    int max_steps;    
    double Area;
    double roof_height;
    boolean use_epw = false;
    ArrayList<double[]> weatherdata;
    int nL;
    String input_file;
    String epw_file;
    int[] general_inf;
    double[] outflow_data ;
    double surface_sensor_depth = 0.01;
    
	public void readData()
	{
		
	    double SiteWindBLHeight;
	    double SiteWindExp;    
	// 


	//%%             LLENAR ESTO             %%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    //input_file = "INPUTSChicago_short.xlsx";
	    //input_file = "INPUTSLIVE_v03.xlsx";
	    //input_file = "INPUTSLIVE_v07.xlsx";
	    String input_file = "INPUTSLIVE_Cv10.xlsx";
	    //input_file = "INPUTSChicago_revgf.xlsx";

	    
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//% PARSE GENERAL DATA

	// Read
	//[numbers,files] = xlsread(input_file,'General','K4:K14');

	// Parse Files
//	out_file = files(1);
//	if(files(2) ~= "")    
//	    use_epw = true;    
//	    epw_file = files(2);
//	    disp("EPW file '" +epw_file+"' will be used");
//	else
//	    use_epw = false;
//	end
//	model_to_use = files(4);
//	location = files(7);
//	clear files
	    



//	// parse numbers
//	meteo_dt = numbers(1); //minutes
//	sim_dt = numbers(2); //minutes
//	event_dt = numbers(3); //hours
//	n_layers = numbers(4);
//	max_steps = numbers(7);
//	Area = numbers(9);
//	roof_height = numbers(10);
	
    meteo_dt= 5;
    sim_dt = 1;
    event_dt =6;
    n_layers = 5;
    max_steps = 3865;
    
    Area = 10;
    roof_height = 3;

	// These values correspond to OCEAN
    String location = "City";

    
	if (location.equals("Country"))
	{
	    SiteWindBLHeight = 270;
	    SiteWindExp = 0.14;    
	}
	else if (location.equals("Suburbs"))
	{
	    SiteWindBLHeight = 370;
	    SiteWindExp = 0.22;
	}
	else if (location.equals("City"))
	{
	    SiteWindBLHeight = 460;
	    SiteWindExp = 0.33;    
	}
	else if (location.equals("Ocean"    ))
	{
	    SiteWindBLHeight = 210;
	    SiteWindExp = 0.1;    
	}
	else if (location.equals("Urban"))
	{
	    SiteWindBLHeight = 370;
	    SiteWindExp = 0.22;   
	}
	else
	{
	    System.out.println("Not recognized location");
	    SiteWindBLHeight = 0;
	    SiteWindExp = 0.;  
	}

	// GF
	int[] general_inf=new int[]{meteo_dt,sim_dt,event_dt,n_layers};
	nL=n_layers; 
//	clear numbers

	//GF
	double Dt=meteo_dt/60.; // Dt expressed in hours                      
	double dt=sim_dt/60.; // Dt expressed in hours  

	//// PARSE OUTFLOW
	//GF
//	outflow_data=xlsread(input_file,'Outflow','K4:K11');
	outflow_data = new double[]{0.01,
			1,
			1,
			0.01,
			0,
			5.56E-03,
			5.56E-03,
			1
		};

	//TODO don't use indoor for now
	//// PARSE OUTDOOR 
//	interior_temperature = xlsread(input_file,"TODO","A:A")+273.15;

	
	//TODO for now, just hard code a bit of data
//	if(use_epw)
//	{
//	    // Load data from EPW
//	    weatherdata = csvread(char(epw_file),8,6);
//	    weatherdata = [weatherdata(end,:);weatherdata];
//	}
//	else
//	{
//	    // Load data from excel
//	    weatherdata = xlsread(input_file,"Outdoors","F:K");
//	}
	
//	OUTDOORS INFORMATION										
//	Year 	Month	Day	Hour	Minute	Precipitation	Solar radiation 	External temperature	External relative humidity	Wind speed at 2m	Irrigation
//						[mm]	[W/m²]	[°C]	[%]	[m/s]	[mm]
//	YYYY	MM	DD	HH	MM	P	Sol rad	T	RH	u	Irrig 
//	2017	12	16	0	0	0.0	-2.192	19.6	53.65	0.1	0
//	2017	12	16	0	5	0.0	-2.460	19.9	54.01	0.3	0
//	2017	12	16	0	10	0.0	-2.412	19.9	54.01	0.2	0
//	2017	12	16	0	15	0.0	-2.676	19.8	53.70	0.3	0
	weatherdata = new ArrayList<double[]>();
	weatherdata.add(new double[]{2017,	12,	16,	0,	0,	0.0,	-2.192,	19.6,	53.65,	0.1,	0});
	weatherdata.add(new double[]{2017,	12,	16,	0,	5,	0.0,	-2.460,	19.9,	54.01,	0.3,	0});
	weatherdata.add(new double[]{2017,	12,	16,	0,	10,	0.0,	-2.412,	19.9,	54.01,	0.2,	0});
	weatherdata.add(new double[]{2017,	12,	16,	0,	15,	0.0,	-2.676,	19.8,	53.70,	0.3,	0});
	}

	public void run()
	{
	//// CREATE PLANT
	Plant plant = new Plant();

//	plant_data = xlsread(input_file,'Plants','E5:E14');
	double[] plant_data = new double[]{1,
			1,
			0.2,
			0.9,
			0.5,
			0.1,
			0.4,
			0.4,
			200,
			0.001
		};

	plant.LAI=plant_data[Constants.ONE];        //leaf area index A FEB Tabares  
	plant.rho=plant_data[Constants.THREE];     //sw reflectivity plant
	plant.em=plant_data[Constants.FOUR];    //emissivity planta
	plant.k=plant_data[Constants.FIVE];      // thermal conductivity planta
	plant.height=plant_data[Constants.SIX];   //  height of plant [m]
	plant.ks=plant_data[Constants.SEVEN]; // Extinsion coefficient   
	plant.ks_ir = plant_data[Constants.EIGHT]; // IR extinsion coefficient
	plant.rsmin=plant_data[Constants.NINE];      //[s/m] resistencia estomatica 
	plant.z=2;       // [m] height of wind measurements
	plant.Zog = plant_data[Constants.TEN]; // VARIA SEGUN SMOOOTH
//	clear plant_data

	//GF
//	plant_data2=xlsread(input_file,'Plants','E4:E14');
	double[] plant_data2 = new double[]{100,1,
			1,
			0.2,
			0.9,
			0.5,
			0.1,
			0.4,
			0.4,
			200,
			0.001
		};

	//// CREATE SUBSTRATE
	Substrate substrate = new Substrate();
//	sub_data = xlsread(input_file,'Substrate','E4:E18');
	double[] sub_data = new double[]{0.15,
			0.48,
			0.01,
			0.7,
			1.50,
			0.5,
			782,
			169.0,
			0.150,
			0.550,
			0.18,
			0.0148,
			950,
			1400,
			0.9
		};

	substrate.depth=sub_data[Constants.ONE];     // depth substrate  (ex t)
	double initial_VWC = sub_data[Constants.TWO]; // This is used later.
	substrate.VWCresidual = sub_data[Constants.THREE];
	substrate.VWCsat=sub_data[Constants.FOUR];     // VWC saturated
	substrate.VWCwilting=sub_data[Constants.NINE];    // VWC wilting point
	substrate.VWCfc = sub_data[Constants.TEN];
	substrate.rho=sub_data[Constants.ELEVEN];     //sw reflectivity substrate
	substrate.k=sub_data[Constants.TWELVE]; //W/mK;
	substrate.dens =sub_data[Constants.THIRTEEN]; //kg/m^3
	substrate.Cp=sub_data[Constants.FOURTEEN]; // J/kgK;
	substrate.em = sub_data[Constants.FIFTEEN];
	// substrate.phi=0.85; 
	// substrate.Lsubm=0.05;
	//clear sub_data
	//GF : GM--> Initial VWC is given in cell E5
	//GM: De acuerdo... modificado 
	//substrate.VWC = 0.23; // initial Volumetric Water Content


	//// CREATE ROOF

	Roof roof = new Roof();
//	roof_data = xlsread(input_file,'Support','E4:E7');
	double[] roof_data = new double[]{0.15,
			1.63,
			2400,
			920
			};

	roof.depth = roof_data[Constants.ONE]; // meters
	roof.k = roof_data[Constants.TWO]; //thermal conductivity
	roof.density = roof_data[Constants.THREE]; 
	roof.Cp = roof_data[Constants.FOUR]; // heat capacity

//	clear roof_data

	String model_to_use = "Tabares";
	//// CREATE MODEL
	System.out.println("Using model "+model_to_use);
	TabaresThermalMass model;
	if (model_to_use.equals("Tabares"))
	{
	    model = new TabaresThermalMass(plant,substrate,roof,n_layers,n_layers,Area,sim_dt);
	}
//	else if (model_to_use.equals("Sailor"))
//	{
//	    model = SailorThermalMass_(plant,substrate,roof,n_layers,n_layers,Area,sim_dt);
//	}
	else    
	{
		System.out.println("Unknown model to use");
	}
	
//	model.VWC = initial_VWC*ones(n_layers,1);
	model.VWC = new double[n_layers];
	for (int i=0;i<n_layers;i++)
	{
		model.VWC[i]=initial_VWC;
	}

	//// FIX TIMESTEPS
//	if (use_epw)
//	{
//	    meteo_dt = 60; //EPW has data every 60 minutes
//	}
	int n_sub_tsteps = meteo_dt/sim_dt;
	int total_steps = n_sub_tsteps*(max_steps-1);
	double[] P = new double[weatherdata.size()];
	double[] R = new double[weatherdata.size()];

	//// INTERPOLATION FOR PRECIPITATION P AND IRRIGATION R
	//GF
	for (int i=Constants.ONE;i<weatherdata.size();i++)
	{
//	    if(use_epw)
//	    {
//	        P(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,28)*dt/Dt;
//	        R(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,6)*0;
//	    }
//	    else
//	    {
		double[] weatherline = weatherdata.get(i);
	    P(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherline[Constants.ONE]*dt/Dt;
	    R(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherline[Constants.SIX]*dt/Dt;
//	    }
	}

	//// CREATE RESULTS VECTORS
	double[] result_T_sky = new double[total_steps];
	double[] result_T_out = new double[total_steps];
	double[] result_wind_speed = new double[total_steps];
	double[] result_Rsh = new double[total_steps];
	double[] result_T_plants = new double[total_steps];
	double[] result_T_substrate = new double[total_steps];
	double[] result_T_5cm = new double[total_steps]; //5cm
	double[] result_T_10cm = new double[total_steps]; //10cm
	double[] result_T_15cm = new double[total_steps]; //15cm
	double[] result_T_interior = new double[total_steps];
	double[] result_VWC_surface = new double[total_steps];
	double[] result_VWC_mid = new double[total_steps];
	double[] result_interface_heat_flux = new double[total_steps];
	double[] result_evaporation = new double[total_steps];
	double[] result_transpiration = new double[total_steps];
	double[] result_substrate_convection = new double[total_steps];
	double[] result_plant_convection = new double[total_steps];
	double[] result_Rain = new double[total_steps];
	double[] result_ET = new double[total_steps];
	double[] result_heating_load = new double[total_steps];

	// Plants balance
	double[] result_plant_absorbed_solar = new double[total_steps];
	double[] result_plant_absorbed_ir_sky = new double[total_steps];      

	// Both
	double[] result_Qir = new double[total_steps];

	// Substrate balance
	double[] result_substrate_solar_radiation = new double[total_steps];
	double[] result_substrate_infrared_radiation = new double[total_steps];
	double[] result_substrate_conduction = new double[total_steps];

	//// INITIALISATION OF VARIABLES
	//GF
	double[][] runon=new double[P.length][1];               //Runoff volume entering to the subcatchment (m3)
	double[][] runoff=new double[P.length][1];              //Runoff volume which drains from the subcatchment (m3)
	double[][] Q_out=new double[(int)Math.round(Math.floor(P.length*1.1))][1];      //Outflow from the subcatchment to the street (m3/h) 

	double[][] theta=new double[P.length][nL];              //Soil water content (m3/m3)
	double[][] DthetaDt=new double[P.length][nL];           //Rate of change in soil water content (mm/h)
	double[][] f=new double[P.length][1];                   //Infiltration (mm/h)
	double[][] pe=new double[P.length][nL];                 //Percolation (mm/h)

	double[][] red=new double[P.length][nL];                //Redistribution (mm/h)
	double[][] F=new double[P.length][1];                   //Cumulative infiltration (mm)
	double[][] Ft=new double[P.length][1];                  //Cumulative infiltration to calculate Green Ampt (mm)
	double[][] Peffect=new double[P.length][1];             //Precipitation plus runoff minus interception (mm)
	double[][] AWI=new double[P.length][1];                 //Available water to infiltrate (mm)
	double[][] esc=new double[P.length][1];                 //effective surface runoff (mm)

	//// INITIAL VALUES
	//GF
	double irr_vol=0;                              //Irrigated volume (m3)
	double Ptot_cum_event=0;

	//// SIMULATE
	String sim_desc;
//	if (use_epw)
//	{
//	    sim_desc = "Performing simulation with inputs from "+input_file+" using "+model_to_use+" model and data from "+epw_file;
//	}
//	else    
//	{
	    sim_desc = "Performing simulation with inputs from "+input_file+" using "+model_to_use+" model and custom data";
//	}
//	h = waitbar(0,char(sim_desc));


	for (int main_step = Constants.ONE;main_step<weatherdata.size();main_step++)
	{
	    
	    double[] this_data_line = weatherdata.get(main_step);
	    double[] next_data_line = weatherdata.get(main_step+1);
	    
	    for (int k=Constants.ONE;k<n_sub_tsteps;k++)
	    {
	        
	        int tstep = (main_step-1)*n_sub_tsteps+k;
//	        waitbar( tstep/ (max_steps*n_sub_tsteps));
	            
	        //TODO getting rid of indoor for now
	        // interpolate   
//	        this_inner_t = interior_temperature(main_step);
//	        next_inner_t = interior_temperature(main_step+1);
//	        inner_T = this_inner_t + (k-1)*(next_inner_t - this_inner_t)/n_sub_tsteps;        
//	        model.T_interior = inner_T;
	        
	        
	                
	        
//	        if(use_epw)
//	        {
//	            data_line = EPWWeatherDataLine(this_data_line + (k-1)*(next_data_line - this_data_line)/n_sub_tsteps,roof_height,SiteWindBLHeight,SiteWindExp);
//	        }
//	        else
//	        {
	        //     	   rainfall,R_sh,Tair,RH,U,irrigation
	        //TODO, I think think this is interpolating, but just use the single line for now
//	            double[] data_line = WeatherDataLine(this_data_line + (k-1)*(next_data_line - this_data_line)/n_sub_tsteps);
	        WeatherDataLine data_line = new WeatherDataLine(this_data_line );
//	        }
	        
	        // Advance
	        model = model.moveForward(data_line);
	                        
	        // GF : to be completed by GM
	        result_ET[tstep] = model.et_mm_hour; // Evaporation + Transpiration in mm/hour
	        //E_T=... //GF expressed in mm/h. Must be a matrix (1 line per timestep)
	    
	        // Perform the mass balance and determine the VWC (Volumetric Water
	        // Content)of each layer
	 
	        MassBalance massBalance = new MassBalance();
//	        [theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event]
//	        returnData = 
	        massBalance.MassBalance_(tstep,P,R,runon, result_ET,general_inf, outflow_data, sub_data,plant_data2,theta,DthetaDt,f,pe,red,F,Ft,Peffect, AWI, esc,irr_vol,Ptot_cum_event);
	       
	        // GF : to be modified by GM 
	        model.VWC= theta(tstep,:);            
	        
	        
	        // Retrieve and store data
	        result_T_sky[tstep] = data_line.Tsky;
	        result_T_out[tstep] = data_line.Tair;
	        result_wind_speed[tstep] = data_line.U;
	        result_Rsh[tstep] = data_line.R_sh;
	        result_Rain[tstep] = data_line.rainfall;
	        result_T_plants[tstep] = model.T_plants;
	        result_T_substrate[tstep] = model.getTemperature(surface_sensor_depth);
	        result_T_5cm[tstep] = model.getTemperature(0.05); //5cm
	        result_T_10cm[tstep] = model.getTemperature(0.1); //10cm
	        result_T_15cm[tstep] = model.getTemperature(0.15); //15cm
	        result_T_interior[tstep] = model.T_interior;
	        result_VWC_surface[tstep] = model.VWC[Constants.ONE];
	        result_VWC_mid[tstep] = model.getVWC(substrate.depth/2);
	        result_interface_heat_flux[tstep] = model.interface_heat_flux;
	        result_evaporation[tstep] = model.evaporation;
	        result_transpiration[tstep] = model.transpiration;
	        result_substrate_convection[tstep] = model.substrate_convection;
	        result_plant_convection[tstep] = model.plant_convection;
	        result_heating_load[tstep] = model.heating_load;
	        
	        result_plant_absorbed_solar[tstep] = model.plant_absorbed_solar ;
	        result_plant_absorbed_ir_sky[tstep] = model.plant_absorbed_ir_sky;      
	        result_Qir[tstep] = model.Qir;
	        result_substrate_solar_radiation[tstep] = model.substrate_solar_radiation;
	        result_substrate_infrared_radiation[tstep] = model.substrate_infrared_radiation;
	        result_substrate_conduction[tstep] = model.substrate_conduction;
	        
	    }
	   

	    
	    // break if needed
	    if (main_step >= (max_steps-1))
	    {
	        break;
	    }
	}
//	close(h);


	headers = [
	    "Sky temperature (�C)",
	    "Exterior temperature (�C)",
	    "Wind speed (m/s)",
	    "Global solar radiation (W/m2)",
	    "Substrate Surface Temperature (�C)",
	    "Temperature 5cm deep (�C)",
	    "Temperature 10cm deep (�C)",
	    "Temperature 15cm deep (�C)",
	    "Foliage Temperature (�C)",
	    "VWC surface",
	    "VWC mid depth",
	    "Substrate sensible heat transfer (W/m2)",
	    "Foliage sensible heat transfer (W/m2)",
	    "Foliage Transpiration (W/m2)",
	    "Substrate evaporation (W/m2)",
	    "Rainfall (mm)",
	    "Interior temperature (�C)",
	    "Heat flux below substrate (W/m2)",
	    "Evapotranspiration (mm/hour)",
	    "plant_absorbed_solar",
	    "plant_absorbed_ir_sky",     
	    "Qir",
	    "substrate_solar_radiation",
	    "substrate_infrared_radiation",
	    "substrate_conduction",
	    "heating_load (W/m2)"
	    ];

	results = [
	    result_T_sky-273.15,
	    result_T_out-273.15,
	    result_wind_speed,
	    result_Rsh,
	    result_T_substrate-273.15,
	    result_T_5cm-273.15,
	    result_T_10cm-273.15,
	    result_T_15cm-273.15,
	    result_T_plants-273.15,
	    result_VWC_surface,
	    result_VWC_mid,
	    result_substrate_convection, 
	    result_plant_convection, 
	    result_transpiration, 
	    result_evaporation, 
	    result_Rain, 
	    result_T_interior-273.15,
	    result_interface_heat_flux,
	    result_ET,    
	    result_plant_absorbed_solar,
	    result_plant_absorbed_ir_sky,      
	    result_Qir,
	    result_substrate_solar_radiation,
	    result_substrate_infrared_radiation,
	    result_substrate_conduction,
	    result_heating_load
	];

	xlswrite(char(out_file),[headers;results(1:n_sub_tsteps:end,:)],char(model_to_use));

	}



}
