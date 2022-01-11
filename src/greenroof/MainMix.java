package greenroof;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

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
    
    
	public static void main(String[] args)
	{
		MainMix main = new MainMix();
		main.run();

	}
    
	public void run()
	{
		
	    double SiteWindBLHeight;
	    double SiteWindExp;    
	// 
	    String dataDir = "/home/kerryn/git/2021-11-VTUF-GreenRoof/Green roof model with IHMORS coupling/";

	//%%             LLENAR ESTO             %%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    //input_file = "INPUTSChicago_short.xlsx";
	    //input_file = "INPUTSLIVE_v03.xlsx";
	    //input_file = "INPUTSLIVE_v07.xlsx";
	    String input_file = "INPUTSLIVE_Cv10.xlsx";
	    //input_file = "INPUTSChicago_revgf.xlsx";

	    
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    FileInputStream fis = null;
	    try
		{
			fis = new FileInputStream(new File(dataDir + input_file));
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
	    
	    XSSFWorkbook wb = null;
	    try
		{
			wb = new XSSFWorkbook(fis);
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	    TreeMap<Integer,String> generalValues = new TreeMap<Integer,String>();
	    XSSFSheet generalSheet = wb.getSheet("General");
	    System.out.println(generalSheet.toString());
	    for (int i=4-1;i<14;i++)
	    {
	    	 XSSFRow row = generalSheet.getRow(i);
	    	 XSSFCell cell = row.getCell(10);//K
	    	 String value = cell.getRawValue();
//	    	 String value = cell.getStringCellValue();
	    	 System.out.println(value);
	    	 generalValues.put(i+1, value);
	    }
	   
	   

	//% PARSE GENERAL DATA

	// Read
	//[numbers,files] = xlsread(input_file,'General','K4:K14');

	// Parse Files
//	out_file = files(1); //TODO
//	if(files(2) ~= "")    //TODO
//	    use_epw = true;    
//	    epw_file = files(2);
//	    disp("EPW file '" +epw_file+"' will be used");
//	else
//	    use_epw = false;
//	end
//	model_to_use = files(4);  //TODO
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
	
    meteo_dt = Integer.parseInt(generalValues.get(4));
    sim_dt = Integer.parseInt(generalValues.get(5));
    event_dt =Integer.parseInt(generalValues.get(6));
    n_layers = Integer.parseInt(generalValues.get(7));
    max_steps = Integer.parseInt(generalValues.get(10));
    
    Area = Integer.parseInt(generalValues.get(12));
    roof_height = Integer.parseInt(generalValues.get(13));

	// These values correspond to OCEAN
    String locationStrInt = generalValues.get(14);
    String location = "City";
    if (locationStrInt.equals("137"))
    {
    	location = "City";
    }

    
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
//	outflow_data = new double[]{0.01,
//			1,
//			1,
//			0.01,
//			0,
//			5.56E-03,
//			5.56E-03,
//			1
//		};
	
	
    TreeMap<Integer,Double> outflowValues = new TreeMap<Integer,Double>();
    XSSFSheet outflowSheet = wb.getSheet("Outflow");
    System.out.println(outflowSheet.toString());
    for (int i=4-1;i<11;i++)
    {
    	 XSSFRow row = outflowSheet.getRow(i);
    	 XSSFCell cell = row.getCell(10);//K
    	 String value = cell.getRawValue();
//    	 String value = cell.getStringCellValue();
    	 System.out.println(value);
    	 outflowValues.put(i+1, Double.parseDouble(value));
    }
    outflow_data = new double[]{outflowValues.get(4),
    		outflowValues.get(5),
    		outflowValues.get(6),
    		outflowValues.get(7),
    		outflowValues.get(8),
    		outflowValues.get(9),
    		outflowValues.get(10),
    		outflowValues.get(11)
		};
    


	//// PARSE OUTDOOR 
//	interior_temperature = xlsread(input_file,"TODO","A:A")+273.15;
//    TreeMap<Integer,Double> interior_temperatureValues = new TreeMap<Integer,Double>();
    double[] interior_temperature = new double[max_steps];
    XSSFSheet interior_temperatureSheet = wb.getSheet("TODO");
    System.out.println(interior_temperatureSheet.toString());
    for (int i=0;i<max_steps;i++)
    {
    	 XSSFRow row = interior_temperatureSheet.getRow(i+1);
    	 XSSFCell cell = row.getCell(0);//K
    	 String value = cell.getRawValue();
//    	 String value = cell.getStringCellValue();
//    	 System.out.println(value);
//    	 outflowValues.put(i+1, Double.parseDouble(value));
    	 interior_temperature[i]=Double.parseDouble(value)+273.15;
    }

	
	//TODO not using EPW yet
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
    
    
//    weatherdata = new double[max_steps][6];
    weatherdata = new ArrayList<double[]>();
    XSSFSheet outdoorsSheet = wb.getSheet("Outdoors");
    System.out.println(outdoorsSheet.toString());
    for (int i=0;i<max_steps;i++)
    {
    	 XSSFRow row = outdoorsSheet.getRow(i+5);
    	 double[] weatherRow = new double[6];
    	 for (int j=0;j<6;j++)
    	 {
    		 XSSFCell cell = row.getCell(j+5);//K
        	 String value = cell.getRawValue();
//        	 System.out.println(value);
        	 weatherRow[j]=Double.parseDouble(value);
    	 }
    	 weatherdata.add(weatherRow);
    }
    
	
//	OUTDOORS INFORMATION										
//	Year 	Month	Day	Hour	Minute	Precipitation	Solar radiation 	External temperature	External relative humidity	Wind speed at 2m	Irrigation
//						[mm]	[W/m²]	[°C]	[%]	[m/s]	[mm]
//	YYYY	MM	DD	HH	MM	P	Sol rad	T	RH	u	Irrig 
//	2017	12	16	0	0	0.0	-2.192	19.6	53.65	0.1	0
//	2017	12	16	0	5	0.0	-2.460	19.9	54.01	0.3	0
//	2017	12	16	0	10	0.0	-2.412	19.9	54.01	0.2	0
//	2017	12	16	0	15	0.0	-2.676	19.8	53.70	0.3	0
//	weatherdata = new ArrayList<double[]>();
//	weatherdata.add(new double[]{2017,	12,	16,	0,	0,	0.0,	-2.192,	19.6,	53.65,	0.1,	0});
//	weatherdata.add(new double[]{2017,	12,	16,	0,	5,	0.0,	-2.460,	19.9,	54.01,	0.3,	0});
//	weatherdata.add(new double[]{2017,	12,	16,	0,	10,	0.0,	-2.412,	19.9,	54.01,	0.2,	0});
//	weatherdata.add(new double[]{2017,	12,	16,	0,	15,	0.0,	-2.676,	19.8,	53.70,	0.3,	0});
//	}
//
//	public void run()
//	{
	//// CREATE PLANT
	Plant plant = new Plant();

//	plant_data = xlsread(input_file,'Plants','E5:E14');
//	double[] plant_data = new double[]{1,
//			1,
//			0.2,
//			0.9,
//			0.5,
//			0.1,
//			0.4,
//			0.4,
//			200,
//			0.001
//		};
	double[] plant_data2 = new double[11];
	XSSFSheet plantsSheet = wb.getSheet("Plants");
//	System.out.println(plantsSheet.toString());
	for (int i=0;i<11;i++)
	{
	   XSSFRow row = plantsSheet.getRow(i+3);
	   XSSFCell cell = row.getCell(4);//E
	   String value = cell.getRawValue();
//	   System.out.println(value);
	   plant_data2[i]=Double.parseDouble(value);
	}
	
	
	plant.LAI=plant_data2[5-4];        //leaf area index A FEB Tabares  
	plant.rho=plant_data2[7-4];     //sw reflectivity plant
	plant.em=plant_data2[8-4];    //emissivity planta
	plant.k=plant_data2[9-4];      // thermal conductivity planta
	plant.height=plant_data2[10-4];   //  height of plant [m]
	plant.ks=plant_data2[11-4]; // Extinsion coefficient   
	plant.ks_ir = plant_data2[12-4]; // IR extinsion coefficient
	plant.rsmin=plant_data2[13-4];      //[s/m] resistencia estomatica 
	plant.z=2;       // [m] height of wind measurements
	plant.Zog = plant_data2[14-4]; // VARIA SEGUN SMOOOTH
//	clear plant_data

	//GF
//	plant_data2=xlsread(input_file,'Plants','E4:E14');
//	double[] plant_data2 = new double[]{100,1,
//			1,
//			0.2,
//			0.9,
//			0.5,
//			0.1,
//			0.4,
//			0.4,
//			200,
//			0.001
//		};

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
	
	
//	double[] sub_data = new double[14];
	XSSFSheet substrateSheet = wb.getSheet("Substrate");
//	System.out.println(substrateSheet.toString());
	for (int i=0;i<14;i++)
	{
	   XSSFRow row = substrateSheet.getRow(i+3);
	   XSSFCell cell = row.getCell(4);//E
	   String value = cell.getRawValue();
//	   System.out.println(value);
	   sub_data[i]=Double.parseDouble(value);
	}

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
	
	XSSFSheet SupportSheet = wb.getSheet("Support");
//	System.out.println(SupportSheet.toString());
	for (int i=0;i<4;i++)
	{
	   XSSFRow row = SupportSheet.getRow(i+3);
	   XSSFCell cell = row.getCell(4);//E
	   String value = cell.getRawValue();
//	   System.out.println(value);
	   roof_data[i]=Double.parseDouble(value);
	}

	roof.depth = roof_data[Constants.ONE]; // meters
	roof.k = roof_data[Constants.TWO]; //thermal conductivity
	roof.density = roof_data[Constants.THREE]; 
	roof.Cp = roof_data[Constants.FOUR]; // heat capacity

//	clear roof_data

	String model_to_use = "Tabares";
	//// CREATE MODEL
	System.out.println("Using model "+model_to_use);
	TabaresThermalMass model = null;
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
	for (int i=Constants.ONE;i<weatherdata.size()-1;i++)
	{
//	    if(use_epw)
//	    {
//	        P(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,28)*dt/Dt;
//	        R(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,6)*0;
//	    }
//	    else
//	    {
		double[] weatherline = weatherdata.get(i);
		int start = (int)(Constants.ONE+Math.round((i-Constants.ONE)*Dt/dt));
		int end = (int)(Math.round((i+1)*Dt/dt)-1);
		// this should go from 1-5 then 6-10 // TODO check starts with zero at first item
		for (int weatherNumber=start;weatherNumber<weatherdata.size()-1;weatherNumber++)
		{
		    P[weatherNumber]=weatherline[Constants.ONE]*dt/Dt;
		    R[weatherNumber]=weatherline[Constants.SIX]*dt/Dt;
		}

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
	double[] runon=new double[P.length];               //Runoff volume entering to the subcatchment (m3)
	double[] runoff=new double[P.length];              //Runoff volume which drains from the subcatchment (m3)
	double[][] Q_out=new double[(int)Math.round(Math.floor(P.length*1.1))][1];      //Outflow from the subcatchment to the street (m3/h) 

	double[][] theta=new double[P.length][nL];              //Soil water content (m3/m3)
	double[][] DthetaDt=new double[P.length][nL];           //Rate of change in soil water content (mm/h)
	double[] f=new double[P.length];                   //Infiltration (mm/h)
	double[][] pe=new double[P.length][nL];                 //Percolation (mm/h)

	double[][] red=new double[P.length][nL];                //Redistribution (mm/h)
	double[] F=new double[P.length];                   //Cumulative infiltration (mm)
	double[] Ft=new double[P.length];                  //Cumulative infiltration to calculate Green Ampt (mm)
	double[]Peffect=new double[P.length];             //Precipitation plus runoff minus interception (mm)
	double[]AWI=new double[P.length];                 //Available water to infiltrate (mm)
	double[]esc=new double[P.length];                 //effective surface runoff (mm)

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
//	    sim_desc = "Performing simulation with inputs from "+input_file+" using "+model_to_use+" model and custom data";
//	}
//	h = waitbar(0,char(sim_desc));


	for (int main_step = Constants.ONE;main_step<weatherdata.size();main_step++)
	{
	    
	    double[] this_data_line = weatherdata.get(main_step);
	    double[] next_data_line = weatherdata.get(main_step+1);
	    
	    for (int k=Constants.ONE;k<n_sub_tsteps;k++)
	    {
	        
	        int tstep = (main_step-Constants.ONE)*n_sub_tsteps+k;
//	        waitbar( tstep/ (max_steps*n_sub_tsteps));
	            
	        // interpolate   
	        double this_inner_t = interior_temperature[main_step];
	        double next_inner_t = interior_temperature[main_step+1];
	        double inner_T = this_inner_t + (k-Constants.ONE)*(next_inner_t - this_inner_t)/n_sub_tsteps;        
	        model.T_interior = inner_T;
	        
	        
	                
	        
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
	        model.moveForward(data_line);
	                        
	        // GF : to be completed by GM
	        result_ET[tstep] = model.et_mm_hour; // Evaporation + Transpiration in mm/hour
	        //E_T=... //GF expressed in mm/h. Must be a matrix (1 line per timestep)
	    
	        // Perform the mass balance and determine the VWC (Volumetric Water
	        // Content)of each layer
	 
	        MassBalance massBalance = new MassBalance();
//	        [theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event]
//	        returnData = 
	        TreeMap massBalanceReturn =
	        		massBalance.MassBalance_(tstep,P,R,runon, result_ET,general_inf, outflow_data, sub_data,plant_data2,theta,DthetaDt,f,pe,red,F,Ft,Peffect, AWI, esc,irr_vol,Ptot_cum_event);
	        theta=(double[][]) massBalanceReturn.get("theta");
	        DthetaDt=(double[][]) massBalanceReturn.get("DthetaDt");
		    f=(double[]) massBalanceReturn.get("f");
		    pe=(double[][]) massBalanceReturn.get("pe");
		    red=(double[][]) massBalanceReturn.get("red");
		    F=(double[]) massBalanceReturn.get("F");
		    Ft=(double[]) massBalanceReturn.get("Ft");
		    Peffect=(double[]) massBalanceReturn.get("Peffect");
		    AWI=(double[]) massBalanceReturn.get("AWI");
		    esc=(double[]) massBalanceReturn.get("esc");
		    irr_vol=(double) massBalanceReturn.get("irr_vol");
		    Ptot_cum_event=(double) massBalanceReturn.get("Ptot_cum_event");
	        
	        
	        // GF : to be modified by GM 
	        model.VWC = GreenRoofCommon.get1DSliceOf2D(theta, 1, tstep);
	        
	       
	        
	        
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

//
	String headers = 
	    "Sky temperature (�C)"+","+
	    "Exterior temperature (�C)"+","+
	    "Wind speed (m/s)"+","+
	    "Global solar radiation (W/m2)"+","+
	    "Substrate Surface Temperature (�C)"+","+
	    "Temperature 5cm deep (�C)"+","+
	    "Temperature 10cm deep (�C)"+","+
	    "Temperature 15cm deep (�C)"+","+
	    "Foliage Temperature (�C)"+","+
	    "VWC surface"+","+
	    "VWC mid depth"+","+
	    "Substrate sensible heat transfer (W/m2)"+","+
	    "Foliage sensible heat transfer (W/m2)"+","+
	    "Foliage Transpiration (W/m2)"+","+
	    "Substrate evaporation (W/m2)"+","+
	    "Rainfall (mm)"+","+
	    "Interior temperature (�C)"+","+
	    "Heat flux below substrate (W/m2)"+","+
	    "Evapotranspiration (mm/hour)"+","+
	    "plant_absorbed_solar"+","+
	    "plant_absorbed_ir_sky"+","+  
	    "Qir"+","+
	    "substrate_solar_radiation"+","+
	    "substrate_infrared_radiation"+","+
	    "substrate_conduction"+","+
	    "heating_load (W/m2)"
	    ;
	System.out.println(headers);
	for (int i=0;i<n_sub_tsteps;i++)
	{	
		String results = 
		    (result_T_sky[i]-273.15)+","+
		    (result_T_out[i]-273.15)+","+
		    (result_wind_speed[i])+","+
		    (result_Rsh[i])+","+
		    (result_T_substrate[i]-273.15)+","+
		    (result_T_5cm[i]-273.15)+","+
		    (result_T_10cm[i]-273.15)+","+
		    (result_T_15cm[i]-273.15)+","+
		    (result_T_plants[i]-273.15)+","+
		    (result_VWC_surface[i])+","+
		    (result_VWC_mid[i])+","+
		    (result_substrate_convection[i])+","+
		    (result_plant_convection[i])+","+
		    (result_transpiration[i])+","+
		    (result_evaporation[i])+","+
		    (result_Rain[i])+","+
		    (result_T_interior[i]-273.15)+","+
		    (result_interface_heat_flux[i])+","+
		    (result_ET[i])+","+
		    (result_plant_absorbed_solar[i])+","+
		    (result_plant_absorbed_ir_sky[i])+","+  
		    (result_Qir[i])+","+
		    (result_substrate_solar_radiation[i])+","+
		    (result_substrate_infrared_radiation[i])+","+
		    (result_substrate_conduction[i])+","+
		    (result_heating_load[i]);
		System.out.println(results);
	}
//	xlswrite(char(out_file),[headers;results(1:n_sub_tsteps:end,:)],char(model_to_use));
//
	}



}
