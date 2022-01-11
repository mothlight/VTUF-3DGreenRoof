package greenroof;

import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
//import org.apache.commons.math3.linear.Array2DRowRealMatrix;
//import org.apache.commons.math3.linear.MatrixUtils;
//import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import jml.matlab.Matlab;

public class TabaresThermalMass
{
	

        //% Current state        
        public double T_plants;
        public double T_substrate;
        public double T_interior ;       
        public double[] VWC; // volumetric water content     
        
        public double interface_heat_flux;
        
        // Plants balance
        public double plant_absorbed_solar;
        public double plant_absorbed_ir_sky        ;
        public double transpiration;
        public double plant_convection  ; 
        
        // Both
        public double Qir;
                
        // Substrate balance
        public double substrate_solar_radiation;
        public double substrate_infrared_radiation;
        public double substrate_convection;
        public double evaporation;
        public double substrate_conduction ;
        
        
        public double[] innerT;
        public double dt;
        
        public double et_mm_hour;
        
        public double heating_load;
        public double rs;
        public double fsolar;
        public double fvpd;
        public double fvwc;
        public double ftemp;
        public double ra;
        
        //% Materials
        public Roof roof;
        public Plant plant;
        public Substrate sub;
        
        //% Matrices for thermal mass
        public double[][] C;
        public double[][] K;
        
        //% other                
        public double L1;
        public int n_roof_nodes; // nnodes on the roof and other nnodes in the substrate
        public int n_sub_nodes;
        public double[] samples;
    
    
   public TabaresThermalMass(Plant plant, Substrate sub, Roof roof, int n_substrate_layers,int n_roof_layers, double Area, double dt)
   {
            
            this.dt = dt*60.;
            
            this.roof = roof;
            this.plant = plant;
            this.sub = sub;
            
            //% Initialize
            T_plants = 285.; // K
            T_substrate = 291.; // K
            T_interior = 293.;
            
            
            //% Other      
            n_sub_nodes = n_substrate_layers;
            n_roof_nodes = n_roof_layers;
            L1 = Math.sqrt(Area);            
            
            double dxSubstrate = sub.depth/n_sub_nodes;
            double dxSupport = roof.depth/n_roof_nodes;
            
            ArrayList<Double> samplesArray = new ArrayList<Double>();
            
//            samples = new double[14];
            
            double samples0 = 0;
            double samples1 = dxSubstrate/2.;//0.15
            double samples2 = dxSubstrate;//0.03
            double samples3 = (sub.depth-dxSubstrate/2.);//0.0135
            double samples4 = (sub.depth+dxSupport/2.);//0.165
            double samples5 = dxSupport;//0.03
            double samples6 = (roof.depth+sub.depth-dxSupport/2.);//0.285
            double samples7 = sub.depth+roof.depth;//0.3
            
//            samples[0] = samples0;           //0.0//
//            samples[1] = samples1;           // 0.015000//
//            samples[2] = samples1;// 0.045000
//            samples[3] = samples2;// 0.075000
//            samples[4] = samples3;// 0.10500
//            samples[5] = samples3;           // 0.13500//
//            samples[6] = samples4;           // 0.16500//
//            samples[7] = samples6;// 0.19500
//            samples[8] = samples7;// 0.22500
//            samples[9] = samples2;// 0.25500
//            samples[10] = samples6;         // 0.28500//
//            samples[11] = samples7;         // 0.30000//
            
//            samples[0]=samples0;  
            samplesArray.add(samples0);
            double startLoop = samples1;
            double loopValue = startLoop;
            double endLoop = samples3;
//            int count = 1;
            //this is to deal with the colons, i.e. 1:2:10 -> 1 3 5 7 9
            while (true)
            {
//            	samples[count]=loopValue;
            	samplesArray.add(loopValue);
            	loopValue+=samples2;
//            	count++;
            	//floating point precision, round off
               	if (Math.abs(loopValue-endLoop) < 0.000001)
            	{
            		continue;
            	}
            	if (loopValue > endLoop)
            	{
            		break;
            	}
            }
            
            startLoop = samples4;
            loopValue = startLoop;
            endLoop = samples6;            
            while (true)
            {
//            	samples[count]=loopValue;
            	samplesArray.add(loopValue);
            	loopValue+=samples5;
//            	count++;
            	//floating point precision, round off
            	if (Math.abs(loopValue-endLoop) < 0.000001)
            	{
            		continue;
            	}
            	if (loopValue > endLoop)
            	{
            		break;
            	}
            }
//            samples[count]=samples7;
            samplesArray.add(samples7);
  
              //             0.015         0.03         0.135                     0.165                 0.03        0.285
              //             0.015  
//            samples = [0 dxSubstrate/2:dxSubstrate:(sub.depth-dxSubstrate/2) (sub.depth+dxSupport/2):dxSupport:(roof.depth+sub.depth-dxSupport/2) sub.depth+roof.depth];  
            
            samples = GreenRoofCommon.toArray(samplesArray);
            
//        	disp(   [0 (sub.depth+roof.depth)]  );  // 0.00000   0.30000
//        	disp(  [obj.T_substrate obj.T_interior]  );  //   291   293
//        	disp(  samples(2:end-1)  );
//        	disp( interp1([0 (sub.depth+roof.depth)],[obj.T_substrate obj.T_interior],obj.samples(2:end-1))'  );
            double interp0 = 0.0;
            double interp1 = (sub.depth+roof.depth);
            double interp2 = T_substrate;
            double interp3 = T_interior;
            double[] interp0Array = new double[]{interp0,interp1};
            double[] interp1Array = new double[]{interp2,interp3};
            double[] interp2Array =  GreenRoofCommon.subsetArray(samples, 1, samples.length-2);
            
            //% Initialize inner temperatures by interpolating
//            innerT = interp1([0 (sub.depth+roof.depth)],[T_substrate T_interior],samples(2:end-1))';  
            innerT = GreenRoofCommon.interpLinear(interp0Array, interp1Array, interp2Array);
//            System.out.println(innerT);
            
            // Create matrices
            //obj = setMatrices();
            
   }
        
   public void setMatrices()
   {
            // Update the matrices
            
            //Create empty matrix
            int total_nodes = n_roof_nodes + n_sub_nodes;
            C = new double[total_nodes][total_nodes];
            K = new double[total_nodes][total_nodes];
                     
            // Define parameters            
            int n = Constants.ONE; //Helper
            
            // Roof properties do not change
            double mcRoof = (roof.density*roof.depth/n_roof_nodes)*roof.Cp;            
            double rRoof = roof.depth/roof.k/n_roof_nodes;
            
            // Connect before substrate (connection with T_substrate)
            double rSub = sub.depth/(sub.thermalConductivity(VWC[Constants.ONE])*n_sub_nodes);
            K[Constants.ONE][Constants.ONE] = 2./rSub;

            
            double mcSub=0;
            double rSub1;
            double rSub2;
            
            // Connect within the substate   
            for (int i=Constants.ONE;i<n_sub_nodes-1;i++)
            {
                mcSub = (sub.depth/n_sub_nodes)*sub.rhoCp(VWC[n]);          
                rSub1 = sub.depth/(sub.thermalConductivity(VWC[n])*n_sub_nodes);
                rSub2 = sub.depth/(sub.thermalConductivity(VWC[n+1])*n_sub_nodes);
                
                int[] dim1 = new int[]{n,n+1}; 
                int[] dim2 = new int[]{n,n+1};            	
            	
                C[n][n] = mcSub;
//                K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]./(rSub1/2+rSub2/2);

//                double[][] Ktemp = GreenRoofCommon.subset2DArray(K, dim1[0], dim1[1], dim2[0], dim2[1]);
//                double[][] oneArray = new double[][]{{1,-1},{-1,1}};
//                double[][] KtempPlusOneDivide = GreenRoofCommon.divideArray(oneArray, (rSub1/2.+rSub2/2.));
//                double[][] KtempPlusOne = GreenRoofCommon.addArrays(Ktemp, KtempPlusOneDivide);               
                double[][] KtempPlusOneC = GreenRoofCommon.addArrays(GreenRoofCommon.subset2DArray(K, dim1[0], dim1[1], dim2[0], dim2[1]), 
                										GreenRoofCommon.divideArray(new double[][]{{1,-1},{-1,1}}, (rSub1/2.+rSub2/2.)));                
                K = GreenRoofCommon.insertArray(K, KtempPlusOneC, n, n);          
                
//                K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]/(rSub1/2.+rSub2/2.);
                n = n+1;
            }
            
            // Connect interface between both
            rSub = sub.depth/(sub.thermalConductivity(VWC[n])*n_sub_nodes);
            C[n][n] = mcSub;
            int[] dim1b = new int[]{n,n+1}; 
            int[] dim2b = new int[]{n,n+1};   
//            K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]./(rSub/2. + rRoof/2.);
            double[][] KtempPlusOne = GreenRoofCommon.addArrays(GreenRoofCommon.subset2DArray(K, dim1b[0], dim1b[1], dim2b[0], dim2b[1]), 
					GreenRoofCommon.divideArray(new double[][]{{1,-1},{-1,1}}, (rSub/2.+rRoof/2.)));                
            K = GreenRoofCommon.insertArray(K, KtempPlusOne, n, n); 

            n=n+1;
            
            // Connect within the roof            
            for (int i=Constants.ONE;i<n_roof_nodes-1;i++)
            {
                int[] dim1c = new int[]{n,n+1}; 
                int[] dim2c = new int[]{n,n+1};  
                C[n][n] = mcRoof;
//                K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]./rRoof;
                double[][] KtempPlusOneB = GreenRoofCommon.addArrays(GreenRoofCommon.subset2DArray(K, dim1c[0], dim1c[1], dim2c[0], dim2c[1]), 
    					GreenRoofCommon.divideArray(new double[][]{{1,-1},{-1,1}}, rRoof));                
                K = GreenRoofCommon.insertArray(K, KtempPlusOneB, n, n); 

                n = n+1;
            }
            
            // Connect final one.
            C[n][n] = mcRoof;
            //TODO, I don't think this does anything
//            consts = Constants;
            K[n][n] = K[n][n]+1/(Constants.rsi_roof + rRoof/2.);

   }
//        
        public void update(WeatherDataLine data)
        {
            
            // Load constants

            Air air = new Air();
            
            // INPUTS
            double R_sh=   data.R_sh;  // incoming SW radiation global horizontal
            double Tair=   data.Tair;
            double Tsky=   data.Tsky;
            double RH=     data.RH;   // Relative humidity [%]
            //T_interior =   data.T_interior;   // interior surface temperature
            double Pa=     data.Pa;   // atmospheric pressure [Pa]
            double U=      Math.max(3,data.U);       // Wind speed [m/s]
            //W=      VWC(1);    //VWC            
            
            double[] T = new double[]{T_plants,T_substrate};
            
            // Derived inputs
            double Mu=air.kinematicViscosity(Tair);   // ***** viscosidad cinematica
            double dens_a= air.density(Pa,Tair);   // density air            
            double Tfilm=  (Tair+T[Constants.ONE])/2.;
            double Beta=   1./Tfilm;
            
            
            // OTHER            
            // From table 5.
            double kpor  = sub.phi*Air.k+(1.-sub.phi)*plant.k;            
            // Auxiliar variable
            double Mg=     VWC[Constants.ONE]/sub.VWCsat;        
            // No idea where this comes from.
            double alphapor= kpor/(dens_a*Air.Cp);
            
            // eq. 6 Tabares ... constants may vary?
            double rsub=   34.52*Math.pow(Mg,(-3.2678));
            
            
            //   DIMENSIONLESS NUMBERS
            // Reynolds **note* with L1 not L
            double Re=     dens_a*U*L1/Mu;
            // Prandtl
            double Pr=     Air.Cp*Mu/Air.k;
            // Grashof
            double Gr=     Math.abs(Constants.g*Beta*Math.pow(dens_a,2)*(T[Constants.ONE]-Tair)*(Math.pow(L1,3))/(Math.pow(Mu,2)));
            // Raleigh
            double Ra=     Gr*Pr;
            // Lewis
            double Le=     1.;
            
            // Nusselt. eq. 4 
            double Nu;
            if (Gr<(0.068*Math.pow(Re,2.2)))
            {
                Nu= 3.+1.25*0.025*Math.pow(Re,0.8);
            }
            else if ((Gr>(0.068*Math.pow(Re,2.2)))&&(Gr<(55.3*Math.pow(Re,(5./3.)))))
            {
                Nu= 2.7*(Math.pow((Gr/(Math.pow(Re,2.2))),(1./3.)))*(3.*(15./4.)+(15./16.)*0.0253*Math.pow(Re,0.8));
            }
            else
            {
                Nu= 0.15*(Math.pow(Ra,(1./3.)));
            }
            
            // available in Table 5. Tabares
            double Pe=     0.3*L1*U/alphapor; // CHECK
            double hpor  = kpor*1.128*Math.pow(Pe,0.5)/L1; // CHECK            
            double hconv= 15*Nu*Air.k/L1; // plantBB was 1.5 then 15
            double ras=  dens_a*Air.Cp*(Math.pow(Le,(2./3.)))/hconv;                   
            double ra=dens_a*Air.Cp*(Math.pow(Le,(2./3.)))/hconv;  //ORIGINAL
            ra = ra;
            double hsub =  hpor*hconv/(hpor+hconv); // CHECK
            
            // Vapor Pressures [kPa]
            
            double e_s  = 610.8*Math.exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000.;     
            double e_sf = 610.8*Math.exp(17.27*(T[Constants.ONE]-273.15)/(T[Constants.ONE]-273.15+237.3))/1000.;
            double e_ss = 610.8*Math.exp(17.27*(T[Constants.TWO]-273.15)/(T[Constants.TWO]-273.15+237.3))/1000.;   
            double e_air= e_s*RH;

                                                          
            // SHORT WAVE RADIATION
                
            // eq. 12 Tabares
            plant_absorbed_solar = (1.-plant.foliage_rho()-plant.tau_fsol())*(1.+plant.tau_fsol()*sub.rho)*R_sh;
            
            // eq. 13 Tabares
            substrate_solar_radiation = plant.tau_fsol()*(1.-sub.rho)*R_sh;
            
            // LONG WAVE RADIATION
         
            // eq. 14 Tabares
            plant_absorbed_ir_sky = (1.-plant.tau_fir())*plant.em*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(Tsky,4));
            
            // eq. 15 Tabares
            substrate_infrared_radiation =   -(plant.tau_fir())*sub.em*Constants.SB*((Math.pow(T[Constants.TWO],4)-Math.pow(Tsky,4)));
            
            // eq. 17 Tabares
            double em_1 = (1./sub.em)+(1./plant.em)-1.;
            Qir = (1.-plant.tau_fir())*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(T[Constants.TWO],4))/em_1;
            
            // CONVECTION                            
            // eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            plant_convection=   -1.5* plant.LAI*hconv*(T[Constants.ONE]-Tair);
            
            // eq. 19 Tabares
            substrate_convection=    -hsub*(T[Constants.TWO]-Tair);
            
                
            // EVAPOTRANSPIRATION                               
        
            // eq. 22 Tabares
            double f_sol=  1.+Math.exp(-0.034*(R_sh-3.5));
            
            // eq. 24 Tabares... extended?
            //f_VPD=  1/(1-0.41*log(e_sf-e_air));
            double VPD=e_ss-e_air;
            double VPD_f=e_sf-e_air;
            double f_VPD;

            if (VPD_f > 0)
            {
                f_VPD=  (1.-0.41*Math.log(VPD_f));
            }
            else
            {
                f_VPD= 1.;
            }
            if (f_VPD>1)
            {
                f_VPD=1.;
            }
            else if (f_VPD<0)
            {
                f_VPD=0.05;
            }

            // eq. 25 Tabares
            double f_temp= 1./(1-0.0016* Math.pow((35-(T[Constants.ONE]-273.15)),2) );
            if (f_temp < 0)
            {
                f_temp = 10000.;
            }
            // eq. 23 Tabares
            double rootW = getVWC(sub.depth/2.);
            double f_W;
            if ( rootW >0.7*sub.VWCfc) // W_fc vs W_sat
            {
                f_W=    1.;
            }
            else if ((rootW <= 0.7*sub.VWCfc) && (rootW > sub.VWCresidual))
            {
                f_W=    (0.7*sub.VWCfc-sub.VWCresidual)/(rootW-sub.VWCresidual);
            }
            else
            {
                f_W=    1000.;
            }
            
            
            // eq. 21 Tabares
            //rs = plant.rsmin*f_sol*f_VPD*f_W*f_temp/plant.LAI;
            double f_hum=1./(f_VPD);    
            rs = (plant.rsmin/plant.LAI)*f_sol*f_hum*f_W*f_temp;
            rs = rs;
            fsolar=f_sol;
            fvpd=f_hum;
            fvwc=f_W;
            ftemp=f_temp;
            
            // eq. 20 Tabares
            //Qt= plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            transpiration= -plant.LAI*dens_a*Air.Cp*(e_sf-e_air)/(Air.gamma*(rs+ra)); 
            
            // eq. 5 Tabares
            //evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            evaporation= -dens_a*Air.Cp*VPD/(Air.gamma*(rsub+ras));                               
        
            // Eq. 16 - Sailor 2008... the second part is from EPlus code
            double Tg = T_substrate;
            double Tf = T_plants;
            double Lef = 1.91846e6*Math.pow((Tf/(Tf-33.91)),2) ;
            if(T_plants < 273.15) // Less than 0 C
            {
                Lef = 2.838e6;
            }
            
            double Leg = 1.91846e6*Math.pow((Tg/(Tg-33.91)),2);
            if(T_plants < 273.15) // Less than 0 C
            {
                Lef = 2.838e6;
            }
            et_mm_hour = -(evaporation/Leg + transpiration/Lef)*3600;   
            
            
                                             
            // Update inner temperatures
            // find temperatures within the subtrate/concrete
            double rSub = (sub.depth/sub.thermalConductivity(VWC[Constants.ONE]))/n_sub_nodes;
            double rRoof = (roof.depth/roof.k)/n_roof_nodes;
                        
            double[] f = new double[innerT.length];
            f[Constants.ONE]   = 2.*T_substrate/rSub;
            f[f.length-1] = T_interior/(Constants.rsi_roof+rRoof/2);
            //innerT = inv(C/dt + K)*(C/dt*innerT + f);
            //innerT = innerT + C\(f - K*innerT)*dt;
            int total_nodes = n_roof_nodes + n_sub_nodes;
            
            double[][] Kdt = GreenRoofCommon.multiply(K, dt);
            RealMatrix matrixKdt = new BlockRealMatrix(Kdt); 
            RealMatrix eyeTotalNodes = Matlab.eye(total_nodes,total_nodes);
            RealMatrix matrixC = new BlockRealMatrix(C);
            RealMatrix matrixCDivKdt = Matlab.mldivide(matrixC, matrixKdt);
            RealMatrix matrixR = eyeTotalNodes.add(matrixCDivKdt);
            
            
            
//            double[][] dtK = GreenRoofCommon.multiply(K, dt);
//            RealMatrix matrixC = new Array2DRowRealMatrix(C);
//            RealMatrix matrixdtK = new Array2DRowRealMatrix(dtK);
//            RealMatrix matrixdtKInv = MatrixUtils.inverse(matrixdtK);
//            RealMatrix matrixCDTK = matrixC.multiply(matrixdtKInv);
//            
//            RealMatrix eyeTotalNodes = GreenRoofCommon.eyeMatrix(total_nodes,total_nodes);
//            RealMatrix matrixR = eyeTotalNodes.add(matrixCDTK);
            
//            double[][] R = GreenRoofCommon.eyeMatrix(total_nodes,total_nodes)+(C\(dtK)); 
            
//            double[] fDt = GreenRoofCommon.multiply(f, dt);
//            RealMatrix matrixfDtInv = MatrixUtils.inverse(matrixdtK);
//            RealMatrix matrixCfDt = matrixC.multiply(matrixfDtInv);
            
            
//            double[] fdt = GreenRoofCommon.multiply(f, dt);
//            RealMatrix matrixfdt = new BlockRealMatrix(fdt); 
            
            
            double[] fdt = GreenRoofCommon.multiply(f, dt);
            RealMatrix matrixfdt = new Array2DRowRealMatrix(fdt); 
            RealMatrix matrixCDivfdt = Matlab.mldivide(matrixC, matrixfdt);
            RealMatrix matrixinnerT = new Array2DRowRealMatrix(innerT); 
            RealMatrix innerTPlusCDivfdt = matrixinnerT.add(matrixCDivfdt);
            RealMatrix iT = Matlab.mldivide(matrixR, innerTPlusCDivfdt);
//            System.out.println(iT);
            double[][] innerTData = iT.getData();
            for (int i=0;i<innerT.length;i++)
            {
            	innerT[i]=innerTData[i][0];
            }
//            innerT  = R\(innerT+(C\f*dt)); 
            
            double Qcond = 2.*(T[Constants.TWO] - innerT[Constants.ONE])/rSub;              
            substrate_conduction = -Qcond;
            
            // Update interface_heat_flux            
            rSub = (sub.depth/sub.thermalConductivity(VWC[n_sub_nodes-1]))/n_sub_nodes; // Last node
//            double deltaT = innerT[n_sub_nodes-1]-innerT[n_sub_nodes+1-1]; //one node on each.  
            double deltaT = innerT[n_sub_nodes-1]-innerT[n_sub_nodes+1-1]; //one node on each.   
            interface_heat_flux = 2*deltaT/(rSub+rRoof);
            
            double Tlosa = getTemperature(sub.depth + roof.depth);
            heating_load = ( T_interior - Tlosa)/Constants.rsi_roof;
            
            
        }
        
        public double[] ResFUN(double[] T, WeatherDataLine data)
        {                       
            // Load constants

            Air air = new Air();
            
            // INPUTS
            double R_sh=   data.R_sh;  // incoming SW radiation global horizontal
            double Tair=   data.Tair;
            double Tsky=   data.Tsky;
            double RH=     data.RH;   // Relative humidity [%]
            //T_interior =   data.T_interior;   // interior surface temperature
            double Pa=     data.Pa;   // atmospheric pressure [Pa]
            double U=      Math.max(3,data.U);       // Wind speed [m/s]
            //W=      VWC(1);    //VWC            
            
            // Derived inputs
            double Mu=air.kinematicViscosity(Tair);   // ***** viscosidad cinematica
            double dens_a= air.density(Pa,Tair);   // density air            
            double Tfilm=  (Tair+T[Constants.ONE])/2;
            double Beta=   1./Tfilm;            
            
            // OTHER            
            // From table 5.
            double kpor  = sub.phi*Air.k+(1.-sub.phi)*plant.k;            
            // Auxiliar variable
            double Mg=     VWC[Constants.ONE]/sub.VWCsat;        
            // No idea where this comes from.
            double alphapor= kpor/(dens_a*Air.Cp);
            
            // eq. 6 Tabares ... constants may vary?
            double rsub=   34.52*Math.pow(Mg,(-3.2678));
            
            
            //   DIMENSIONLESS NUMBERS
            // Reynolds **note* with L1 not L
            double Re=     dens_a*U*L1/Mu;
            // Prandtl
            double Pr=     Air.Cp*Mu/Air.k;
            // Grashof
            double Gr=     Math.abs(Constants.g*Beta*Math.pow(dens_a,2)*(T[Constants.ONE]-Tair)*(Math.pow(L1,3))/(Math.pow(Mu,2)));
            // Raleigh
            double Ra=     Gr*Pr;
            // Lewis
            double Le=     1.;
            double Nu;
            
            // Nusselt. eq. 4 
            if (Gr<(0.068*Math.pow(Re,2.2)))
            {
                Nu= 3.+1.25*0.025*Math.pow(Re,0.8);
            }
            else if ((Gr>(0.068*Math.pow(Re,2.2)))&&(Gr<(55.3*Math.pow(Re,(5./3.)))))
            {
                Nu= 2.7*( Math.pow((Gr/(Math.pow(Re,2.2))),(1./.3)) )*(3.*(15./4.)+(15./16.)*0.0253*Math.pow(Re,0.8));
            }
            else
            {
                Nu= 0.15*(Math.pow(Re,(1./3.)));
            }
            
            // available in Table 5. Tabares
            double Pe=     0.3*L1*U/alphapor; // CHECK
            double hpor  = kpor*1.128*Math.pow(Pe,0.5)/L1; // CHECK            
            double hconv= 15*Nu*Air.k/L1; // plantBB was 1.5 then 15
            double ras=  dens_a*Air.Cp*( Math.pow(Le,(2./3.)) )/hconv;                   
            double ra=dens_a*Air.Cp*( Math.pow(Le,(2./3.)) )/hconv;  //ORIGINAL
            double hsub =  hpor*hconv/(hpor+hconv); // CHECK
            
            // Vapor Pressures [kPa]
            
            double e_s  = 610.8*Math.exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000.;     
            double e_sf = 610.8*Math.exp(17.27*(T[Constants.ONE]-273.15)/(T[Constants.ONE]-273.15+237.3))/1000.;
            double e_ss = 610.8*Math.exp(17.27*(T[Constants.TWO]-273.15)/(T[Constants.TWO]-273.15+237.3))/1000.;   
            double e_air= e_s*RH;

                                                          
            // SHORT WAVE RADIATION
            
            // eq. 12 Tabares
            double Rsh_f = (1.-plant.foliage_rho()-plant.tau_fsol())*(1.+plant.tau_fsol()*sub.rho)*R_sh;
            
            // eq. 13 Tabares
            double Rsh_s = plant.tau_fsol()*(1.-sub.rho)*R_sh;
            
            // LONG WAVE RADIATION
            
            // eq. 14 Tabares
            double Qir_f= (1.-plant.tau_fir())*plant.em*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(Tsky,4));
            
            // eq. 15 Tabares
            double Qir_scov=   (plant.tau_fir())*sub.em*Constants.SB*((Math.pow(T[Constants.TWO],4)-Math.pow(Tsky,4)));
            
            // eq. 17 Tabares
            double em_1 = (1./sub.em)+(1/plant.em)-1.;
            double Qir_sp = (1.-plant.tau_fir())*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(T[Constants.TWO],4))/em_1;
            
            // CONVECTION
            
            // eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            plant_convection=   1.5* plant.LAI*hconv*(T[Constants.ONE]-Tair);
            
            // eq. 19 Tabares
            substrate_convection=    hsub*(T[Constants.TWO]-Tair);
            
            // EVAPOTRANSPIRATION
            
            // eq. 22 Tabares
            double f_sol=  1.+Math.exp(-0.034*(R_sh-3.5));
            
            // eq. 24 Tabares... extended?
            //f_VPD=  1/(1-0.41*log(e_sf-e_air));
            double VPD=e_ss-e_air;
            double VPD_f=e_sf-e_air;
            double f_VPD;

            if (VPD_f > 0)
            {
                f_VPD=  (1.-0.41*Math.log(VPD_f));
            }
            else
            {
                f_VPD= 1.;
            }
            if (f_VPD>1)
            {
                f_VPD=1.;
            }
            else if (f_VPD<0)
            {
                f_VPD=0.05;
            }

            // eq. 25 Tabares
            double f_temp= 1./(1.-0.0016* Math.pow((35.-(T[Constants.ONE]-273.15)),2) );
            if (f_temp < 0)
            {
                f_temp = 10000;
            }
            
            // eq. 23 Tabares
            double rootW = getVWC(sub.depth/2.);
            double f_W;
            if (rootW > 0.7*sub.VWCfc) // W_fc vs W_sat
            {
                f_W=    1.;
            }
            else if ((rootW < 0.7*sub.VWCfc) && (rootW > sub.VWCresidual))
            {
                f_W=    (0.7*sub.VWCfc-sub.VWCresidual)/(rootW-sub.VWCresidual);
            }
            else
            {
                f_W=    1000.;
            }
            
            
            // eq. 21 Tabares
            //rs = plant.rsmin*f_sol*f_VPD*f_W*f_temp/plant.LAI;
            double f_hum=1./(f_VPD);    
            rs = (plant.rsmin/plant.LAI)*f_sol*f_hum*f_W*f_temp;

            // eq. 20 Tabares
            //Qt= plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            double Qt= plant.LAI*dens_a*Air.Cp*(e_sf-e_air)/(Air.gamma*(rs+ra)); 
            
            // eq. 5 Tabares
            //evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            evaporation= dens_a*Air.Cp*VPD/(Air.gamma*(rsub+ras));    
            
            // CONDUCTION
                                    
            // eq. 8 Tabares modified to only consider the conduction
            // through the first layer of substrate. 
            
            // find temperatures within the subtrate/concrete            
            double rSub = sub.depth/sub.thermalConductivity(VWC[Constants.ONE])/n_sub_nodes; // Using thermal conductivity of the first node.
            double rRoof = roof.depth/roof.k/n_roof_nodes;
            
            double[] f = new double[innerT.length];
            f[Constants.ONE] =  2.*T[Constants.TWO]/rSub;
            f[f.length-1]= T_interior/(Constants.rsi_roof+rRoof/2.);
           
            //iT = inv(C/dt - K)*((C/dt)*innerT - f);
            //iT = innerT + C\(f - K*innerT)*dt;
            int total_nodes = n_roof_nodes + n_sub_nodes;
            
//            System.out.println(C);
//            System.out.println(dt);
//            System.out.println(K);
//            System.out.println();
            
//            double[][] R = GreenRoofCommon.eye(total_nodes,total_nodes)+(C\(dt*K));   
            
          
            double[][] Kdt = GreenRoofCommon.multiply(K, dt);
            RealMatrix matrixKdt = new BlockRealMatrix(Kdt);            
            RealMatrix eyeTotalNodes = Matlab.eye(total_nodes,total_nodes);
            RealMatrix matrixC = new BlockRealMatrix(C);
            RealMatrix matrixCDivKdt = Matlab.mldivide(matrixC, matrixKdt);
            RealMatrix matrixR = eyeTotalNodes.add(matrixCDivKdt);
//            System.out.println(matrixCDivKdt);
//            System.out.println(matrixR);
            
//            double[][] CDivKdt = GreenRoofCommon.divideArrays(C, Kdt);
//            RealMatrix matrixCDivKdt = new Array2DRowRealMatrix(CDivKdt);
//            System.out.println(matrixCDivKdt);
//            
//            RealMatrix matrixKdt = new Array2DRowRealMatrix(Kdt);
//            RealMatrix matrixC = new Array2DRowRealMatrix(C);
//            RealMatrix matrixKdtInv = MatrixUtils.inverse(matrixKdt);
//            RealMatrix matrixCKdtInv = matrixC.multiply(matrixKdtInv);
//            RealMatrix matrixR = eyeTotalNodes.add(matrixCKdtInv);
            
            
//            double iT = R\(innerT+(C\f*dt)); 
            double[] fdt = GreenRoofCommon.multiply(f, dt);
            RealMatrix matrixfdt = new Array2DRowRealMatrix(fdt); 
            RealMatrix matrixCDivfdt = Matlab.mldivide(matrixC, matrixfdt);
            RealMatrix matrixinnerT = new Array2DRowRealMatrix(innerT); 
            RealMatrix innerTPlusCDivfdt = matrixinnerT.add(matrixCDivfdt);
            RealMatrix iT = Matlab.mldivide(matrixR, innerTPlusCDivfdt);
//            System.out.println(iT);
            double[][] innerTData = iT.getData();
            for (int i=0;i<innerT.length;i++)
            {
            	innerT[i]=innerTData[i][0];
            }
            double iT0 = innerTData[Constants.ONE][0];
            
            double Qcond = 2.*(T[Constants.TWO] - iT0)/rSub;
            
            // ENERGY (& MASS) BALANCE
            
            // Foliage energy balance
            double Ef= - Rsh_f + plant_convection + Qir_f + Qt + Qir_sp; // Eq 14 Camilo
            
            // Substrate energy balance
            double Es= - Rsh_s - Qir_sp + substrate_convection + Qir_scov + evaporation + Qcond; // Eq 15 Camilo
                    
           
            // RETURN
            double[] Res=new double[]{Ef, Es} ;
            return Res;
  
            
        } // end resFUN
        
   
        
        
        public void moveForward(WeatherDataLine data)
        {
            setMatrices();
            
//            options=optimset('Display','off','TolFun',1e-6);
            double[] Tguess = new double[]{T_plants, T_substrate};                                                      
//            Tss = (fsolve(@(TT) ResFUN(TT,data),Tguess,options));
            
            //testing value
//            double[] T = new double[]{285,291};
            this.a = Tguess[0];
            this.b = Tguess[1];
            double[] Tss;            
            double adjust =0.1;
            int count = 0;            
            double center = (Tguess[0]+Tguess[1])/2.0;
            Tguess[0]=center;
            Tguess[1]=center;    
            ArrayList<Double> previousDifferences = new ArrayList<Double>();
            previousDifferences.add(99999.);
            previousDifferences.add(99998.);
            boolean oss = false;
            int ossCount = 0;
            //this is all working yet, skip for now.
            while (true)
            {
            	if (true)
            	{
            		Tss = new double[]{292.61   ,287.96};
            		break;
            	}
            	
            	double returnValue = calcExpression( Tguess, data);
            	previousDifferences.add(returnValue);

            	System.out.println(returnValue + " " + Tguess[0] + " " + Tguess[1] ); 
            	if (returnValue < 1)
            	{
            		adjust = getAdjust(returnValue);
            	}                   	
            	if (returnValue < previousDifferences.get(previousDifferences.size()-2) &&   previousDifferences.get(previousDifferences.size()-2) < previousDifferences.get(previousDifferences.size()-3))
            	{
            		if (oss)
            		{
            			Tguess[1] = Tguess[1] - adjust;
            		}
            		else
            		{
            			Tguess[0] = Tguess[0] + adjust;
            		}              	                	
            	}
            	else
            	{
            		if (oss)
            		{
            			Tguess[1] = Tguess[1] + adjust;
            		}
            		else
            		{
            			Tguess[0] = Tguess[0] - adjust;
            		}                               	
            	} 
  

            	double prev3 = previousDifferences.get(previousDifferences.size()-3) ;
            	double prev1 = previousDifferences.get(previousDifferences.size()-1);
            	if (prev3 == prev1)
            	{
            		ossCount ++;
            		System.out.println("oss");
            		if (oss)
            		{
//            			oss=false;
//            			adjust = getAdjust(returnValue);
            		}
            		else
            		{
            			oss=true;
//            			adjust = getAdjust(returnValue);
            		}
            	}
            	
            	Tss = Tguess;
            	count++;
            	
                if (returnValue < 0.01)
                {
                	break;
                }
            }    
            
            
//            Tss = bisectionalSolver(Tguess, data);
           

            
//            double[] Tss = (fsolve(@(TT) ResFUN(TT,data),Tguess,options)); //TODO, the mess above still doesn't really work
            T_plants=(Tss[Constants.ONE]);
            T_substrate=(Tss[Constants.TWO]); 
            update(data);   
            
        }
        
        public double getAdjust(double returnValue)
        {
        	double adjust;
        	
        	if (returnValue < 0.1)
        	{
        		adjust = Math.random()*0.0001;
        	} 
         	else if (returnValue < 1)
        	{
        		adjust = Math.random()*0.001;
        	}   
        	else if (returnValue < 5)
        	{
        		adjust = Math.random()*0.01;
        	}  
        	else 
        	{
        		adjust = Math.random()*0.1;
        	} 
        	return adjust;
        }
        
        public double calcExpression(double[] Tguess, WeatherDataLine data)
        {
        	double[] Tss = ResFUN(Tguess, data);
        	double difference = Math.abs(Tss[0]-Tss[1]);
//        	System.out.println(Tss[0]+" "+Tss[1] + " " + difference);
        	return difference;
        }
        
        double a;
        double b;
    	public static final int ITERATIONS = 10000;
    	boolean DEBUG=true;
    	boolean SKIP=true;
    	public static final double[] ERROR_RETURN = {-9999.};
    	FSolve fsolve = new FSolve();    	
    	double TOL = 1e-6;
        
//    	// adapted from https://x-engineer.org/bisection-method/
//    	public double[] bisectionalSolver(double[] Tguess, WeatherDataLine data)
//    	{
//    		double localA1 = Tguess[0];
//    		double localA2 = Tguess[0];
//    		double localB1 = Tguess[1];
//    		double localB2 = Tguess[1];
//    		
//    		double NMAX = ITERATIONS; /* maximum number of iterations */
//    		double c1 = 0; /* estimated root */
//    		double c2 = 0; /* estimated root */
//    		int index = 0; /* index */
////    		int stuckCount = 0;
//    		 
//    		c1 = (localA1 + localB1)/2.0; /* calculate the midpoint */
//    		c2 = (localA2 + localB2)/2.0; /* calculate the midpoint */
//    		 /* Evaluate loop until the result is less than the tolerance
//    		 * maximum number of iterations is not yet reached*/
//    		
////    		 double[] parametersC = new double[]{D_theta ,psi, K, dt, F0, c};
////    		 double[] parametersA = new double[]{D_theta ,psi, K, dt, F0, localA};
//    		
//    		double[] parametersA1 = new double[]{localA1,localB1};
//    		double[] parametersC1 = new double[]{localA1,c1};
//    		
//    		 if (calcExpression(parametersC,data) == 0)
//    		 {
//    		  /* If the first midpoint gives f(c) = 0, c is the root */
//    	//		 System.out.println("root is " + c);
//    		 }
//    		 else
//    		 {
//    			 double calEx2PCInitial = Math.abs(calcExpression(parametersC,data));
//    			 while ( (calEx2PCInitial > TOL || Double.isNaN(calEx2PCInitial) ) && (index<=NMAX))
//    			 {
//    				 double cex2PC = fsolve.sign(calcExpression(parametersC,data));
//    				 double cex2PA = fsolve.sign(calcExpression(parametersA,data));
//    				 if (cex2PC == cex2PA)
//    				 {
//    					 /* f(c) has same sign as f(a) */
//    					 localA = c;
//    				 }
//    				 else
//    				 {
//    					 /* f(c) has same sign as f(b) */
//    					 localB = c;
//    				 } 
//    				 c = (localA+localB)/2.0; /* midpoint update */
//    				 double difference = Math.abs(localA-localB);
//    				 if (DEBUG)
//    				 {
//    					 System.out.println( (index+1) + " | " + cex2PC + " | " + cex2PA + " | " + localA + " | " + localB + " | " + difference + " | " + c );
//    				 }
//    				
// 
////    				 parametersA = new double[]{D_theta ,psi, K, dt, F0, localA};
////    				 parametersC = new double[]{D_theta ,psi, K, dt, F0, c};
//    		    	parametersA = new double[]{localA,localB};
//    		    	parametersC = new double[]{localA,c};
//    	
//    				 index++; /* index increment */
//    				 calEx2PCInitial = Math.abs(calcExpression(parametersC,data));
//    			 }
//    		 } 
//    		 
//    		 if (index>=NMAX)
//    		 {
//    			 if (DEBUG)
//    			 {
//    				 System.out.println("Root not found " + c + " after " + index + " iterations");	
//    			 }
//    			 
//    			 return ERROR_RETURN;
//    		 }
//    		 
//    		 /* Display results */
//    		 if (DEBUG)
//    		 {
//    			 System.out.println("Root is " + c + " found after " + index + " iterations");		 
//    		 }
//    		
//    		 return new double[]{localA,c};
//    	}
        
        public double getTemperature(double depth)
        {
        	double T=0;
            if (depth > (roof.depth + sub.depth) )
            {
                T = T_interior;
            }
            else
            {
     
                double Rroof = roof.depth/roof.k/n_roof_nodes/2;
                double Ts=(Rroof*T_interior+Constants.rsi_roof*innerT[innerT.length-1])/(Constants.rsi_roof+Rroof);
//                double interp2 = T_substrate;                
//                double interp3 = Ts;
//                double interp4 = depth;
//                double[] interp1Array = new double[]{interp2,interp3};
                double[] interp1Array = new double[innerT.length+2];
                for (int i=0;i<innerT.length;i++)
                {
                	interp1Array[i+1]=innerT[i];
                }
                interp1Array[0]=T_substrate;
                interp1Array[interp1Array.length-1]=Ts;
//                double[] interp2Array = new double[]{depth};
                double[] interReturn = GreenRoofCommon.interpLinear(samples, interp1Array, new double[]{depth});  
                T = interReturn[0];
//                T = interp1([samples],[T_substrate innerT' Ts],depth);  
                                       
                                       
            }
            return T;
        }
        
        public double getVWC(double depth)
        {
        	double vwc=0;
            double dx = (sub.depth)/n_sub_nodes;
                        
            if (depth > sub.depth)
            {
                vwc = 0;
            }
            else
            {
            	 //this is to deal with the colons, i.e. 1:2:10 -> 1 3 5 7 9
            	//  dx/2:dx:(sub.depth-dx/2)
            	// start = dx/2   //  0.015000
            	// interval = dx   // 0.030000
            	// end = sub.depth-dx/2  //  0.13500
            	double[] interp1Array = colonArray(dx/2., dx, (sub.depth-dx/2.));
            	double[] interp2Array = VWC;
            	double[] interReturn = GreenRoofCommon.interpLinear(interp1Array, interp2Array, new double[]{depth}); 
            	vwc = interReturn[0];
 
//                vwc = interp1([dx/2:dx:(sub.depth-dx/2)],[VWC],depth);                
            }
            return vwc;
        }
        
        public double[] colonArray(double start, double interval, double end)
        {	
        	ArrayList<Double> samplesArray = new  ArrayList<Double>();
        	
            
            double startLoop = start;
            double loopValue = startLoop;
            double endLoop = end;
//            int count = 1;
            //this is to deal with the colons, i.e. 1:2:10 -> 1 3 5 7 9
            while (true)
            {
            	samplesArray.add(loopValue);
            	loopValue+=interval;
            	//floating point precision, round off
               	if (Math.abs(loopValue-endLoop) < 0.000001)
            	{
            		continue;
            	}
            	if (loopValue > endLoop)
            	{
            		break;
            	}
            }
        	
        	double[] returnArray = new double[samplesArray.size()];
        	for (int i=0;i<samplesArray.size();i++)
        	{
        		returnArray[i]=samplesArray.get(i);
        	}
        	return returnArray;	
        }
        
     // end methods section
    
 // end class
}
