package greenroof;

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
        public double samples;
    
    
   public TabaresThermalMass(Plant plant, Substrate sub, Roof roof, int n_substrate_layers,int n_roof_layers, double Area, double dt)
   {
            
            dt = dt*60;
            
            this.roof = roof;
            this.plant = plant;
            this.sub = sub;
            
            //% Initialize
            T_plants = 285; // K
            T_substrate = 291; // K
            T_interior = 293;
            
            
            //% Other      
            n_sub_nodes = n_substrate_layers;
            n_roof_nodes = n_roof_layers;
            L1 = Math.sqrt(Area);            
            
            double dxSubstrate = sub.depth/n_sub_nodes;
            double dxSupport = roof.depth/n_roof_nodes;
            
            
            samples = [0 dxSubstrate/2:dxSubstrate:(sub.depth-dxSubstrate/2) (sub.depth+dxSupport/2):dxSupport:(roof.depth+sub.depth-dxSupport/2) sub.depth+roof.depth];            
            
            //% Initialize inner temperatures by interpolating
            innerT = interp1([0 (sub.depth+roof.depth)],[T_substrate T_interior],samples(2:end-1))';           
            
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
            K[Constants.ONE][Constants.ONE] = 2/rSub;
            
            // Connect within the substate   
            for (int i=Constants.ONE;i<n_sub_nodes-1;i++)
            {
                mcSub = (sub.depth/n_sub_nodes)*sub.rhoCp(VWC[n]);          
                rSub1 = sub.depth/(sub.thermalConductivity(VWC[n])*n_sub_nodes);
                rSub2 = sub.depth/(sub.thermalConductivity(VWC[n+1])*n_sub_nodes);
            
                C[n][n] = mcSub;
                K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]./(rSub1/2+rSub2/2);
                n = n+1;
            }
            
            // Connect interface between both
            rSub = sub.depth/(sub.thermalConductivity(VWC[n])*n_sub_nodes);
            C[n][n] = mcSub;
            K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]./(rSub/2 + rRoof/2);
            n=n+1;
            
            // Connect within the roof            
            for (int i=1;i<n_roof_nodes-1;i++)
            {
                C[n][n] = mcRoof;
                K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+[1,-1;-1,1]./rRoof;
                n = n+1;
            }
            
            // Connect final one.
            C[n][n] = mcRoof;
            consts = Constants;
            K[n][n] = K[n][n]+1/(Constants.rsi_roof + rRoof/2);
   }
        
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
            double Tfilm=  (Tair+T[Constants.ONE])/2;
            double Beta=   1/Tfilm;
            
            
            // OTHER            
            // From table 5.
            double kpor  = sub.phi*air.k+(1-sub.phi)*plant.k;            
            // Auxiliar variable
            double Mg=     VWC[Constants.ONE]/sub.VWCsat;        
            // No idea where this comes from.
            double alphapor= kpor/(dens_a*air.Cp);
            
            // eq. 6 Tabares ... constants may vary?
            double rsub=   34.52*Math.pow(Mg,(-3.2678));
            
            
            //   DIMENSIONLESS NUMBERS
            // Reynolds **note* with L1 not L
            double Re=     dens_a*U*L1/Mu;
            // Prandtl
            double Pr=     air.Cp*Mu/air.k;
            // Grashof
            double Gr=     Math.abs(Constants.g*Beta*Math.pow(dens_a,2)*(T[Constants.ONE]-Tair)*(Math.pow(L1,3))/(Math.pow(Mu,2)));
            // Raleigh
            double Ra=     Gr*Pr;
            // Lewis
            double Le=     1;
            
            // Nusselt. eq. 4 
            double Nu;
            if (Gr<(0.068*Math.pow(Re,2.2)))
            {
                Nu= 3+1.25*0.025*Math.pow(Re,0.8);
            }
            else if ((Gr>(0.068*Math.pow(Re,2.2)))&&(Gr<(55.3*Math.pow(Re,(5/3)))))
            {
                Nu= 2.7*(Math.pow((Gr/(Math.pow(Re,2.2))),(1/3)))*(3*(15/4)+(15/16)*0.0253*Math.pow(Re,0.8));
            }
            else
            {
                Nu= 0.15*(Math.pow(Ra,(1/3)));
            }
            
            // available in Table 5. Tabares
            double Pe=     0.3*L1*U/alphapor; // CHECK
            double hpor  = kpor*1.128*Math.pow(Pe,0.5)/L1; // CHECK            
            double hconv= 15*Nu*air.k/L1; // plantBB was 1.5 then 15
            double ras=  dens_a*air.Cp*(Math.pow(Le,(2/3)))/hconv;                   
            double ra=dens_a*air.Cp*(Math.pow(Le,(2/3)))/hconv;  //ORIGINAL
            ra = ra;
            double hsub =  hpor*hconv/(hpor+hconv); // CHECK
            
            // Vapor Pressures [kPa]
            
            double e_s  = 610.8*Math.exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            double e_sf = 610.8*Math.exp(17.27*(T[Constants.ONE]-273.15)/(T[Constants.ONE]-273.15+237.3))/1000;
            double e_ss = 610.8*Math.exp(17.27*(T[Constants.TWO]-273.15)/(T[Constants.TWO]-273.15+237.3))/1000;   
            double e_air= e_s*RH;

                                                          
            // SHORT WAVE RADIATION
                
            // eq. 12 Tabares
            plant_absorbed_solar = (1-plant.foliage_rho()-plant.tau_fsol())*(1+plant.tau_fsol()*sub.rho)*R_sh;
            
            // eq. 13 Tabares
            substrate_solar_radiation = plant.tau_fsol()*(1-sub.rho)*R_sh;
            
            // LONG WAVE RADIATION
         
            // eq. 14 Tabares
            plant_absorbed_ir_sky = (1-plant.tau_fir())*plant.em*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(Tsky,4));
            
            // eq. 15 Tabares
            substrate_infrared_radiation =   -(plant.tau_fir())*sub.em*Constants.SB*((Math.pow(T[Constants.TWO],4)-Math.pow(Tsky,4)));
            
            // eq. 17 Tabares
            double em_1 = (1/sub.em)+(1/plant.em)-1;
            Qir = (1-plant.tau_fir())*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(T[Constants.TWO],4))/em_1;
            
            // CONVECTION                            
            // eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            plant_convection=   -1.5* plant.LAI*hconv*(T[Constants.ONE]-Tair);
            
            // eq. 19 Tabares
            substrate_convection=    -hsub*(T[Constants.TWO]-Tair);
            
                
            // EVAPOTRANSPIRATION                               
        
            // eq. 22 Tabares
            double f_sol=  1+Math.exp(-0.034*(R_sh-3.5));
            
            // eq. 24 Tabares... extended?
            //f_VPD=  1/(1-0.41*log(e_sf-e_air));
            double VPD=e_ss-e_air;
            double VPD_f=e_sf-e_air;
            double f_VPD;

            if (VPD_f > 0)
            {
                f_VPD=  (1-0.41*Math.log(VPD_f));
            }
            else
            {
                f_VPD= 1;
            }
            if (f_VPD>1)
            {
                f_VPD=1;
            }
            else if (f_VPD<0)
            {
                f_VPD=0.05;
            }

            // eq. 25 Tabares
            double f_temp= 1/(1-0.0016* Math.pow((35-(T[Constants.ONE]-273.15)),2) );
            if (f_temp < 0)
            {
                f_temp = 10000;
            }
            // eq. 23 Tabares
            double rootW = getVWC(sub.depth/2);
            double f_W;
            if ( rootW >0.7*sub.VWCfc) // W_fc vs W_sat
            {
                f_W=    1;
            }
            else if ((rootW <= 0.7*sub.VWCfc) && (rootW > sub.VWCresidual))
            {
                f_W=    (0.7*sub.VWCfc-sub.VWCresidual)/(rootW-sub.VWCresidual);
            }
            else
            {
                f_W=    1000;
            }
            
            
            // eq. 21 Tabares
            //rs = plant.rsmin*f_sol*f_VPD*f_W*f_temp/plant.LAI;
            double f_hum=1/(f_VPD);    
            rs = (plant.rsmin/plant.LAI)*f_sol*f_hum*f_W*f_temp;
            rs = rs;
            fsolar=f_sol;
            fvpd=f_hum;
            fvwc=f_W;
            ftemp=f_temp;
            
            // eq. 20 Tabares
            //Qt= plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            transpiration= -plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ra)); 
            
            // eq. 5 Tabares
            //evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            evaporation= -dens_a*air.Cp*VPD/(air.gamma*(rsub+ras));                               
        
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
            f[Constants.ONE]   = 2*T_substrate/rSub;
            f[f.length-1] = T_interior/(Constants.rsi_roof+rRoof/2);
            //innerT = inv(C/dt + K)*(C/dt*innerT + f);
            //innerT = innerT + C\(f - K*innerT)*dt;
            int total_nodes = n_roof_nodes + n_sub_nodes;
            R = eye(total_nodes,total_nodes)+(C\(dt*K)); 
            innerT  = R\(innerT+(C\f*dt)); 
            
            double Qcond = 2*(T[Constants.TWO] - innerT[Constants.ONE])/rSub;              
            substrate_conduction = -Qcond;
            
            // Update interface_heat_flux            
            rSub = (sub.depth/sub.thermalConductivity(VWC(n_sub_nodes)))/n_sub_nodes; // Last node
            double deltaT = innerT(n_sub_nodes)-innerT(n_sub_nodes+1); //one node on each.            
            interface_heat_flux = 2*deltaT/(rSub+rRoof);
            
            double Tlosa = getTemperature(sub.depth + roof.depth);
            heating_load = ( T_interior - Tlosa)/Constants.rsi_roof;
            
            
        }
        
        public void ResFUN(double[] T, WeatherDataLine data)
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
            double Beta=   1/Tfilm;
            
            
            // OTHER            
            // From table 5.
            double kpor  = sub.phi*air.k+(1-sub.phi)*plant.k;            
            // Auxiliar variable
            double Mg=     VWC[Constants.ONE]/sub.VWCsat;        
            // No idea where this comes from.
            double alphapor= kpor/(dens_a*air.Cp);
            
            // eq. 6 Tabares ... constants may vary?
            double rsub=   34.52*Math.pow(Mg,(-3.2678));
            
            
            //   DIMENSIONLESS NUMBERS
            // Reynolds **note* with L1 not L
            double Re=     dens_a*U*L1/Mu;
            // Prandtl
            double Pr=     air.Cp*Mu/air.k;
            // Grashof
            double Gr=     Math.abs(Constants.g*Beta*Math.pow(dens_a,2)*(T[Constants.ONE]-Tair)*(Math.pow(L1,3))/(Math.pow(Mu,2)));
            // Raleigh
            double Ra=     Gr*Pr;
            // Lewis
            double Le=     1;
            double Nu;
            
            // Nusselt. eq. 4 
            if (Gr<(0.068*Math.pow(Re,2.2)))
            {
                Nu= 3+1.25*0.025*Math.pow(Re,0.8);
            }
            else if ((Gr>(0.068*Math.pow(Re,2.2)))&&(Gr<(55.3*Math.pow(Re,(5/3)))))
            {
                Nu= 2.7*( Math.pow((Gr/(Math.pow(Re,2.2))),(1/3)) )*(3*(15/4)+(15/16)*0.0253*Math.pow(Re,0.8));
            }
            else
            {
                Nu= 0.15*(Math.pow(Re,(1/3)));
            }
            
            // available in Table 5. Tabares
            double Pe=     0.3*L1*U/alphapor; // CHECK
            double hpor  = kpor*1.128*Math.pow(Pe,0.5)/L1; // CHECK            
            double hconv= 15*Nu*air.k/L1; // plantBB was 1.5 then 15
            double ras=  dens_a*air.Cp*( Math.pow(Le,(2/3)) )/hconv;                   
            double ra=dens_a*air.Cp*( Math.pow(Le,(2/3)) )/hconv;  //ORIGINAL
            double hsub =  hpor*hconv/(hpor+hconv); // CHECK
            
            // Vapor Pressures [kPa]
            
            double e_s  = 610.8*Math.exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            double e_sf = 610.8*Math.exp(17.27*(T[Constants.ONE]-273.15)/(T[Constants.ONE]-273.15+237.3))/1000;
            double e_ss = 610.8*Math.exp(17.27*(T[Constants.TWO]-273.15)/(T[Constants.TWO]-273.15+237.3))/1000;   
            double e_air= e_s*RH;

                                                          
            // SHORT WAVE RADIATION
            
            // eq. 12 Tabares
            double Rsh_f = (1-plant.foliage_rho()-plant.tau_fsol())*(1+plant.tau_fsol()*sub.rho)*R_sh;
            
            // eq. 13 Tabares
            double Rsh_s = plant.tau_fsol()*(1-sub.rho)*R_sh;
            
            // LONG WAVE RADIATION
            
            // eq. 14 Tabares
            double Qir_f= (1-plant.tau_fir())*plant.em*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(Tsky,4));
            
            // eq. 15 Tabares
            double Qir_scov=   (plant.tau_fir())*sub.em*Constants.SB*((Math.pow(T[Constants.TWO],4)-Math.pow(Tsky,4)));
            
            // eq. 17 Tabares
            double em_1 = (1/sub.em)+(1/plant.em)-1;
            double Qir_sp = (1-plant.tau_fir())*Constants.SB*(Math.pow(T[Constants.ONE],4)-Math.pow(T[Constants.TWO],4))/em_1;
            
            // CONVECTION
            
            // eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            plant_convection=   1.5* plant.LAI*hconv*(T[Constants.ONE]-Tair);
            
            // eq. 19 Tabares
            substrate_convection=    hsub*(T[Constants.TWO]-Tair);
            
            // EVAPOTRANSPIRATION
            
            // eq. 22 Tabares
            double f_sol=  1+Math.exp(-0.034*(R_sh-3.5));
            
            // eq. 24 Tabares... extended?
            //f_VPD=  1/(1-0.41*log(e_sf-e_air));
            double VPD=e_ss-e_air;
            double VPD_f=e_sf-e_air;
            double f_VPD;

            if (VPD_f > 0)
            {
                f_VPD=  (1-0.41*Math.log(VPD_f));
            }
            else
            {
                f_VPD= 1;
            }
            if (f_VPD>1)
            {
                f_VPD=1;
            }
            else if (f_VPD<0)
            {
                f_VPD=0.05;
            }

            // eq. 25 Tabares
            double f_temp= 1/(1-0.0016* Math.pow((35-(T[Constants.ONE]-273.15)),2) );
            if (f_temp < 0)
            {
                f_temp = 10000;
            }
            
            // eq. 23 Tabares
            double rootW = getVWC(sub.depth/2);
            double f_W;
            if (rootW > 0.7*sub.VWCfc) // W_fc vs W_sat
            {
                f_W=    1;
            }
            else if ((rootW < 0.7*sub.VWCfc) && (rootW > sub.VWCresidual))
            {
                f_W=    (0.7*sub.VWCfc-sub.VWCresidual)/(rootW-sub.VWCresidual);
            }
            else
            {
                f_W=    1000;
            }
            
            
            // eq. 21 Tabares
            //rs = plant.rsmin*f_sol*f_VPD*f_W*f_temp/plant.LAI;
            double f_hum=1/(f_VPD);    
            rs = (plant.rsmin/plant.LAI)*f_sol*f_hum*f_W*f_temp;

            // eq. 20 Tabares
            //Qt= plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            double Qt= plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ra)); 
            
            // eq. 5 Tabares
            //evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            evaporation= dens_a*air.Cp*VPD/(air.gamma*(rsub+ras));    
            
            // CONDUCTION
                                    
            // eq. 8 Tabares modified to only consider the conduction
            // through the first layer of substrate. 
            
            // find temperatures within the subtrate/concrete            
            double rSub = sub.depth/sub.thermalConductivity(VWC[Constants.ONE])/n_sub_nodes; // Using thermal conductivity of the first node.
            double rRoof = roof.depth/roof.k/n_roof_nodes;
            
            double[] f = new double[innerT.length];
            f[Constants.ONE] =  2*T[Constants.TWO]/rSub;
            f[f.length-1]= T_interior/(Constants.rsi_roof+rRoof/2);
           
            //iT = inv(C/dt - K)*((C/dt)*innerT - f);
            //iT = innerT + C\(f - K*innerT)*dt;
            int total_nodes = n_roof_nodes + n_sub_nodes;
            R = eye(total_nodes,total_nodes)+(C\(dt*K)); 
            double iT = R\(innerT+(C\f*dt)); 
            
            double Qcond = 2*(T[Constants.TWO] - iT[Constants.ONE])/rSub;
            
            // ENERGY (& MASS) BALANCE
            
            // Foliage energy balance
            double Ef= - Rsh_f + plant_convection + Qir_f + Qt + Qir_sp; // Eq 14 Camilo
            
            // Substrate energy balance
            double Es= - Rsh_s - Qir_sp + substrate_convection + Qir_scov + evaporation + Qcond; // Eq 15 Camilo
                    
           
            // RETURN
            Res=[Ef Es] ;
            
        } // end resFUN
        
        
        public void moveForward(WeatherDataLine data)
        {
            setMatrices();
            
            options=optimset('Display','off','TolFun',1e-6);
            Tguess = [T_plants T_substrate];                                                      
            Tss = (fsolve(@(TT) ResFUN(TT,data),Tguess,options));
            T_plants=(Tss[Constants.ONE]);
            T_substrate=(Tss[Constants.TWO]); 
            update(data);   
            
        }
        
        public double getTemperature(double depth)
        {
        	double T;
            if (depth > (roof.depth + sub.depth) )
            {
                T = T_interior;
            }
            else
            {
                double Rroof = roof.depth/roof.k/n_roof_nodes/2;
                double Ts=(Rroof*T_interior+Constants.rsi_roof*innerT(end))/(Constants.rsi_roof+Rroof);
                T = interp1([samples],[T_substrate innerT' Ts],depth);                
            }
            return T;
        }
        
        public double getVWC(double depth)
        {
        	double vwc;
            double dx = (sub.depth)/n_sub_nodes;
                        
            if (depth > sub.depth)
            {
                vwc = 0;
            }
            else
            {
                vwc = interp1([dx/2:dx:(sub.depth-dx/2)],[VWC],depth);                
            }
            return vwc;
        }
        
     // end methods section
    
 // end class
}
