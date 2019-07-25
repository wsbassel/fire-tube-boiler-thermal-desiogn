//function [comb_prod_mass]=boiler_thermal_design ()
getd("BOILER")

    xc=.85  //carbon mass per kg of feul 
    xh=.15  // hidrogen mass per kg of feul
    excess=.15    //excess air
    LC_val=37000  // lower calorfic value of feul
[comb_prod_mass, a0_mix, a1_mix, a2_mix, a3_mix, T_ad_comb, con_eta, exp_eta, con_lamda, exp_lamda  ]=combustion(xc, xh,excess, LC_val )
disp 'cp=a0_mix +a1_mix*T+a2_mix*T^2+ a3_mix*T^3' 
a0_mix  // pqrqmeter of cp equation
disp (a0_mix)
a1_mix // pqrqmeter of cp equation
disp (a1_mix)
a2_mix // pqrqmeter of cp equation
disp (a2_mix)
a3_mix // pqrqmeter of cp equation
disp (a3_mix)
disp  'eta=con_eta*T^exp_eta'
disp 'lamd=con_lamda*T^exp_lamda'
con_eta
disp (con_eta)
disp  (exp_eta)
disp (con_lamda)
disp (exp_lamda)
disp 'adiabatic comustion temperature'

disp (T_ad_comb) // adiabatic comustion temperature

disp 'comb_prod_mass per kg of feul'
disp (comb_prod_mass)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                       RADIATION HEAT TRANSFER
    T_air_in=30 // inlet air temperature
  Q_boiler=2500  //boiler power
 
 [T_ad_comb, Q_fur, T_inlet_next, Tc, M]=temp_radiation (a0_mix, a1_mix, a2_mix, a3_mix,comb_prod_mass,T_air_in, Q_boiler, LC_val)
  Q_fur

Tc
T_inlet_next  
    
disp 'Q_fur  furnace radiant heat tyransfer'

disp (Q_fur)

disp 'T_inlet_next is the outlet temperature of furnace which the inlet tempeature of pass2'
disp (T_inlet_next)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//            CONVECTION HEAT TRANSFER pass2

L_t_pass=4.5 // length of tube of the pass2            
  N_t_pass=48    // number of tubes in pass2  
  d_in_t=.05      // internal diameter of tube       
        
    T_sat=160+273.15  // saturation temperature of steam
    n_div=20          // number of divisions
[T_inlet_next,Q_convec, Dp_tot]=fire_tube_conv_rev2(M,a0_mix,a1_mix,a2_mix,a3_mix,...
con_eta,exp_eta,con_lamda,exp_lamda,L_t_pass,N_t_pass,d_in_t,T_inlet_next,T_sat,n_div)
disp 'T_inlet_next the inlet temperature of next pass'
disp (T_inlet_next)
Q_convec_pass(2)=Q_convec
disp 'Q_tot_convec_pass(2) convective heat transfered' 
disp (Q_convec_pass(2))
Dp_tot_pass(2)=Dp_tot
disp 'pressure_loss_Pa of the pass2'
disp (Dp_tot_pass(2))
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//             convection heat transfer pass3
T_inlet_next
L_t_pass=5         // length of tube of the pass3  
  N_t_pass=30     // number of tubes in pass3 
  d_in_t=.05         // internal diameter of tube          
        
    T_sat=160+273.15 // saturation temperature of steam
    n_div=20         // number of divisions
[T_inlet_next,Q_convec, Dp_tot]=fire_tube_conv_rev2(M,a0_mix,a1_mix,a2_mix,a3_mix,...
con_eta,exp_eta,con_lamda,exp_lamda,L_t_pass,N_t_pass,d_in_t,T_inlet_next,T_sat,n_div)
disp T_inlet_next
disp (T_inlet_next)
Q_convec_pass(3)=Q_convec
Dp_tot_pass(3)=Dp_tot
disp Q_tot_convec_pass(3)
disp (Q_convec_pass(3))
disp pressure_loss_Pa
disp (Dp_tot_pass(3))
Dp_total=Dp_tot_pass(2)+Dp_tot_pass(3)

Q_fur_convec=Q_fur+Q_convec_pass(2)+Q_convec_pass(3)
disp Dp_total
disp (Dp_total)
disp heat_transfered 
disp (Q_fur_convec)
