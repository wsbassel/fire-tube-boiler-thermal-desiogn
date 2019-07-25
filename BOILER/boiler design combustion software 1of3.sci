function [comb_prod_mass, a0_mix, a1_mix, a2_mix, a3_mix, T_ad_comb, con_eta, exp_eta, con_lamda, exp_lamda  ]=combustion(xc, xh,excess, lc_val )
                          
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
   //  SOFTWARE DONE IN  6 JUNE 2019 
   // boiler design combustion software 1of3
    
    //+++++++++++++++++++++++++++++++++++++++++++++++=
    
    n2_o2_ratio=78.084/20.947
    ar_o2_ratio=0.934/20.947
    mol_o2=32
    mol_n2=28
    mol_ar=40
    mol_co2=44
    mol_h2o=18
    
    mol_o2_stoich=xc/12+(xh/4)
   // mol_air_stoich=(xc/12+xh/4)MOL O2 +n2_o2_ratio * (xc/12+xh/4)MOLN@+ar_o2_ratio * (xc/12+xh/4) MOILAR // o2 , n2 , argon
    mol_air_stoich=(xc/12+xh/4) +n2_o2_ratio * (xc/12+xh/4)+ar_o2_ratio * (xc/12+xh/4)  // o2 , n2 , argon
    mass_air_stoich=(xc/12+xh/4)*mol_o2 +n2_o2_ratio * (xc/12+xh/4)*mol_n2+ar_o2_ratio * (xc/12+xh/4)*mol_ar  // o2 , n2 , argon
    mass_air_excess=mass_air_stoich*(1+excess)
    //disp mass_air
    
comb_prod_mol= xc/12+xh/2+(xc/12+xh/4)*excess+ ...
                          n2_o2_ratio * (xc/12+xh/4)*(1+excess)+ ar_o2_ratio * (xc/12+xh/4)*(1+excess)

//disp (comb_prod_mol)
comb_prod_mass= xc/12*mol_co2+xh/2*mol_h2o+(xc/12+xh/4)*mol_o2*excess+ ...
                          n2_o2_ratio * (xc/12+xh/4)*mol_n2*(1+excess)+ ar_o2_ratio * (xc/12+xh/4)*mol_ar*(1+excess)
                          
                    //                  COMBUSTION PRODUCTS    
                    mol_frac_co2=xc/12 / comb_prod_mol   
                    mol_frac_h2o=xh/2/comb_prod_mol 
                    mol_frac_o2=(xc/12+xh/4)*excess/comb_prod_mol 
                    mol_frac_n2=n2_o2_ratio * (xc/12+xh/4)*(1+excess)/comb_prod_mol 
                    mol_frac_ar=ar_o2_ratio * (xc/12+xh/4)*(1+excess)/comb_prod_mol 
                    
      
      sum_mol_frac=mol_frac_co2+mol_frac_h2o+mol_frac_o2+ mol_frac_n2+mol_frac_ar
      disp mol_frac_co2
      disp (mol_frac_co2)
      disp (mol_frac_h2o)
      
                 mass_frac_co2= xc/12*mol_co2/comb_prod_mass
                 mass_frac_h2o=xh/2*mol_h2o/comb_prod_mass
                 mass_frac_o2=(xc/12+xh/4)*mol_o2*excess/comb_prod_mass
                 mass_frac_n2=n2_o2_ratio * (xc/12+xh/4)*mol_n2*(1+excess)/comb_prod_mass
                 mass_frac_ar=ar_o2_ratio * (xc/12+xh/4)*mol_ar*(1+excess)/comb_prod_mass
                 
      sum_mass_frac=mass_frac_co2+mass_frac_h2o+mass_frac_o2+ mass_frac_n2+mass_frac_ar
      
a0_h2o= 1.8038844
a1_h2o=  0.0053782
a2_h2o= 0.0061945
a3_h2o=-0.0001961

//cp_h2o=a0_h2o+a1_h2o*T*1e-2+a2_h2o*T^2*1e-4+a3_h2o*T^3*1e-6


 
  a0_n2=0.9831813
  a1_n2=  0.0183156
  a2_n2= 0.0000475
  a3_n2=  -0.0000099


  
  //cp_n2=a0_n2+a1_n2*T*1e-2+a2_n2*T^2*1e-4+a3_n2*T^3*1e-6

a0_co2= 0.4969898
a1_co2=   0.1385139
a2_co2=  -0.0081433
a3_co2=   0.0001714

 
   
  // cp_co2=a0_co2+a1_co2*T*1e-2+a2_co2*T^2*1e-4+a3_co2*T^3*1e-6


   
   
  a0_o2=0.8339398
  a1_o2=0.0329738
  a2_o2=-0.00079
  a3_o2=-0.0000006
//cp_o2=a0_o2+a1_o2*T*1e-2+a2_o2*T^2*1e-4+a3_o2*T^3*1e-6

a0_ar=.523
Tref=273.15
a0_mix=a0_h2o*mass_frac_h2o+a0_co2*mass_frac_co2+a0_n2*mass_frac_n2+a0_o2*mass_frac_o2+a0_ar*mass_frac_ar
a1_mix=a1_h2o*mass_frac_h2o+a1_co2*mass_frac_co2+a1_n2*mass_frac_n2+a1_o2*mass_frac_o2
a2_mix=a2_h2o*mass_frac_h2o+a2_co2*mass_frac_co2+a2_n2*mass_frac_n2+a2_o2*mass_frac_o2
a3_mix=a3_h2o*mass_frac_h2o+a3_co2*mass_frac_co2+a3_n2*mass_frac_n2+a3_o2*mass_frac_o2

h_mix_ref=a0_mix*Tref+a1_mix*Tref^2/2*1e-2+a2_mix*Tref^3/3*1e-4+a3_mix*Tref^4/4*1e-6

S_mix_ref=a0_mix*log(Tref)+a1_mix*Tref*1e-2+a2_mix*Tref^2/2*1e-4+a3_mix*Tref^3/3*1e-6


T=1500

T1=T
iter=1
f=a0_mix*T+a1_mix*T^2/2*1e-2+a2_mix*T^3/3*1e-4+a3_mix*T^4/4*1e-6-h_mix_ref-LC_val/comb_prod_mass
//for iter=1:100
while abs(f) >1e-6
dfdt=a0_mix+a1_mix*T*1e-2+a2_mix*T^2*1e-4+a3_mix*T^3*1e-6
dt=f/dfdt
T=T-dt

f=a0_mix*T+a1_mix*T^2/2*1e-2+a2_mix*T^3/3*1e-4+a3_mix*T^4/4*1e-6-h_mix_ref-LC_val/comb_prod_mass
iter=iter+1
end
T_ad_comb=T
disp convergence 
disp (f)

T=1500

for N =1:2
// viscosoty  ln eta=a*ln T+B/T+C/T^2+D
X_mol_n2=mol_frac_n2*sqrt(28)
X_mol_o2=mol_frac_o2*sqrt(32)
X_mol_co2=mol_frac_co2*sqrt(44)
X_mol_h2o=mol_frac_h2o*sqrt(18)
sum_X_mol=X_mol_n2+X_mol_o2+X_mol_co2+X_mol_h2o
eta_exp_n2=0.6124009

 lamda_exp_n2=0.6495569

 con_eta_n2=0.5789397

 con_amda_n2=0.0007182
 eta_n2=con_eta_n2*T^eta_exp_n2
 lamda_n2=con_amda_n2*T^lamda_exp_n2
 eta_exp_o2=0.6795293

 lamda_exp_o2=0.8313437

 con_eta_o2=0.4380801

 con_lamda_o2=0.000238
 eta_o2=con_eta_o2*T^eta_exp_o2
 lamda_o2=con_lamda_o2*T^lamda_exp_o2
 
 eta_exp_co2=0.8432338

 lamda_exp_co2=1.3797723

 con_eta_co2= 0.1218704

 con_lamda_co2=0.0000063
eta_co2=con_eta_co2*T^eta_exp_co2
 lamda_co2=con_lamda_co2*T^lamda_exp_co2
 
 eta_exp_h2o=1.0541399

 lamda_exp_h2o=1.1821355

 con_eta_h2o=0.0242389

 con_lamda_h2o=0.0000219
 
 eta_h2o=con_eta_h2o*T^eta_exp_h2o
 lamda_h2o=con_lamda_h2o*T^lamda_exp_h2o
 
 
 eta_mix(N)=(eta_n2*X_mol_n2+eta_o2*X_mol_o2+eta_co2*X_mol_co2+eta_h2o*X_mol_h2o)/sum_X_mol
 
 lamda_mix(N)=(lamda_n2*X_mol_n2+lamda_o2*X_mol_o2+lamda_co2*X_mol_co2+lamda_h2o*X_mol_h2o)/sum_X_mol
 
 T=300
 T2=T
end

 LT1=log(T1)
 LT2=log(T2)
 A1=zeros(2,2)
 A1=[1 LT1
    1  LT2]
   
    Leta1=log(eta_mix(1))
    Leta2=log(eta_mix(2))
 B1=[Leta1
     Leta2]
 C1=A1\B1
 con_eta=exp(C1(1))
 exp_eta=C1(2)
 
 A2=[1 LT1
    1  LT2]
    Llamd1=log(lamda_mix(1))
    Llamd2=log(lamda_mix(2))
 B2=[Llamd1
     Llamd2]
 C2=A2\B2
con_lamda=exp(C2(1)) 
exp_lamda=C2(2)


 
 
 
 

endfunction

//xc=.85
//    xh=.15
//    excess=.15
 //   LC_val=39000
//[comb_prod_mass, a0_mix, a1_mix, a2_mix, a3_mix, T_ad_comb, con_eta, exp_eta, con_lamda, exp_lamda  ]=combustion(xc, xh,excess, LC_val )
//disp ao_mix
//disp (a0_mix)
//disp a1_mix
//disp (a1_mix)
//disp a2_mix
//disp (a2_mix)
//disp a3_mix
//disp (a3_mix)
//disp con_eta
//disp (con_eta)
//disp  exp_eta
//disp(exp_eta)
//disp con_lamda
//disp (con_lamda)
//disp exp_lamda
//disp (exp_lamda)
//disp T_ad_comb
//disp (T_ad_comb)
//disp comb_prod_mass
//disp (comb_prod_mass)
