function [T_ad_comb, Q_fur, T2, Tc, M]=temp_radiation (a0_mix, a1_mix, a2_mix, a3_mix,comb_prod_mass,T_air_in, Q_boiler, LC_val)
     
      
  //T=2000
  H0 =(LC_val+(comb_prod_mass-1)*T_air_in)/ comb_prod_mass

  //Q_boiler=2500
  M=Q_boiler/LC_val*comb_prod_mass
  //disp M
  //disp(M)
  T=2000
 // disp H0
//disp (H0)
   f=a0_mix*T+a1_mix*T^2/2*1e-2+a2_mix*T^3/3*1e-4+a3_mix*T^4/4*1e-6-273.15-H0
 for iter =1:50
   
   df_dt=a0_mix+a1_mix*T*1e-2+a2_mix*T^2*1e-4+a3_mix*T^3*1e-6
   delta=df_dt\f
   T=T-delta
    f=a0_mix*T+a1_mix*T^2/2*1e-2+a2_mix*T^3/3*1e-4+a3_mix*T^4/4*1e-6-273.15-H0
end
//disp T
//disp (T)
 //a0_mix=0.950939
 // a1_mix=0.0381212
 // a2_mix=-0.000923
 // a3_mix=0.0000077
 // con_eta=0.367437
   
 //exp_eta= 0.6788303
  
 //con_lamda=0.000257

//exp_lamda=0.8117109

T_ad_comb =T
eps_g=.07
eps_wall=.8
A_fur=%pi*1.1*4.5
T_sat=160+273.15
Q_fur=250
//Q_aux=(1/(1/eps_g+1/eps_wall-1))*5.67*((1500/100)^4-(T_sat/100)^4)*A_fur*1e-3

Tc=1600
T2=1400
H2=(a0_mix*T2+a1_mix*T2^2/2*1e-2+a2_mix*T2^3/3*1e-4+a3_mix*T2^4/4*1e-6-273.15)

//Q_fur=M*(H0-(a0_mix*T2+a1_mix*T2^2/2*1e-2+a2_mix*T2^3/3*1e-4+a3_mix*T2^4/4*1e-6-273.15) )

f(1)=-Q_fur+(1/(1/eps_g+1/eps_wall-1))*5.67*((Tc/100)^4-(T_sat/100)^4)*A_fur*1e-3
f(2)=-Q_fur+M*(H0-(a0_mix*T2+a1_mix*T2^2/2*1e-2+a2_mix*T2^3/3*1e-4+a3_mix*T2^4/4*1e-6-273.15) )
//f(3)=T_ad_comb-T2-1.2*(T_ac_comb-Tc)
f(3)=-.2*T_ad_comb+1.2*Tc-T2
  //f(3)=-.2*T_flame+1.2*Tc-T2
  for iter=1 :500
 
  //while max(abs(f))> 1e-6
 j=zeros(3,3)
 j(1,1)=-1
 j(1,2)=(1/(1/eps_g+1/eps_wall-1))*5.67*((Tc)^3*4/(100)^4)*A_fur*1e-3
 j(2,1)=-1
 j(2,3)=M*(-(a0_mix*+a1_mix*T2*1e-2+a2_mix*T2^2*1e-4+a3_mix*T2^3*1e-6 ))
 j(3,2)=1.2
 j(2,3)=-1
 
 // j=1/(1/eps_g+1/eps_wall-1)*5.67*(Tc^3*4/(100^4))*A_fur*1e-3
  
  
  //j(3,2)=1.2
 // j(3,3)=-1
  delta=j\f
  nf=10
  Q_fur=Q_fur-delta(1)/nf
  Tc=Tc-delta(2)/nf
  T2=T2-delta(3)/nf
  f(1)=-Q_fur+(1/(1/eps_g+1/eps_wall-1))*5.67*((Tc/100)^4-(T_sat/100)^4)*A_fur*1e-3
f(2)=-Q_fur+M*(H0-(a0_mix*T2+a1_mix*T2^2/2*1e-2+a2_mix*T2^3/3*1e-4+a3_mix*T2^4/4*1e-6-273.15) )
f(3)=-.2*T_ad_comb+1.2*Tc-T2
//disp (f)
end

T_inlet_next=T2
disp iter
disp (iter)
 // Tc=Tc-delta/10
 // T2=T2-delta(3)/10
  //f=-Q_fur+(1/(1/eps_g+1/eps_wall-1))*5.67*((Tc/100)^4-(T_sat/100)^4)*A_fur*1e-3
//f(2)=-Q_fur+M*(H0-(a0_mix*T2+a1_mix*T2^2/2*1e-2+a2_mix*T2^3/3*1e-4+a3_mix*T2^4/4)*1e-6-273.15)
 // f(3)=-.2*T_flame+1.2*Tc-T2
  fmax=max(abs(f))
disp convergence
disp (fmax)


  
  
  
  
endfunction
 //a0_mix=0.950939
 // a1_mix=0.0381212
 // a2_mix=-0.000923
 // a3_mix=0.0000077
 // con_eta=0.367437
   
 //exp_eta= 0.6788303
  
 //con_lamda=0.000257

//exp_lamda=0.8117109
 // T_air_in=30
 // Q_boiler=2500
 // comb_prod_mass=18.22
 // LC_val=37000
  //T_ad_comb=1960
 //[T_ad_comb, Q_fur, T_inlet_next, Tc, M]=temp_radiation (a0_mix, a1_mix, a2_mix, a3_mix,comb_prod_mass,T_air_in, Q_boiler, LC_val)
// disp T_ad_comb
 //   disp (T_ad_comb)
  //  disp M
   // disp (M)
    
//disp Q_fur
//disp (Q_fur)
//disp Tc
//disp (Tc)
//disp T_inlet_next
//disp (T_inlet_next)

    
