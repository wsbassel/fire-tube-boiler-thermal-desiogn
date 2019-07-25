function [T_inlet_next,Q_convec, Dp_tot]=fire_tube_conv_rev2(M,a0_mix,a1_mix,a2_mix,a3_mix,...
con_eta,exp_eta,con_lamda,exp_lamda,L_t_pass,N_t_pass,d_in_t,T_inlet_next,T_sat,n_div)


     p=1e5
     T2=1564.2287
   //for Iter=1:2
   
   L_t=L_t_pass
    L_div=L_t/n_div
   
   
   
    Asr=L_t/n_div*%pi*d_in_t
    T_ref=273.15
    h_mix_ref=a0_mix*T_ref+a1_mix*T_ref^2/2*1e-2+a2_mix*T_ref^3/3*1e-4+a3_mix*T_ref^4/4*1e-6
    alfa=.07
    
    

   
   // T_fur=666.38448
   //T_in_div=T_1iter
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //                calculations of (n_div-1)  divisions
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   m_t=M/N_t_pass
   n_tube=N_t_pass
   
   // n_div=5
    qx=1
   T_in_div(1)=T_inlet_next
    for n =1 : n_div-1
      //  if n > 1 then
          //  T_in_div(n)=Tx_out
        //    disp (Tx_out)
            //end 
        
   // disp (T_in_div(n))
    h_in(n)=a0_mix*T_in_div(n)+a1_mix*T_in_div(n)^2/2*1e-2+a2_mix*T_in_div(n)^3/3*1e-4+a3_mix*T_in_div(n)^4/4*1e-6-h_mix_ref
    eta(n)=con_eta*T_in_div(n)^exp_eta*1e-6                     
lamd(n)=con_lamda*T_in_div(n)^exp_lamda
 cp(n)=a0_mix+a1_mix*T_in_div(n)*1e-2+a2_mix*T_in_div(n)^2*1e-4+a3_mix*T_in_div(n)^3*1e-6
 
Pr(n)=cp(n)*eta(n)/lamd(n)*1000
mol_mass=29.16
r_j=8314/mol_mass
ro(n)=p/(r_j*T_in_div(n))
v(n)=m_t/ro(n)/(%pi/4*d_in_t^2)

G=m_t/(%pi/4*d_in_t^2)
Re(n)=G*d_in_t/eta(n)

//if Re <= 1e5 then
 //   Dp=.316/Re^.25*L_div/d_in_t*ro*v^2/2
//else
 //   Dp=.184/Re^.2*L_div/d_in_t*ro*v^2/2
//end
 f1_bare=1.98*log10(Re(n)/6.9)
        f_bare=1/f1_bare^2
 Dp(n)=f_bare*L_div/d_in_t*ro(n)*v(n)^2/2

alfa(n)=.023*Re(n)^.8*Pr(n)^.3*lamd(n)/d_in_t*1e-3
Tx_in=T_in_div(n)
T_out_div(n)=T_in_div(n)*.99
hx_in=h_in(n)
//qx=q(n)
//Tx_in=T_in_div(n)
//Tx_out=T_out_div(n)
//T_out_div(n)=Tx_out

Tx_out=T_out_div(n)
//T_in_div(n)=Tx_out
alfax=alfa(n)
 f(1)=-qx+alfax*Asr*(.5*Tx_in+.5*Tx_out-T_sat)
f(2)=-m_t*(a0_mix*Tx_out+a1_mix*Tx_out^2/2*1e-2+a2_mix*Tx_out^3/3*1e-4+a3_mix* ...
 Tx_out^4/4*1e-6-h_mix_ref)+m_t*hx_in-qx

 j=zeros(2,2)
 //for iter=1:30
 while max(abs(f))> 1e-8
 j(1,1)=-1
 j(1,2)=alfax*Asr*.5
 j(2,1)=-1
 j(2,2)=-m_t*(a0_mix+a1_mix*Tx_out*1e-2+a2_mix*Tx_out^2*1e-4+a3_mix*Tx_out^3*1e-6)
 delta=j\f
 nf=1
 qx=qx-delta(1)/nf
 
 Tx_out=Tx_out-delta(2)/nf
 
f(1)=-qx+alfax*Asr*(.5*Tx_in+.5*Tx_out-T_sat)
 f(2)=-m_t*(a0_mix*Tx_out+a1_mix*Tx_out^2/2*1e-2+a2_mix*Tx_out^3/3*1e-4+a3_mix* ...
 Tx_out^4/4*1e-6-h_mix_ref)+m_t*hx_in-qx
 f1_N(n)=f(1)
 f2_N(n)=f(2)
 //disp (j)
end

//disp (Tx_out)
q(n)=qx
T_out_div(n)=Tx_out
T_in_div(n+1)=Tx_out
end

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//            calculation of the (n_div ) division on;y one division


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//q(n)=qx

 h_in(n_div)=a0_mix*T_in_div(n_div)+a1_mix*T_in_div(n_div)^2/2*1e-2+a2_mix*T_in_div(n_div)^3/3*1e-4+ ...
 a3_mix*T_in_div(n_div)^4/4*1e-6-h_mix_ref
Tx_in=T_in_div(n_div)
//Tx_out=T_out_div(n) 
 eta(n_div)=con_eta*T_in_div(n_div)^exp_eta*1e-6                     
lamd(n_div)=con_lamda*T_in_div(n_div)^exp_lamda
 cp(n_div)=a0_mix+a1_mix*T_in_div(n_div)*1e-2+a2_mix*T_in_div(n_div)^2*1e-4+a3_mix*T_in_div(n_div)^3*1e-6
 
Pr(n_div)=cp(n_div)*eta(n_div)/lamd(n_div)*1000
mol_mass=29.16
r_j=8314/mol_mass
ro(n_div)=p/(r_j*T_in_div(n_div))
v(n_div)=m_t/ro(n_div)/(%pi/4*d_in_t^2)

G=m_t/(%pi/4*d_in_t^2)
Re(n_div)=G*d_in_t/eta(n_div)

//if Re <= 1e5 then
 //   Dp=.316/Re^.25*L_div/d_in_t*ro*v^2/2
//else
 //   Dp=.184/Re^.2*L_div/d_in_t*ro*v^2/2
//end
 f1_bare=1.98*log10(Re(n_div)/6.9)
        f_bare=1/f1_bare^2
 Dp(n_div)=f_bare*L_div/d_in_t*ro(n_div)*v(n_div)^2/2

alfa(n_div)=.023*Re(n_div)^.8*Pr(n_div)^.3*lamd(n_div)/d_in_t*1e-3
alfax=alfa(n_div)
hx_in=h_in(n_div)
 f(1)=-qx+alfax*Asr*(.5*Tx_in+.5*Tx_out-T_sat)
 f(2)=-m_t*(a0_mix*Tx_out+a1_mix*Tx_out^2/2*1e-2+a2_mix*Tx_out^3/3*1e-4+a3_mix* ...
 Tx_out^4/4*1e-6-h_mix_ref)+m_t*hx_in-qx
 j=zeros(2,2)
 //for iter=1:30
 while max(abs(f))>1e-8
 j(1,1)=-1
 j(1,2)=alfax*Asr*.5
 j(2,1)=-1
 j(2,2)=-m_t*(a0_mix+a1_mix*Tx_out*1e-2+a2_mix*Tx_out^2*1e-4+a3_mix*Tx_out^3*1e-6)
 delta=j\f
 //disp (delta)
 nf=1
 qx=qx-delta(1)/nf
 
 Tx_out=Tx_out-delta(2)/nf
 
f(1)=-qx+alfax*Asr*(.5*Tx_in+.5*Tx_out-T_sat)
 f(2)=-m_t*(a0_mix*Tx_out+a1_mix*Tx_out^2/2*1e-2+a2_mix*Tx_out^3/3*1e-4+a3_mix* ...
 Tx_out^4/4*1e-6-h_mix_ref)+m_t*hx_in-qx
 f1_N(n_div)=f(1)
 f2_N(n_div)=f(2)
end
disp convergencef1
f1_N_max=max(f1_N)
disp (f1_N_max)
disp convergencef2
f2_N_max=max(f2_N)
disp (f2_N_max)

Dp_tot=sum(Dp)
T_out_div(n_div)=Tx_out
T_inlet_next=Tx_out

q(n_div)=qx

Q_tube=sum (q)
Q_convec=Q_tube*n_tube
//disp Q_TOT
//disp (Q_tot)
//disp exit_temperature
//disp (T_inlet_next)
//disp inlet_temp
//disp (T_in_div)

endfunction
 //M=1.2305611
 //boiler_power=2500
  //  Low_cal_val=37000

    // a0_mix=0.950939
  //a1_mix=0.0381212
  //a2_mix=-0.000923
  //a3_mix=0.0000077
  //con_eta=0.367437
   
 //exp_eta= 0.6788303
  
 //con_lamda=0.000257

//exp_lamda=0.8117109
 // L_t_pass=4.5             
 // N_t_pass=48      
 // d_in_t=.05             
 // T_fur=1490        
  //  T_sat=160+273.15
   // n_div=20
    //T_inlet_next,Q_tot]=fire_tube_conv_rev2(M,a0_mix,a1_mix,a2_mix,a3_mix,...
//con_eta,exp_eta,con_lamda,exp_lamda,L_t_pass,N_t_pass,d_in_t,T_inlet_next,T_sat,n_div )

