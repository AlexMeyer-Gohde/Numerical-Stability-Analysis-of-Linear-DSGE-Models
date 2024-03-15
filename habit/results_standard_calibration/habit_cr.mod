
var c           
    k           
    z          
 ;

varexo e        
        ;


parameters   
    h
beta
delta
alpha
sigma
rho
K_aminus1
CK
omega
    ;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
 load('Y')
 Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
 %Y=Y_1;
h=Y(1);
beta=Y(2);
delta=Y(3);
alpha=Y(4);
sigma=Y(5);
rho=Y(6);
omega=Y(7);
K_aminus1=(1/beta-1+delta)/alpha;
CK=K_aminus1-delta;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model(linear);
0=-sigma/(1-h)*c(+1)+sigma*(1+h)/(1-h)*c+ alpha*beta*K_aminus1*(alpha-1)*k-sigma*h/(1-h)*c(-1)+beta*alpha*K_aminus1*rho*z;
0=-CK*c -k+(1-delta+alpha*K_aminus1)*k(-1)+K_aminus1*z;
z=rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

steady_state_model;
 c=0;
k=0;
z=0;
end;



shocks;
var e = omega^2; 
end;

steady;
check;
stoch_simul(order = 1,irf=0,noprint,dr=cycle_reduction,dr_cycle_reduction_tol=1e-16); 

BETA_ZS_dynare_cr=oo_.dr.ghu(oo_.dr.inv_order_var,:);BETA_ZS_dynare_cr=BETA_ZS_dynare_cr(1:2,:);
ALPHA_ZS_dynare_cr=oo_.dr.ghx(oo_.dr.inv_order_var,oo_.dr.inv_order_var);ALPHA_ZS_dynare_cr=ALPHA_ZS_dynare_cr(1:2,1:2);

rp_dynare_cr=4*100*sigma/(1-h)*BETA_ZS_dynare_cr(1)*alpha*K_aminus1*beta*(omega)^2
AA=[ALPHA_ZS_dynare_cr(1,1)-1  ALPHA_ZS_dynare_cr(1,2) BETA_ZS_dynare_cr(1)];
BB=[ALPHA_ZS_dynare_cr BETA_ZS_dynare_cr;0 0 rho ];
DD=[0;0;1];
log_c_std_dynare_cr=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)