dynare_location='C:\dynare\4.6.1\matlab';
 save('dynare_location','dynare_location')

clear all
 load('Y')
 Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
 %Y=Y_1;
syms x_11 x_12 x_21 x_22 t_1 t_2
X=[x_11 x_12; x_21 x_22];
T=[t_1; t_2];
h=sym('h');beta=sym('beta');delta=sym('delta');alpha=sym('alpha');sigma=sym('sigma');rho=sym('rho');omega=sym('omega');
K_aminus1=(1/beta-1+delta)/alpha;
CK=K_aminus1-delta;
A=[-sigma/(1-h) 0; 0 0];
B=[sigma*(1+h)/(1-h) alpha*beta*K_aminus1*(alpha-1);-CK -1];
C=[-sigma*h/(1-h) 0;0 (1-delta+alpha*K_aminus1)];
D=[beta*alpha*K_aminus1*rho; K_aminus1];


R=A*X^2+B*X+C;
s_1=solve(R(1,1), R(2,1), R(1,2), R(2,2), x_11, x_12, x_21, x_22);
X_stable=[s_1.x_11(1) s_1.x_12(1);s_1.x_21(1) s_1.x_22(1)];
X_stable=subs(X_stable,{h beta delta alpha sigma rho},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6)});
Q=(A*X_stable+A*rho+B)*T+D;
z_1=solve(Q(1,1), Q(2,1),  t_1, t_2);
Q_stable=[z_1.t_1; z_1.t_2];

Q_stable=subs(Q_stable,{h beta delta alpha sigma rho},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6)});
rp_exact=4*100*sigma/(1-h)*Q_stable(1)*alpha*K_aminus1*beta*(omega)^2;
rp_exact=subs(rp_exact,{h beta delta alpha sigma rho omega},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6) Y(7)});


AA=[X_stable(1,1)-1  X_stable(1,2) Q_stable(1)];
BB=[X_stable Q_stable;0 0 rho ];
DD=[0;0;1];

log_c_std_exact=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2);
log_c_std_exact=subs(log_c_std_exact,{h beta delta alpha sigma rho omega},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6) Y(7)});

X_exact=vpa(X_stable,100);
Q_exact=vpa(Q_stable,100);
rp_exact=vpa(rp_exact,100)
log_c_std_exact=vpa(log_c_std_exact,100)
eig_s=eig(X_exact)
X_unstable=-inv(B+A*X_exact)*A;
X_unstable=subs(X_unstable,{h beta delta alpha sigma rho},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6)});
X_unstable=vpa(X_unstable,100);
eig_u=eig(X_unstable);
eig_u=1./eig_u
eig_seperation=eig_u-eig_s';[a b]=min(abs(eig_seperation(:)));eig_seperation_exact=eig_seperation(b)
eigval=[eig_s;eig_u];%eigval=vpa(eigval,100);
[sortabs,sortindex] = sort(abs(eigval));
eigval_exact=eigval(sortindex)
A_d=subs(A,{h beta delta alpha sigma rho omega},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6) Y(7)});
B_d=subs(B,{h beta delta alpha sigma rho omega},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6) Y(7)});
C_d=subs(C,{h beta delta alpha sigma rho omega},{Y(1) Y(2) Y(3) Y(4) Y(5) Y(6) Y(7)});
inputs.A=double(A_d);
inputs.B=double(B_d);
inputs.C=double(C_d);
inputs.X=double(X_stable);
[errors_exact] = matrix_quadratic_backward_errors(inputs)
 save('exact_results','rp_exact','X_exact','Q_exact','log_c_std_exact','eig_seperation_exact','eigval_exact','errors_exact')



clear all
 addpath('uhlig')
 load('Y')
 Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
 %Y=Y_1;
example_habit;
ALPHA_ZS_uhlig=PP;
BETA_ZS_uhlig=QQ;
eig_seperation=Xi_sortval(abs(Xi_sortval)>1)-Xi_sortval(abs(Xi_sortval)<1)';[a b]=min(abs(eig_seperation(:)));eig_seperation_uhlig=eig_seperation(b)
eigval_uhlig=Xi_sortval 
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
A=[-sigma/(1-h) 0; 0 0];
B=[sigma*(1+h)/(1-h) alpha*beta*K_aminus1*(alpha-1);-CK -1];
C=[-sigma*h/(1-h) 0;0 (1-delta+alpha*K_aminus1)];
D=[beta*alpha*K_aminus1*rho; K_aminus1];
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_uhlig;
[errors_uhlig] = matrix_quadratic_backward_errors(inputs)
save('uhlig_results','rp_uhlig','ALPHA_ZS_uhlig','BETA_ZS_uhlig','log_c_std_uhlig','eig_seperation_uhlig','eigval_uhlig','errors_uhlig')
rmpath('uhlig')



clear all
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
A=[-sigma/(1-h) 0; 0 0];
B=[sigma*(1+h)/(1-h) alpha*beta*K_aminus1*(alpha-1);-CK -1];
C=[-sigma*h/(1-h) 0;0 (1-delta+alpha*K_aminus1)];
D=[beta*alpha*K_aminus1*rho; K_aminus1];
% [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=QZ_solve(A,B,C,zeros(2,1),D,rho,2,2,1,1);
% max(max(abs(X_1-ALPHA_ZS)))
% abs(T_1-BETA_ZS)
% 100*(BETA_ZS(1)-T_1(1))/T_1(1)
% 
% rp_amg=4*100*sigma/(1-h)*BETA_ZS(1)*alpha*K_aminus1*beta*(0.0071)^2

X.A=A;
X.B=B;
X.C=C;
figure;ps_quad(X,10000,1e-6);
figure;ps_quad(X,10000,2e-5);

[ALPHA_ZS_aim,BETA_ZS_aim,sorted_eigenvalues,existence_uniqueness]=...
    sp_solve(A,B,C,zeros(2,1),D,rho,2,2,1,1);

AA=[ALPHA_ZS_aim(1,1)-1  ALPHA_ZS_aim(1,2) BETA_ZS_aim(1)];
BB=[ALPHA_ZS_aim BETA_ZS_aim;0 0 rho ];
DD=[0;0;1];
log_c_std_aim=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)
rp_aim=4*100*sigma/(1-h)*BETA_ZS_aim(1)*alpha*K_aminus1*beta*(omega)^2
eig_seperation=sorted_eigenvalues(abs(sorted_eigenvalues)>1)-sorted_eigenvalues(abs(sorted_eigenvalues)<1)';[a b]=min(abs(eig_seperation(:)));eig_seperation_aim=eig_seperation(b)
eigval_aim=sorted_eigenvalues(sorted_eigenvalues~=rho);
[sortabs,sortindex] = sort(abs(eigval_aim));
eigval_aim=eigval_aim(sortindex);eigval_aim=[eigval_aim; 1/0]
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_aim;
[errors_aim] = matrix_quadratic_backward_errors(inputs)
save('aim_results','rp_aim','ALPHA_ZS_aim','BETA_ZS_aim','log_c_std_aim','eig_seperation_aim','eigval_aim','errors_aim')


 
 
 [BETA_ZS_solab,ll,ff,ALPHA_ZS_solab]  = solab([ zeros(2,2) -A;  eye(2,2) zeros(2,2)],[C B;   zeros(2,2) eye(2,2)],rho,[D; zeros(2,1)],2);
 AA=[ALPHA_ZS_solab(1,1)-1  ALPHA_ZS_solab(1,2) BETA_ZS_solab(1)];
BB=[ALPHA_ZS_solab BETA_ZS_solab;0 0 rho ];
DD=[0;0;1];
log_c_std_solab=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)
rp_solab=4*100*sigma/(1-h)*BETA_ZS_solab(1)*alpha*K_aminus1*beta*(omega)^2
 eig_s=eig(ALPHA_ZS_solab)
X_unstable=-(B+A*ALPHA_ZS_solab)\A;
eig_u=eig(X_unstable);
eig_u=1./eig_u
eig_seperation=eig_u-eig_s';[a b]=min(abs(eig_seperation(:)));eig_seperation_solab=eig_seperation(b)
 eigval=[eig_s;eig_u];
[sortabs,sortindex] = sort(abs(eigval));
eigval_solab=eigval(sortindex)
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_solab;
[errors_solab] = matrix_quadratic_backward_errors(inputs)
 save('solab_results','rp_solab','ALPHA_ZS_solab','BETA_ZS_solab','log_c_std_solab','eig_seperation_solab','eigval_solab','errors_solab')
 
 
g0=[[A(1,1); 0],  B D; 0 1 0 0; 0 0 0 1];
g1=[ zeros(2,1), -C zeros(2,1); 1 0 0 0; 0 0 0 rho];
c=zeros(4,1);
psi=[zeros(3,1); 1];
pi=zeros(4,4); pi(3,3)=-1;
addpath('gensys')
[G1,Cgensys,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi);
ALPHA_ZS_gensys=G1(2:3,2:3);
BETA_ZS_gensys=impact(2:3,1);

AA=[ALPHA_ZS_gensys(1,1)-1  ALPHA_ZS_gensys(1,2) BETA_ZS_gensys(1)];
BB=[ALPHA_ZS_gensys BETA_ZS_gensys;0 0 rho ];
DD=[0;0;1];
log_c_std_gensys=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)
rp_gensys=4*100*sigma/(1-h)*impact(2)*alpha*K_aminus1*beta*(omega)^2
eig_s=eig(ALPHA_ZS_gensys)
X_unstable=-(B+A*ALPHA_ZS_gensys)\A;
eig_u=eig(X_unstable);
eig_u=1./eig_u
eig_seperation=eig_u-eig_s';[a b]=min(abs(eig_seperation(:)));eig_seperation_gensys=eig_seperation(b)
eigval=[eig_s;eig_u];
[sortabs,sortindex] = sort(abs(eigval));
eigval_gensys=eigval(sortindex)
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_gensys;
[errors_gensys] = matrix_quadratic_backward_errors(inputs)
 save('gensys_results','rp_gensys','ALPHA_ZS_gensys','BETA_ZS_gensys','log_c_std_gensys','eig_seperation_gensys','eigval_gensys','errors_gensys')
rmpath('gensys') 
 
 
[ALPHA_ZS_bp,BETA_ZS_bp]=bp(A,-B,C,D,zeros(2,1),rho,80e5,1e-10);
AA=[ALPHA_ZS_bp(1,1)-1  ALPHA_ZS_bp(1,2) BETA_ZS_bp(1)];
BB=[ALPHA_ZS_bp BETA_ZS_bp;0 0 rho ];
DD=[0;0;1];
log_c_std_bp=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)
rp_bp=4*100*sigma/(1-h)*BETA_ZS_bp(1)*alpha*K_aminus1*beta*(omega)^2
eig_s=eig(ALPHA_ZS_bp)
X_unstable=-(B+A*ALPHA_ZS_bp)\A;
eig_u=eig(X_unstable);
eig_u=1./eig_u
eig_seperation=eig_u-eig_s';[a b]=min(abs(eig_seperation(:)));eig_seperation_bp=eig_seperation(b)
eigval=[eig_s;eig_u];
[sortabs,sortindex] = sort(abs(eigval));
eigval_bp=eigval(sortindex)
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_bp;
[errors_bp] = matrix_quadratic_backward_errors(inputs)
 save('bp_results','rp_bp','ALPHA_ZS_bp','BETA_ZS_bp','log_c_std_bp','eig_seperation_bp','eigval_bp', 'errors_bp')
 
 
 
 clear all
 p = path; 
 load('dynare_location')
 addpath(dynare_location)
 dynare habit noclearall;
 eigenvalues=oo_.dr.eigval(abs(oo_.dr.eigval-rho)>1e-14);
 eig_seperation=eigenvalues(abs(eigenvalues)>1)-eigenvalues(abs(eigenvalues)<1)';[a b]=min(abs(eig_seperation(:)));eig_seperation_dynare=eig_seperation(b)
[sortabs,sortindex] = sort(abs(eigenvalues));
eigval_dynare=eigenvalues(sortindex); eigval_dynare=[eigval_dynare;1/0]
 Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
  h=Y(1);
beta=Y(2);
delta=Y(3);
alpha=Y(4);
sigma=Y(5);
rho=Y(6);
omega=Y(7);
K_aminus1=(1/beta-1+delta)/alpha;
CK=K_aminus1-delta;
A=[-sigma/(1-h) 0; 0 0];
B=[sigma*(1+h)/(1-h) alpha*beta*K_aminus1*(alpha-1);-CK -1];
C=[-sigma*h/(1-h) 0;0 (1-delta+alpha*K_aminus1)];
D=[beta*alpha*K_aminus1*rho; K_aminus1];
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_dynare;
[errors_dynare] = matrix_quadratic_backward_errors(inputs)
 save('dynare_results','rp_dynare','ALPHA_ZS_dynare','BETA_ZS_dynare','log_c_std_dynare','eig_seperation_dynare','eigval_dynare','errors_dynare')
     path(p); 
 clear all
 p = path; 
 load('dynare_location')
 addpath(dynare_location)
 dynare habit_cr noclearall;
 eigenvalues=oo_.dr.eigval(abs(oo_.dr.eigval-rho)>1e-14);
 eig_seperation=eigenvalues(abs(eigenvalues)>1)-eigenvalues(abs(eigenvalues)<1)';[a b]=min(abs(eig_seperation(:)));eig_seperation_dynare_cr=eig_seperation(b)
[sortabs,sortindex] = sort(abs(eigenvalues));
eigval_dynare_cr=eigenvalues(sortindex); eigval_dynare_cr=[eigval_dynare_cr;1/0]


 Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
 h=Y(1);
beta=Y(2);
delta=Y(3);
alpha=Y(4);
sigma=Y(5);
rho=Y(6);
omega=Y(7);
K_aminus1=(1/beta-1+delta)/alpha;
CK=K_aminus1-delta;
A=[-sigma/(1-h) 0; 0 0];
B=[sigma*(1+h)/(1-h) alpha*beta*K_aminus1*(alpha-1);-CK -1];
C=[-sigma*h/(1-h) 0;0 (1-delta+alpha*K_aminus1)];
D=[beta*alpha*K_aminus1*rho; K_aminus1];
inputs.A=A;
inputs.B=B;
inputs.C=C;
inputs.X=ALPHA_ZS_dynare_cr;
[errors_dynare_cr] = matrix_quadratic_backward_errors(inputs)
 save('dynare_cr_results','rp_dynare_cr','ALPHA_ZS_dynare_cr','BETA_ZS_dynare_cr','log_c_std_dynare_cr','eig_seperation_dynare_cr','eigval_dynare_cr','errors_dynare_cr')
  
    path(p);
 clear all
 
 
results=nan(5,9);
 eig_results=nan(5,8);
 num_results=nan(4,8);
 results([1 3],1)=[7.8; 100*0.00566];
 load('exact_results')
 results([1 3],2)=[rp_exact; 100*log_c_std_exact];
 eig_results(:,1)=[eig_seperation_exact; eigval_exact];
 num_results(:,1)=[eig_seperation_exact; errors_exact.conditioning; errors_exact.lower_backward_error; errors_exact.upper_backward_error];
 load('solab_results')
 results(:,3)=[rp_solab; rp_exact-rp_solab; 100*log_c_std_solab;max(max(abs(X_exact-ALPHA_ZS_solab)));max(max(abs(Q_exact-BETA_ZS_solab)))];
 eig_results(:,2)=[eig_seperation_solab; eigval_exact-eigval_solab];
 num_results(:,2)=[eig_seperation_solab; errors_solab.conditioning; errors_solab.lower_backward_error; errors_solab.upper_backward_error];
 load('gensys_results')
 results(:,4)=[rp_gensys; rp_exact-rp_gensys; 100*log_c_std_gensys;max(max(abs(X_exact-ALPHA_ZS_gensys)));max(max(abs(Q_exact-BETA_ZS_gensys)))];
 eig_results(:,3)=[eig_seperation_gensys; eigval_exact-eigval_gensys];
 num_results(:,3)=[eig_seperation_gensys; errors_gensys.conditioning; errors_gensys.lower_backward_error; errors_gensys.upper_backward_error];
load('uhlig_results')
  results(:,5)=[rp_uhlig; rp_exact-rp_uhlig; 100*log_c_std_uhlig;max(max(abs(X_exact-ALPHA_ZS_uhlig)));max(max(abs(Q_exact-BETA_ZS_uhlig)))];
eig_results(:,4)=[eig_seperation_uhlig; eigval_exact-eigval_uhlig];
 num_results(:,4)=[eig_seperation_uhlig; errors_uhlig.conditioning; errors_uhlig.lower_backward_error; errors_uhlig.upper_backward_error];
load('dynare_results')
  results(:,6)=[rp_dynare; rp_exact-rp_dynare; 100*log_c_std_dynare;max(max(abs(X_exact-ALPHA_ZS_dynare)));max(max(abs(Q_exact-BETA_ZS_dynare)))];
eig_results(:,5)=[eig_seperation_dynare; eigval_exact-eigval_dynare];
 num_results(:,5)=[eig_seperation_dynare; errors_dynare.conditioning; errors_dynare.lower_backward_error; errors_dynare.upper_backward_error];
load('aim_results')
  results(:,7)=[rp_aim; rp_exact-rp_aim;100*log_c_std_aim;max(max(abs(X_exact-ALPHA_ZS_aim)));max(max(abs(Q_exact-BETA_ZS_aim)))];
eig_results(:,6)=[eig_seperation_aim; eigval_exact-eigval_aim];
 num_results(:,6)=[eig_seperation_aim; errors_aim.conditioning; errors_aim.lower_backward_error; errors_aim.upper_backward_error];
load('bp_results') 
results(:,8)=[rp_bp; rp_exact-rp_bp;100*log_c_std_bp;max(max(abs(X_exact-ALPHA_ZS_bp)));max(max(abs(Q_exact-BETA_ZS_bp)))];
eig_results(:,7)=[eig_seperation_bp; eigval_exact-eigval_bp];
 num_results(:,7)=[eig_seperation_bp; errors_bp.conditioning; errors_bp.lower_backward_error; errors_bp.upper_backward_error];
load('dynare_cr_results')
  results(:,9)=[rp_dynare_cr; rp_exact-rp_dynare_cr; 100*log_c_std_dynare_cr;max(max(abs(X_exact-ALPHA_ZS_dynare_cr)));max(max(abs(Q_exact-BETA_ZS_dynare_cr)))];
eig_results(:,8)=[eig_seperation_dynare_cr; eigval_exact-eigval_dynare_cr];
 num_results(:,8)=[eig_seperation_dynare_cr; errors_dynare_cr.conditioning; errors_dynare_cr.lower_backward_error; errors_dynare_cr.upper_backward_error];

 load('Y')
 Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
%final_results=[results(1:3,:);NaN eig_results(1,:); NaN NaN max(abs(eig_results(2:4,2:end)));results(4:5,:)];
final_results=[results(1:3,:);results(4:5,:);NaN eig_results(1,:); NaN NaN max(abs(eig_results(2:4,2:end))); NaN(4,1) num_results];
final_results=[results(1:5,:);NaN NaN max(abs(eig_results(2:4,2:end))); NaN(4,1) num_results];
 save('tables','results','eig_results','Y','final_results','num_results')
 disp(compose('%1.3g',final_results))
 one_vector=[1 1 0 1 0 1 0];
 disp(compose('%1.3e',Y'-one_vector))

 