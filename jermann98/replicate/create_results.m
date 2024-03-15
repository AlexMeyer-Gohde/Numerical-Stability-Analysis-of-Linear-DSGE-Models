dynare_location='C:\dynare\5.1\matlab';
 save('dynare_location','dynare_location')

addpath('..\..\algorithm\')
addpath('..\..\solution_methods\')

clear all
 p = path; 
 load('dynare_location')
 addpath(dynare_location)
 dynare jermann98_moments_comparison noclearall;
 
 
 %%%%%%%%%%%%%%%%%%
[ALPHA_ZS_aim,BETA_ZS_aim,sorted_eigenvalues,existence_uniqueness]=...
    sp_solve(matrix_quadratic.A,matrix_quadratic.B,matrix_quadratic.C,zeros(M_.endo_nbr,M_.exo_nbr),matrix_quadratic.D,zeros(M_.exo_nbr,M_.exo_nbr),M_.endo_nbr,M_.endo_nbr,M_.exo_nbr,1);

%[vx_aim, u] =lyapunov_symm(ALPHA_ZS_aim,BETA_ZS_aim*M_.Sigma_e*BETA_ZS_aim',options_.lyapunov_fixed_point_tol,1+eps,options_.lyapunov_complex_threshold,[],options_.debug);

oo_aim.dr=oo_.dr;
oo_aim.dr.ghu=BETA_ZS_aim(oo_.dr.order_var,:);
oo_aim.dr.ghx=ALPHA_ZS_aim(oo_.dr.order_var,oo_.dr.order_var);
oo_aim.dr.ghx=oo_aim.dr.ghx(:,nstatic+1:end-nfwrd);
option=options_;
option.varlist=M_.endo_names;
option.qz_criterium=1+eps;
moments_aim=compute_model_moments(oo_aim.dr,M_,option);
oo_aim.var=moments_aim{1};

std_Y_aim=(oo_aim.var(Y_gr_index,Y_gr_index))^(1/2);
std_C_aim=(oo_aim.var(C_gr_index,C_gr_index)/oo_aim.var(Y_gr_index,Y_gr_index))^(1/2);
std_I_aim=(oo_aim.var(I_gr_index,I_gr_index)/oo_aim.var(Y_gr_index,Y_gr_index))^(1/2);
rp_aim=-400*oo_aim.var(R_index,M_index);

rf_aim=-400*(oo_.steady_state(M_index)+0.5*oo_aim.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_aim.dr.ghu(oo_.dr.inv_order_var(M_index),:)');

actual_aim=[std_Y_aim std_C_aim std_I_aim rf_aim rp_aim];

matrix_quadratic.X=ALPHA_ZS_aim;
matrix_quadratic.P=ALPHA_ZS_aim;
matrix_quadratic.Q=BETA_ZS_aim;
%quad_error=matrix_quadratic_backward_errors_short(matrix_quadratic);
[errors_aim] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_aim=matrix_quadratic;
 save('aim_results','targets','actual_aim','matrix_quadratic_aim','errors_aim')
%%%%%%%%%%%%%%%%%%%%%

[BETA_ZS_solab,ll,ff,ALPHA_ZS_solab,eigval_solab]  = solab([ zeros(M_.endo_nbr,M_.endo_nbr) -matrix_quadratic.A;  eye(M_.endo_nbr,M_.endo_nbr) zeros(M_.endo_nbr,M_.endo_nbr)],[matrix_quadratic.C matrix_quadratic.B;   zeros(M_.endo_nbr,M_.endo_nbr) eye(M_.endo_nbr,M_.endo_nbr)],zeros(M_.exo_nbr,M_.exo_nbr),[matrix_quadratic.D; zeros(M_.endo_nbr,M_.exo_nbr)],M_.endo_nbr);

oo_solab.dr=oo_.dr;
oo_solab.dr.ghu=BETA_ZS_solab(oo_.dr.order_var,:);
oo_solab.dr.ghx=ALPHA_ZS_solab(oo_.dr.order_var,oo_.dr.order_var);
oo_solab.dr.ghx=oo_solab.dr.ghx(:,nstatic+1:end-nfwrd);
option=options_;
option.varlist=M_.endo_names;
option.qz_criterium=1+eps;
moments_solab=compute_model_moments(oo_solab.dr,M_,option);
oo_solab.var=moments_solab{1};

std_Y_solab=(oo_solab.var(Y_gr_index,Y_gr_index))^(1/2);
std_C_solab=(oo_solab.var(C_gr_index,C_gr_index)/oo_solab.var(Y_gr_index,Y_gr_index))^(1/2);
std_I_solab=(oo_solab.var(I_gr_index,I_gr_index)/oo_solab.var(Y_gr_index,Y_gr_index))^(1/2);
rp_solab=-400*oo_solab.var(R_index,M_index);

rf_solab=-400*(oo_.steady_state(M_index)+0.5*oo_solab.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_solab.dr.ghu(oo_.dr.inv_order_var(M_index),:)');

actual_solab=[std_Y_solab std_C_solab std_I_solab rf_solab rp_solab];

matrix_quadratic.X=ALPHA_ZS_solab;
matrix_quadratic.P=ALPHA_ZS_solab;
matrix_quadratic.Q=BETA_ZS_solab;
[errors_solab] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_solab=matrix_quadratic;
save('solab_results','targets','actual_solab','matrix_quadratic_solab','errors_solab')

%%%%%%%%%%%%%%%%%%%%%
addpath('gensys')
g0=[matrix_quadratic.A, matrix_quadratic.B;zeros(size(matrix_quadratic.A)) eye(size(matrix_quadratic.A))];
%g1=[zeros(size(matrix_quadratic.A)), -matrix_quadratic.C; eye(size(matrix_quadratic.A)) zeros(size(matrix_quadratic.A))];
g1=[zeros(size(matrix_quadratic.A)), -matrix_quadratic.C; eye(size(matrix_quadratic.A)) zeros(size(matrix_quadratic.A))];
c=zeros(2*size(matrix_quadratic.A,1),1);
psi=[-matrix_quadratic.D; zeros(size(matrix_quadratic.D))];
%pi=[matrix_quadratic.A; zeros(size(matrix_quadratic.A))];
pi=[zeros(size(matrix_quadratic.A)); eye(size(matrix_quadratic.A))];

[G1,Cgensys,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi);
ALPHA_ZS_gensys=G1(size(matrix_quadratic.A,1)+1:end,size(matrix_quadratic.A,1)+1:end);
BETA_ZS_gensys=impact(size(matrix_quadratic.A,1)+1:end,:);

oo_gensys.dr=oo_.dr;
oo_gensys.dr.ghu=BETA_ZS_gensys(oo_.dr.order_var,:);
oo_gensys.dr.ghx=ALPHA_ZS_gensys(oo_.dr.order_var,oo_.dr.order_var);
oo_gensys.dr.ghx=oo_gensys.dr.ghx(:,nstatic+1:end-nfwrd);
option=options_;
option.varlist=M_.endo_names;
option.qz_criterium=1+eps;
moments_gensys=compute_model_moments(oo_gensys.dr,M_,option);
oo_gensys.var=moments_gensys{1};

std_Y_gensys=(oo_gensys.var(Y_gr_index,Y_gr_index))^(1/2);
std_C_gensys=(oo_gensys.var(C_gr_index,C_gr_index)/oo_gensys.var(Y_gr_index,Y_gr_index))^(1/2);
std_I_gensys=(oo_gensys.var(I_gr_index,I_gr_index)/oo_gensys.var(Y_gr_index,Y_gr_index))^(1/2);
rp_gensys=-400*oo_gensys.var(R_index,M_index);

rf_gensys=-400*(oo_.steady_state(M_index)+0.5*oo_gensys.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_gensys.dr.ghu(oo_.dr.inv_order_var(M_index),:)');

actual_gensys=[std_Y_gensys std_C_gensys std_I_gensys rf_gensys rp_gensys];

matrix_quadratic.P=ALPHA_ZS_gensys;
matrix_quadratic.Q=BETA_ZS_gensys;
matrix_quadratic.X=ALPHA_ZS_gensys;
[errors_gensys] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_gensys=matrix_quadratic;

rmpath('gensys')
save('gensys_results','targets','actual_gensys','matrix_quadratic_gensys','errors_gensys')

%%%%%%%%%%%%%%%%%%%%%
addpath('uhlig')

AA=[]; BB=[]; CC=[]; DD=[];
FF=matrix_quadratic.A;
GG=matrix_quadratic.B;
HH=matrix_quadratic.C;
JJ=[]; KK=[];
MM=matrix_quadratic.D;
LL=zeros(size(MM));
NN=zeros(size(matrix_quadratic.D,2),size(matrix_quadratic.D,2));

[l_equ,m_states] = size(FF);
[l_equ,n_endog ] = size(FF);
[l_equ,k_exog  ] = size(MM);
warnings=[];
options;
solve;
ALPHA_ZS_uhlig=PP;
BETA_ZS_uhlig=QQ;

oo_uhlig.dr=oo_.dr;
oo_uhlig.dr.ghu=BETA_ZS_uhlig(oo_.dr.order_var,:);
oo_uhlig.dr.ghx=ALPHA_ZS_uhlig(oo_.dr.order_var,oo_.dr.order_var);
oo_uhlig.dr.ghx=oo_uhlig.dr.ghx(:,nstatic+1:end-nfwrd);
option=options_;
option.varlist=M_.endo_names;
option.qz_criterium=1+eps;
moments_uhlig=compute_model_moments(oo_uhlig.dr,M_,option);
oo_uhlig.var=moments_uhlig{1};

std_Y_uhlig=(oo_uhlig.var(Y_gr_index,Y_gr_index))^(1/2);
std_C_uhlig=(oo_uhlig.var(C_gr_index,C_gr_index)/oo_uhlig.var(Y_gr_index,Y_gr_index))^(1/2);
std_I_uhlig=(oo_uhlig.var(I_gr_index,I_gr_index)/oo_uhlig.var(Y_gr_index,Y_gr_index))^(1/2);
rp_uhlig=-400*oo_uhlig.var(R_index,M_index);

rf_uhlig=-400*(oo_.steady_state(M_index)+0.5*oo_uhlig.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_uhlig.dr.ghu(oo_.dr.inv_order_var(M_index),:)');

actual_uhlig=[std_Y_uhlig std_C_uhlig std_I_uhlig rf_uhlig rp_uhlig];

eig_seperation=Xi_sortval(abs(Xi_sortval)>1)-Xi_sortval(abs(Xi_sortval)<1)';[a b]=min(abs(eig_seperation(:)));eig_seperation_uhlig=eig_seperation(b);
eigval_uhlig=Xi_sortval ;
matrix_quadratic.P=ALPHA_ZS_uhlig;
matrix_quadratic.Q=BETA_ZS_uhlig;
matrix_quadratic.X=ALPHA_ZS_uhlig;
[errors_uhlig] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_uhlig=matrix_quadratic;
rmpath('uhlig')
save('uhlig_results','targets','actual_uhlig','matrix_quadratic_uhlig','errors_uhlig')
%%%%%%%%%%%%%%%%%%%%%
[ALPHA_ZS_bp,BETA_ZS_bp]=bp(matrix_quadratic.A,-matrix_quadratic.B,matrix_quadratic.C,matrix_quadratic.D,zeros(M_.endo_nbr,M_.exo_nbr),zeros(M_.exo_nbr,M_.exo_nbr),80e5,1e-10);

oo_bp.dr=oo_.dr;
oo_bp.dr.ghu=BETA_ZS_bp(oo_.dr.order_var,:);
oo_bp.dr.ghx=ALPHA_ZS_bp(oo_.dr.order_var,oo_.dr.order_var);
oo_bp.dr.ghx=oo_bp.dr.ghx(:,nstatic+1:end-nfwrd);
option=options_;
option.varlist=M_.endo_names;
option.qz_criterium=1+eps;
moments_bp=compute_model_moments(oo_bp.dr,M_,option);
oo_bp.var=moments_bp{1};

std_Y_bp=(oo_bp.var(Y_gr_index,Y_gr_index))^(1/2);
std_C_bp=(oo_bp.var(C_gr_index,C_gr_index)/oo_bp.var(Y_gr_index,Y_gr_index))^(1/2);
std_I_bp=(oo_bp.var(I_gr_index,I_gr_index)/oo_bp.var(Y_gr_index,Y_gr_index))^(1/2);
rp_bp=-400*oo_bp.var(R_index,M_index);

rf_bp=-400*(oo_.steady_state(M_index)+0.5*oo_bp.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_bp.dr.ghu(oo_.dr.inv_order_var(M_index),:)');

actual_bp=[std_Y_bp std_C_bp std_I_bp rf_bp rp_bp];

matrix_quadratic.P=ALPHA_ZS_bp;
matrix_quadratic.Q=BETA_ZS_bp;
matrix_quadratic.X=ALPHA_ZS_bp;
[errors_bp] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_bp=matrix_quadratic;
save('bp_results','targets','actual_bp','matrix_quadratic_bp','errors_bp')

% results=nan(6,8);
% num_results=nan(4,7);
% load('aim_results')
% results([1 3:6],1)=targets([5 4 1 2 3])';
% load('solab_results')
% results([1 3:6],2)=actual_solab([5 4 1 2 3])';
% results(2,2)=actual_aim(5)-actual_solab(5);
% num_results(:,1)=[eig_seperation_solab; errors_solab(3); errors_solab(1:2)];
% load('gensys_results')
% results([1 3:6],3)=actual_gensys([5 4 1 2 3])';
% results(2,3)=actual_aim(5)-actual_gensys(5);
% num_results(:,2)=[eig_seperation_gensys; errors_gensys(3); errors_gensys(1:2)];
% load('uhlig_results')
% results([1 3:6],4)=actual_uhlig([5 4 1 2 3])';
% results(2,4)=actual_aim(5)-actual_uhlig(5);
% num_results(:,3)=[eig_seperation_uhlig; errors_uhlig(3); errors_uhlig(1:2)];
% load('dynare_results')
% results([1 3:6],5)=actual([5 4 1 2 3])';
% results(2,5)=actual_aim(5)-actual(5);
% num_results(:,4)=[eig_seperation; errors(3); errors(1:2)];
% load('aim_results')
% results([1 3:6],6)=actual_aim([5 4 1 2 3])';
% results(2,6)=NaN;
% num_results(:,5)=[eig_seperation_aim; errors_aim(3); errors_aim(1:2)];
% load('dynare_cr_results')
% results([1 3:6],7)=actual_dynare_cr([5 4 1 2 3])';
% results(2,7)=actual_aim(5)-actual_dynare_cr(5);
% num_results(:,6)=[eig_seperation_dynare_cr; errors_dynare_cr(3); errors_dynare_cr(1:2)];
% 
% final_results=[results;[nan(4,1) num_results]];
%  disp(compose('%1.3g',final_results))