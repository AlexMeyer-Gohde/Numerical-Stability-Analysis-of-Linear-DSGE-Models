%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'jermann98_rep_def_R_rr_compare';
M_.dynare_version = '5.1';
oo_.dynare_version = '5.1';
options_.dynare_version = '5.1';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'e_z'};
M_.exo_names_tex(1) = {'e\_z'};
M_.exo_names_long(1) = {'e_z'};
M_.exo_names(2) = {'e_a'};
M_.exo_names_tex(2) = {'e\_a'};
M_.exo_names_long(2) = {'e_a'};
M_.endo_names = cell(17,1);
M_.endo_names_tex = cell(17,1);
M_.endo_names_long = cell(17,1);
M_.endo_names(1) = {'k'};
M_.endo_names_tex(1) = {'k'};
M_.endo_names_long(1) = {'k'};
M_.endo_names(2) = {'c'};
M_.endo_names_tex(2) = {'c'};
M_.endo_names_long(2) = {'c'};
M_.endo_names(3) = {'i'};
M_.endo_names_tex(3) = {'i'};
M_.endo_names_long(3) = {'i'};
M_.endo_names(4) = {'y'};
M_.endo_names_tex(4) = {'y'};
M_.endo_names_long(4) = {'y'};
M_.endo_names(5) = {'N'};
M_.endo_names_tex(5) = {'N'};
M_.endo_names_long(5) = {'N'};
M_.endo_names(6) = {'C_gr'};
M_.endo_names_tex(6) = {'C\_gr'};
M_.endo_names_long(6) = {'C_gr'};
M_.endo_names(7) = {'Y_gr'};
M_.endo_names_tex(7) = {'Y\_gr'};
M_.endo_names_long(7) = {'Y_gr'};
M_.endo_names(8) = {'I_gr'};
M_.endo_names_tex(8) = {'I\_gr'};
M_.endo_names_long(8) = {'I_gr'};
M_.endo_names(9) = {'mu'};
M_.endo_names_tex(9) = {'mu'};
M_.endo_names_long(9) = {'mu'};
M_.endo_names(10) = {'cost'};
M_.endo_names_tex(10) = {'cost'};
M_.endo_names_long(10) = {'cost'};
M_.endo_names(11) = {'cost_der'};
M_.endo_names_tex(11) = {'cost\_der'};
M_.endo_names_long(11) = {'cost_der'};
M_.endo_names(12) = {'epsilon_a'};
M_.endo_names_tex(12) = {'epsilon\_a'};
M_.endo_names_long(12) = {'epsilon_a'};
M_.endo_names(13) = {'z'};
M_.endo_names_tex(13) = {'z'};
M_.endo_names_long(13) = {'z'};
M_.endo_names(14) = {'w'};
M_.endo_names_tex(14) = {'w'};
M_.endo_names_long(14) = {'w'};
M_.endo_names(15) = {'d'};
M_.endo_names_tex(15) = {'d'};
M_.endo_names_long(15) = {'d'};
M_.endo_names(16) = {'R'};
M_.endo_names_tex(16) = {'R'};
M_.endo_names_long(16) = {'R'};
M_.endo_names(17) = {'M'};
M_.endo_names_tex(17) = {'M'};
M_.endo_names_long(17) = {'M'};
M_.endo_partitions = struct();
M_.param_names = cell(12,1);
M_.param_names_tex = cell(12,1);
M_.param_names_long = cell(12,1);
M_.param_names(1) = {'tau'};
M_.param_names_tex(1) = {'tau'};
M_.param_names_long(1) = {'tau'};
M_.param_names(2) = {'b'};
M_.param_names_tex(2) = {'b'};
M_.param_names_long(2) = {'b'};
M_.param_names(3) = {'beta_star'};
M_.param_names_tex(3) = {'beta\_star'};
M_.param_names_long(3) = {'beta_star'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'a_bar'};
M_.param_names_tex(5) = {'a\_bar'};
M_.param_names_long(5) = {'a_bar'};
M_.param_names(6) = {'rho_z'};
M_.param_names_tex(6) = {'rho\_z'};
M_.param_names_long(6) = {'rho_z'};
M_.param_names(7) = {'sigma_a_bar'};
M_.param_names_tex(7) = {'sigma\_a\_bar'};
M_.param_names_long(7) = {'sigma_a_bar'};
M_.param_names(8) = {'sigma_z_bar'};
M_.param_names_tex(8) = {'sigma\_z\_bar'};
M_.param_names_long(8) = {'sigma_z_bar'};
M_.param_names(9) = {'cost_xi'};
M_.param_names_tex(9) = {'cost\_xi'};
M_.param_names_long(9) = {'cost_xi'};
M_.param_names(10) = {'delta'};
M_.param_names_tex(10) = {'delta'};
M_.param_names_long(10) = {'delta'};
M_.param_names(11) = {'N_ss'};
M_.param_names_tex(11) = {'N\_ss'};
M_.param_names_long(11) = {'N_ss'};
M_.param_names(12) = {'Scac'};
M_.param_names_tex(12) = {'Scac'};
M_.param_names_long(12) = {'Scac'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 17;
M_.param_nbr = 12;
M_.orig_endo_nbr = 17;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 17;
M_.eq_nbr = 17;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 8 0;
 2 9 25;
 3 10 0;
 4 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 5 16 26;
 0 17 0;
 6 18 0;
 0 19 27;
 7 20 0;
 0 21 0;
 0 22 0;
 0 23 28;
 0 24 0;]';
M_.nstatic = 8;
M_.nfwrd   = 2;
M_.npred   = 5;
M_.nboth   = 2;
M_.nsfwrd   = 4;
M_.nspred   = 7;
M_.ndynamic   = 9;
M_.dynamic_tmp_nbr = [9; 5; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'mu' ;
  2 , 'name' , '2' ;
  3 , 'name' , '3' ;
  4 , 'name' , '4' ;
  5 , 'name' , '5' ;
  6 , 'name' , '6' ;
  7 , 'name' , 'cost' ;
  8 , 'name' , 'cost_der' ;
  9 , 'name' , '9' ;
  10 , 'name' , '10' ;
  11 , 'name' , '11' ;
  12 , 'name' , '12' ;
  13 , 'name' , 'C_gr' ;
  14 , 'name' , 'Y_gr' ;
  15 , 'name' , 'I_gr' ;
  16 , 'name' , 'epsilon_a' ;
  17 , 'name' , 'z' ;
};
M_.mapping.k.eqidx = [3 4 7 8 11 ];
M_.mapping.c.eqidx = [1 6 13 ];
M_.mapping.i.eqidx = [3 6 7 8 10 11 15 ];
M_.mapping.y.eqidx = [4 6 9 10 11 14 ];
M_.mapping.N.eqidx = [4 5 9 10 ];
M_.mapping.C_gr.eqidx = [13 ];
M_.mapping.Y_gr.eqidx = [14 ];
M_.mapping.I_gr.eqidx = [15 ];
M_.mapping.mu.eqidx = [1 2 12 ];
M_.mapping.cost.eqidx = [3 7 11 ];
M_.mapping.cost_der.eqidx = [8 11 ];
M_.mapping.epsilon_a.eqidx = [1 2 3 4 7 8 11 12 13 14 15 16 ];
M_.mapping.z.eqidx = [4 17 ];
M_.mapping.w.eqidx = [9 10 ];
M_.mapping.d.eqidx = [10 ];
M_.mapping.R.eqidx = [11 12 ];
M_.mapping.M.eqidx = [2 ];
M_.mapping.e_z.eqidx = [17 ];
M_.mapping.e_a.eqidx = [16 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 2 3 4 9 11 13 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(17, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(12, 1);
M_.endo_trends = struct('deflator', cell(17, 1), 'log_deflator', cell(17, 1), 'growth_factor', cell(17, 1), 'log_growth_factor', cell(17, 1));
M_.NNZDerivatives = [69; -1; -1; ];
M_.static_tmp_nbr = [8; 5; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load('Y_asset_test')
M_.params(1) = X_0(1);
tau = M_.params(1);
M_.params(2) = X_0(2);
b = M_.params(2);
M_.params(4) = 0.36;
alpha = M_.params(4);
M_.params(5) = 0.005;
a_bar = M_.params(5);
M_.params(3) = X_0(3);
beta_star = M_.params(3);
M_.params(10) = 0.025;
delta = M_.params(10);
M_.params(9) = X_0(4);
cost_xi = M_.params(9);
M_.params(6) = 0.99;
rho_z = M_.params(6);
M_.params(7) = 0;
sigma_a_bar = M_.params(7);
M_.params(8) = X_0(5);
sigma_z_bar = M_.params(8);
M_.params(11) = 1;
N_ss = M_.params(11);
M_.params(12) = 1;
Scac = M_.params(12);
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
options_.drop = 0;
options_.irf = 0;
options_.noprint = true;
options_.order = 1;
options_.periods = 0;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
options_.drop = 0;
options_.irf = 0;
options_.noprint = true;
options_.order = 1;
options_.periods = 0;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
C_gr_index=find(strcmp(M_.endo_names,'C_gr'));
Y_gr_index=find(strcmp(M_.endo_names,'Y_gr'));
I_gr_index=find(strcmp(M_.endo_names,'I_gr'));
M_index=find(strcmp(M_.endo_names,'M'));
R_index=find(strcmp(M_.endo_names,'R'));
Rf_index=find(strcmp(M_.endo_names,'Rf'));
std_Y=(oo_.var(Y_gr_index,Y_gr_index))^(1/2);
std_C=(oo_.var(C_gr_index,C_gr_index)/oo_.var(Y_gr_index,Y_gr_index))^(1/2);
std_I=(oo_.var(I_gr_index,I_gr_index)/oo_.var(Y_gr_index,Y_gr_index))^(1/2);
rp=-400*oo_.var(R_index,M_index);
rf=-400*(oo_.steady_state(M_index)+0.5*oo_.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_.dr.ghu(oo_.dr.inv_order_var(M_index),:)');
100*(4*oo_.var(R_index,R_index))^(1/2);
100*(4*oo_.var(Rf_index,Rf_index))^(1/2);
targets=[0.01 0.51 2.65 0.8 6.18];
actual=[std_Y std_C std_I rf rp];
create_matrix_quadratic_from_dynare
ALPHA_ZS_dynare=[zeros(n,nstatic) oo_.dr.ghx zeros(n,nfwrd)];
ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
matrix_quadratic.X=ALPHA_ZS_dynare;
BETA_ZS_dynare=oo_.dr.ghu; BETA_ZS_dynare=BETA_ZS_dynare(oo_.dr.inv_order_var,:);
matrix_quadratic.W=BETA_ZS_dynare;
[errors_dynare] = matrix_quadratic_backward_errors_vector(matrix_quadratic)
[impact_errors_dynare] = matrix_impact_backward_errors_vector(matrix_quadratic)
eigval=oo_.dr.eigval;
l_1=eigval(abs(eigval)>1);
l_2=eigval(abs(eigval)<1);
chordal_dist=abs(l_1-l_2')./((abs(l_1).^2+1).^(1/2).*(abs(l_2').^2+1).^(1/2));
eig_seperation=min(min(chordal_dist));
save('dynare_results','targets','actual','matrix_quadratic','errors_dynare','eig_seperation','eigval','impact_errors_dynare')
options_.dr_cycle_reduction = true;
options_.dr_cycle_reduction_tol = 1e-16;
options_.irf = 0;
options_.noprint = true;
options_.order = 1;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
std_Y=(oo_.var(Y_gr_index,Y_gr_index))^(1/2);
std_C=(oo_.var(C_gr_index,C_gr_index)/oo_.var(Y_gr_index,Y_gr_index))^(1/2);
std_I=(oo_.var(I_gr_index,I_gr_index)/oo_.var(Y_gr_index,Y_gr_index))^(1/2);
rp=-400*oo_.var(R_index,M_index);
rf=-400*(oo_.steady_state(M_index)+0.5*oo_.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_.dr.ghu(oo_.dr.inv_order_var(M_index),:)');
100*(4*oo_.var(R_index,R_index))^(1/2);
100*(4*oo_.var(Rf_index,Rf_index))^(1/2);
actual_dynare_cr=[std_Y std_C std_I rf rp];
ALPHA_ZS_dynare=[zeros(n,nstatic) oo_.dr.ghx zeros(n,nfwrd)];
ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
matrix_quadratic.X=ALPHA_ZS_dynare;
BETA_ZS_dynare=oo_.dr.ghu; BETA_ZS_dynare=BETA_ZS_dynare(oo_.dr.inv_order_var,:);
matrix_quadratic.W=BETA_ZS_dynare;
[errors_dynare_cr] = matrix_quadratic_backward_errors_vector(matrix_quadratic)
[impact_errors_dynare_cr] = matrix_impact_backward_errors_vector(matrix_quadratic)
matrix_quadratic_dynare_cr=matrix_quadratic;
eigval=oo_.dr.eigval;
eigval_dynare_cr=oo_.dr.eigval;
l_1=eigval(abs(eigval)>1);
l_2=eigval(abs(eigval)<1);
chordal_dist=abs(l_1-l_2')./((abs(l_1).^2+1).^(1/2).*(abs(l_2').^2+1).^(1/2));
eig_seperation_dynare_cr=min(min(chordal_dist));
save('dynare_cr_results','targets','actual_dynare_cr','matrix_quadratic_dynare_cr','errors_dynare_cr','eig_seperation_dynare_cr','eigval_dynare_cr','impact_errors_dynare_cr')


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_rep_def_R_rr_compare_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
