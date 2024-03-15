%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'jermann98_moments_comparison';
M_.dynare_version = '5.1';
oo_.dynare_version = '5.1';
options_.dynare_version = '5.1';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'e_z'};
M_.exo_names_tex(1) = {'e\_z'};
M_.exo_names_long(1) = {'e_z'};
M_.endo_names = cell(13,1);
M_.endo_names_tex = cell(13,1);
M_.endo_names_long = cell(13,1);
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
M_.endo_names(5) = {'C_gr'};
M_.endo_names_tex(5) = {'C\_gr'};
M_.endo_names_long(5) = {'C_gr'};
M_.endo_names(6) = {'Y_gr'};
M_.endo_names_tex(6) = {'Y\_gr'};
M_.endo_names_long(6) = {'Y_gr'};
M_.endo_names(7) = {'I_gr'};
M_.endo_names_tex(7) = {'I\_gr'};
M_.endo_names_long(7) = {'I_gr'};
M_.endo_names(8) = {'mu'};
M_.endo_names_tex(8) = {'mu'};
M_.endo_names_long(8) = {'mu'};
M_.endo_names(9) = {'cost'};
M_.endo_names_tex(9) = {'cost'};
M_.endo_names_long(9) = {'cost'};
M_.endo_names(10) = {'cost_der'};
M_.endo_names_tex(10) = {'cost\_der'};
M_.endo_names_long(10) = {'cost_der'};
M_.endo_names(11) = {'z'};
M_.endo_names_tex(11) = {'z'};
M_.endo_names_long(11) = {'z'};
M_.endo_names(12) = {'R'};
M_.endo_names_tex(12) = {'R'};
M_.endo_names_long(12) = {'R'};
M_.endo_names(13) = {'M'};
M_.endo_names_tex(13) = {'M'};
M_.endo_names_long(13) = {'M'};
M_.endo_partitions = struct();
M_.param_names = cell(16,1);
M_.param_names_tex = cell(16,1);
M_.param_names_long = cell(16,1);
M_.param_names(1) = {'tau'};
M_.param_names_tex(1) = {'tau'};
M_.param_names_long(1) = {'tau'};
M_.param_names(2) = {'b'};
M_.param_names_tex(2) = {'b'};
M_.param_names_long(2) = {'b'};
M_.param_names(3) = {'beta'};
M_.param_names_tex(3) = {'beta'};
M_.param_names_long(3) = {'beta'};
M_.param_names(4) = {'beta_star'};
M_.param_names_tex(4) = {'beta\_star'};
M_.param_names_long(4) = {'beta_star'};
M_.param_names(5) = {'alpha'};
M_.param_names_tex(5) = {'alpha'};
M_.param_names_long(5) = {'alpha'};
M_.param_names(6) = {'a_bar'};
M_.param_names_tex(6) = {'a\_bar'};
M_.param_names_long(6) = {'a_bar'};
M_.param_names(7) = {'rho_z'};
M_.param_names_tex(7) = {'rho\_z'};
M_.param_names_long(7) = {'rho_z'};
M_.param_names(8) = {'sigma_a_bar'};
M_.param_names_tex(8) = {'sigma\_a\_bar'};
M_.param_names_long(8) = {'sigma_a_bar'};
M_.param_names(9) = {'sigma_z_bar'};
M_.param_names_tex(9) = {'sigma\_z\_bar'};
M_.param_names_long(9) = {'sigma_z_bar'};
M_.param_names(10) = {'cost_xi'};
M_.param_names_tex(10) = {'cost\_xi'};
M_.param_names_long(10) = {'cost_xi'};
M_.param_names(11) = {'cost_b'};
M_.param_names_tex(11) = {'cost\_b'};
M_.param_names_long(11) = {'cost_b'};
M_.param_names(12) = {'cost_c'};
M_.param_names_tex(12) = {'cost\_c'};
M_.param_names_long(12) = {'cost_c'};
M_.param_names(13) = {'delta'};
M_.param_names_tex(13) = {'delta'};
M_.param_names_long(13) = {'delta'};
M_.param_names(14) = {'N_ss'};
M_.param_names_tex(14) = {'N\_ss'};
M_.param_names_long(14) = {'N_ss'};
M_.param_names(15) = {'expi_expk'};
M_.param_names_tex(15) = {'expi\_expk'};
M_.param_names_long(15) = {'expi_expk'};
M_.param_names(16) = {'Scac'};
M_.param_names_tex(16) = {'Scac'};
M_.param_names_long(16) = {'Scac'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 13;
M_.param_nbr = 16;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
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
M_.orig_eq_nbr = 13;
M_.eq_nbr = 13;
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
 2 9 21;
 3 10 0;
 4 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 5 15 22;
 0 16 0;
 6 17 0;
 7 18 0;
 0 19 23;
 0 20 0;]';
M_.nstatic = 5;
M_.nfwrd   = 1;
M_.npred   = 5;
M_.nboth   = 2;
M_.nsfwrd   = 3;
M_.nspred   = 7;
M_.ndynamic   = 8;
M_.dynamic_tmp_nbr = [4; 5; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'mu' ;
  2 , 'name' , '2' ;
  3 , 'name' , '3' ;
  4 , 'name' , '4' ;
  5 , 'name' , '5' ;
  6 , 'name' , 'cost' ;
  7 , 'name' , 'cost_der' ;
  8 , 'name' , '8' ;
  9 , 'name' , '9' ;
  10 , 'name' , 'C_gr' ;
  11 , 'name' , 'Y_gr' ;
  12 , 'name' , 'I_gr' ;
  13 , 'name' , 'z' ;
};
M_.mapping.k.eqidx = [3 4 6 7 8 ];
M_.mapping.c.eqidx = [1 5 10 ];
M_.mapping.i.eqidx = [3 5 6 7 8 12 ];
M_.mapping.y.eqidx = [4 5 8 11 ];
M_.mapping.C_gr.eqidx = [10 ];
M_.mapping.Y_gr.eqidx = [11 ];
M_.mapping.I_gr.eqidx = [12 ];
M_.mapping.mu.eqidx = [1 2 9 ];
M_.mapping.cost.eqidx = [3 6 8 ];
M_.mapping.cost_der.eqidx = [7 8 ];
M_.mapping.z.eqidx = [4 13 ];
M_.mapping.R.eqidx = [8 9 ];
M_.mapping.M.eqidx = [2 ];
M_.mapping.e_z.eqidx = [13 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 2 3 4 8 10 11 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(16, 1);
M_.endo_trends = struct('deflator', cell(13, 1), 'log_deflator', cell(13, 1), 'growth_factor', cell(13, 1), 'log_growth_factor', cell(13, 1));
M_.NNZDerivatives = [45; -1; -1; ];
M_.static_tmp_nbr = [4; 3; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load('Y_calibrate')
M_.params(1) = X_0(1);
tau = M_.params(1);
M_.params(2) = X_0(2);
b = M_.params(2);
M_.params(5) = 0.36;
alpha = M_.params(5);
M_.params(6) = 0.005;
a_bar = M_.params(6);
M_.params(4) = X_0(3);
beta_star = M_.params(4);
M_.params(13) = 0.025;
delta = M_.params(13);
M_.params(10) = X_0(4);
cost_xi = M_.params(10);
M_.params(7) = 0.99;
rho_z = M_.params(7);
M_.params(9) = X_0(5);
sigma_z_bar = M_.params(9);
M_.params(14) = 1;
N_ss = M_.params(14);
M_.params(16) = 1;
Scac = M_.params(16);
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
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
matrix_quadratic.P=ALPHA_ZS_dynare;
matrix_quadratic.Q=BETA_ZS_dynare;
[errors] = dsge_backward_errors_condition_full(matrix_quadratic);
save('dynare_results','targets','actual','matrix_quadratic','errors')
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
matrix_quadratic.P=ALPHA_ZS_dynare;
matrix_quadratic.Q=BETA_ZS_dynare;
[errors_dynare_cr] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_dynare_cr=matrix_quadratic;
save('dynare_cr_results','targets','actual_dynare_cr','matrix_quadratic_dynare_cr','errors_dynare_cr')


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'jermann98_moments_comparison_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
