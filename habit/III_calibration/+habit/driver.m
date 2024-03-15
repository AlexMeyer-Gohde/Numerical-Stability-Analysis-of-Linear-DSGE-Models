%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'habit';
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
M_.exo_names(1) = {'e'};
M_.exo_names_tex(1) = {'e'};
M_.exo_names_long(1) = {'e'};
M_.endo_names = cell(3,1);
M_.endo_names_tex = cell(3,1);
M_.endo_names_long = cell(3,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'k'};
M_.endo_names_tex(2) = {'k'};
M_.endo_names_long(2) = {'k'};
M_.endo_names(3) = {'z'};
M_.endo_names_tex(3) = {'z'};
M_.endo_names_long(3) = {'z'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'h'};
M_.param_names_tex(1) = {'h'};
M_.param_names_long(1) = {'h'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'sigma'};
M_.param_names_tex(5) = {'sigma'};
M_.param_names_long(5) = {'sigma'};
M_.param_names(6) = {'rho'};
M_.param_names_tex(6) = {'rho'};
M_.param_names_long(6) = {'rho'};
M_.param_names(7) = {'K_aminus1'};
M_.param_names_tex(7) = {'K\_aminus1'};
M_.param_names_long(7) = {'K_aminus1'};
M_.param_names(8) = {'CK'};
M_.param_names_tex(8) = {'CK'};
M_.param_names_long(8) = {'CK'};
M_.param_names(9) = {'omega'};
M_.param_names_tex(9) = {'omega'};
M_.param_names_long(9) = {'omega'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 3;
M_.param_nbr = 9;
M_.orig_endo_nbr = 3;
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
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 3;
M_.eq_nbr = 3;
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
 1 4 7;
 2 5 0;
 3 6 0;]';
M_.nstatic = 0;
M_.nfwrd   = 0;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 1;
M_.nspred   = 3;
M_.ndynamic   = 3;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'z' ;
};
M_.mapping.c.eqidx = [1 2 ];
M_.mapping.k.eqidx = [1 2 ];
M_.mapping.z.eqidx = [1 2 3 ];
M_.mapping.e.eqidx = [3 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 2 3 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(3, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(3, 1), 'log_deflator', cell(3, 1), 'growth_factor', cell(3, 1), 'log_growth_factor', cell(3, 1));
M_.NNZDerivatives = [12; 0; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load('Y')
Y=Y_1;
M_.params(1) = Y(1);
h = M_.params(1);
M_.params(2) = Y(2);
beta = M_.params(2);
M_.params(3) = Y(3);
delta = M_.params(3);
M_.params(4) = Y(4);
alpha = M_.params(4);
M_.params(5) = Y(5);
sigma = M_.params(5);
M_.params(6) = Y(6);
rho = M_.params(6);
M_.params(9) = Y(7);
omega = M_.params(9);
M_.params(7) = (1/M_.params(2)-1+M_.params(3))/M_.params(4);
K_aminus1 = M_.params(7);
M_.params(8) = M_.params(7)-M_.params(3);
CK = M_.params(8);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(9)^2;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 0;
options_.noprint = true;
options_.order = 1;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
BETA_ZS_dynare=oo_.dr.ghu(oo_.dr.inv_order_var,:);BETA_ZS_dynare=BETA_ZS_dynare(1:2,:);
ALPHA_ZS_dynare=oo_.dr.ghx(oo_.dr.inv_order_var,oo_.dr.inv_order_var);ALPHA_ZS_dynare=ALPHA_ZS_dynare(1:2,1:2);
rp_dynare=4*100*sigma/(1-h)*BETA_ZS_dynare(1)*alpha*K_aminus1*beta*(omega)^2
AA=[ALPHA_ZS_dynare(1,1)-1  ALPHA_ZS_dynare(1,2) BETA_ZS_dynare(1)];
BB=[ALPHA_ZS_dynare BETA_ZS_dynare;0 0 rho ];
DD=[0;0;1];
log_c_std_dynare=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'habit_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
