%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Smets_Wouters_2007_rep';
M_.dynare_version = '5.1';
oo_.dynare_version = '5.1';
options_.dynare_version = '5.1';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(7,1);
M_.exo_names_tex = cell(7,1);
M_.exo_names_long = cell(7,1);
M_.exo_names(1) = {'ea'};
M_.exo_names_tex(1) = {'ea'};
M_.exo_names_long(1) = {'ea'};
M_.exo_names(2) = {'eb'};
M_.exo_names_tex(2) = {'eb'};
M_.exo_names_long(2) = {'eb'};
M_.exo_names(3) = {'eqs'};
M_.exo_names_tex(3) = {'eqs'};
M_.exo_names_long(3) = {'eqs'};
M_.exo_names(4) = {'eg'};
M_.exo_names_tex(4) = {'eg'};
M_.exo_names_long(4) = {'eg'};
M_.exo_names(5) = {'em'};
M_.exo_names_tex(5) = {'em'};
M_.exo_names_long(5) = {'em'};
M_.exo_names(6) = {'epinf'};
M_.exo_names_tex(6) = {'epinf'};
M_.exo_names_long(6) = {'epinf'};
M_.exo_names(7) = {'ew'};
M_.exo_names_tex(7) = {'ew'};
M_.exo_names_long(7) = {'ew'};
M_.endo_names = cell(43,1);
M_.endo_names_tex = cell(43,1);
M_.endo_names_long = cell(43,1);
M_.endo_names(1) = {'labobs'};
M_.endo_names_tex(1) = {'labobs'};
M_.endo_names_long(1) = {'labobs'};
M_.endo_names(2) = {'robs'};
M_.endo_names_tex(2) = {'robs'};
M_.endo_names_long(2) = {'robs'};
M_.endo_names(3) = {'pinfobs'};
M_.endo_names_tex(3) = {'pinfobs'};
M_.endo_names_long(3) = {'pinfobs'};
M_.endo_names(4) = {'dy'};
M_.endo_names_tex(4) = {'dy'};
M_.endo_names_long(4) = {'dy'};
M_.endo_names(5) = {'dc'};
M_.endo_names_tex(5) = {'dc'};
M_.endo_names_long(5) = {'dc'};
M_.endo_names(6) = {'dinve'};
M_.endo_names_tex(6) = {'dinve'};
M_.endo_names_long(6) = {'dinve'};
M_.endo_names(7) = {'dw'};
M_.endo_names_tex(7) = {'dw'};
M_.endo_names_long(7) = {'dw'};
M_.endo_names(8) = {'ewma'};
M_.endo_names_tex(8) = {'ewma'};
M_.endo_names_long(8) = {'ewma'};
M_.endo_names(9) = {'epinfma'};
M_.endo_names_tex(9) = {'epinfma'};
M_.endo_names_long(9) = {'epinfma'};
M_.endo_names(10) = {'zcapf'};
M_.endo_names_tex(10) = {'zcapf'};
M_.endo_names_long(10) = {'zcapf'};
M_.endo_names(11) = {'rkf'};
M_.endo_names_tex(11) = {'rkf'};
M_.endo_names_long(11) = {'rkf'};
M_.endo_names(12) = {'kf'};
M_.endo_names_tex(12) = {'kf'};
M_.endo_names_long(12) = {'kf'};
M_.endo_names(13) = {'pkf'};
M_.endo_names_tex(13) = {'pkf'};
M_.endo_names_long(13) = {'pkf'};
M_.endo_names(14) = {'cf'};
M_.endo_names_tex(14) = {'cf'};
M_.endo_names_long(14) = {'cf'};
M_.endo_names(15) = {'invef'};
M_.endo_names_tex(15) = {'invef'};
M_.endo_names_long(15) = {'invef'};
M_.endo_names(16) = {'yf'};
M_.endo_names_tex(16) = {'yf'};
M_.endo_names_long(16) = {'yf'};
M_.endo_names(17) = {'labf'};
M_.endo_names_tex(17) = {'labf'};
M_.endo_names_long(17) = {'labf'};
M_.endo_names(18) = {'wf'};
M_.endo_names_tex(18) = {'wf'};
M_.endo_names_long(18) = {'wf'};
M_.endo_names(19) = {'rrf'};
M_.endo_names_tex(19) = {'rrf'};
M_.endo_names_long(19) = {'rrf'};
M_.endo_names(20) = {'mc'};
M_.endo_names_tex(20) = {'mc'};
M_.endo_names_long(20) = {'mc'};
M_.endo_names(21) = {'zcap'};
M_.endo_names_tex(21) = {'zcap'};
M_.endo_names_long(21) = {'zcap'};
M_.endo_names(22) = {'rk'};
M_.endo_names_tex(22) = {'rk'};
M_.endo_names_long(22) = {'rk'};
M_.endo_names(23) = {'k'};
M_.endo_names_tex(23) = {'k'};
M_.endo_names_long(23) = {'k'};
M_.endo_names(24) = {'pk'};
M_.endo_names_tex(24) = {'pk'};
M_.endo_names_long(24) = {'pk'};
M_.endo_names(25) = {'c'};
M_.endo_names_tex(25) = {'c'};
M_.endo_names_long(25) = {'c'};
M_.endo_names(26) = {'inve'};
M_.endo_names_tex(26) = {'inve'};
M_.endo_names_long(26) = {'inve'};
M_.endo_names(27) = {'y'};
M_.endo_names_tex(27) = {'y'};
M_.endo_names_long(27) = {'y'};
M_.endo_names(28) = {'lab'};
M_.endo_names_tex(28) = {'lab'};
M_.endo_names_long(28) = {'lab'};
M_.endo_names(29) = {'pinf'};
M_.endo_names_tex(29) = {'pinf'};
M_.endo_names_long(29) = {'pinf'};
M_.endo_names(30) = {'w'};
M_.endo_names_tex(30) = {'w'};
M_.endo_names_long(30) = {'w'};
M_.endo_names(31) = {'r'};
M_.endo_names_tex(31) = {'r'};
M_.endo_names_long(31) = {'r'};
M_.endo_names(32) = {'a'};
M_.endo_names_tex(32) = {'a'};
M_.endo_names_long(32) = {'a'};
M_.endo_names(33) = {'b'};
M_.endo_names_tex(33) = {'b'};
M_.endo_names_long(33) = {'b'};
M_.endo_names(34) = {'g'};
M_.endo_names_tex(34) = {'g'};
M_.endo_names_long(34) = {'g'};
M_.endo_names(35) = {'qs'};
M_.endo_names_tex(35) = {'qs'};
M_.endo_names_long(35) = {'qs'};
M_.endo_names(36) = {'ms'};
M_.endo_names_tex(36) = {'ms'};
M_.endo_names_long(36) = {'ms'};
M_.endo_names(37) = {'spinf'};
M_.endo_names_tex(37) = {'spinf'};
M_.endo_names_long(37) = {'spinf'};
M_.endo_names(38) = {'sw'};
M_.endo_names_tex(38) = {'sw'};
M_.endo_names_long(38) = {'sw'};
M_.endo_names(39) = {'kpf'};
M_.endo_names_tex(39) = {'kpf'};
M_.endo_names_long(39) = {'kpf'};
M_.endo_names(40) = {'kp'};
M_.endo_names_tex(40) = {'kp'};
M_.endo_names_long(40) = {'kp'};
M_.endo_names(41) = {'pinf4'};
M_.endo_names_tex(41) = {'pinf4'};
M_.endo_names_long(41) = {'pinf4'};
M_.endo_names(42) = {'AUX_ENDO_LAG_28_1'};
M_.endo_names_tex(42) = {'AUX\_ENDO\_LAG\_28\_1'};
M_.endo_names_long(42) = {'AUX_ENDO_LAG_28_1'};
M_.endo_names(43) = {'AUX_ENDO_LAG_28_2'};
M_.endo_names_tex(43) = {'AUX\_ENDO\_LAG\_28\_2'};
M_.endo_names_long(43) = {'AUX_ENDO_LAG_28_2'};
M_.endo_partitions = struct();
M_.param_names = cell(54,1);
M_.param_names_tex = cell(54,1);
M_.param_names_long = cell(54,1);
M_.param_names(1) = {'curvw'};
M_.param_names_tex(1) = {'curvw'};
M_.param_names_long(1) = {'curvw'};
M_.param_names(2) = {'cgy'};
M_.param_names_tex(2) = {'cgy'};
M_.param_names_long(2) = {'cgy'};
M_.param_names(3) = {'curvp'};
M_.param_names_tex(3) = {'curvp'};
M_.param_names_long(3) = {'curvp'};
M_.param_names(4) = {'constelab'};
M_.param_names_tex(4) = {'constelab'};
M_.param_names_long(4) = {'constelab'};
M_.param_names(5) = {'constepinf'};
M_.param_names_tex(5) = {'constepinf'};
M_.param_names_long(5) = {'constepinf'};
M_.param_names(6) = {'constebeta'};
M_.param_names_tex(6) = {'constebeta'};
M_.param_names_long(6) = {'constebeta'};
M_.param_names(7) = {'cmaw'};
M_.param_names_tex(7) = {'cmaw'};
M_.param_names_long(7) = {'cmaw'};
M_.param_names(8) = {'cmap'};
M_.param_names_tex(8) = {'cmap'};
M_.param_names_long(8) = {'cmap'};
M_.param_names(9) = {'calfa'};
M_.param_names_tex(9) = {'calfa'};
M_.param_names_long(9) = {'calfa'};
M_.param_names(10) = {'czcap'};
M_.param_names_tex(10) = {'czcap'};
M_.param_names_long(10) = {'czcap'};
M_.param_names(11) = {'cbeta'};
M_.param_names_tex(11) = {'cbeta'};
M_.param_names_long(11) = {'cbeta'};
M_.param_names(12) = {'csadjcost'};
M_.param_names_tex(12) = {'csadjcost'};
M_.param_names_long(12) = {'csadjcost'};
M_.param_names(13) = {'ctou'};
M_.param_names_tex(13) = {'ctou'};
M_.param_names_long(13) = {'ctou'};
M_.param_names(14) = {'csigma'};
M_.param_names_tex(14) = {'csigma'};
M_.param_names_long(14) = {'csigma'};
M_.param_names(15) = {'chabb'};
M_.param_names_tex(15) = {'chabb'};
M_.param_names_long(15) = {'chabb'};
M_.param_names(16) = {'ccs'};
M_.param_names_tex(16) = {'ccs'};
M_.param_names_long(16) = {'ccs'};
M_.param_names(17) = {'cinvs'};
M_.param_names_tex(17) = {'cinvs'};
M_.param_names_long(17) = {'cinvs'};
M_.param_names(18) = {'cfc'};
M_.param_names_tex(18) = {'cfc'};
M_.param_names_long(18) = {'cfc'};
M_.param_names(19) = {'cindw'};
M_.param_names_tex(19) = {'cindw'};
M_.param_names_long(19) = {'cindw'};
M_.param_names(20) = {'cprobw'};
M_.param_names_tex(20) = {'cprobw'};
M_.param_names_long(20) = {'cprobw'};
M_.param_names(21) = {'cindp'};
M_.param_names_tex(21) = {'cindp'};
M_.param_names_long(21) = {'cindp'};
M_.param_names(22) = {'cprobp'};
M_.param_names_tex(22) = {'cprobp'};
M_.param_names_long(22) = {'cprobp'};
M_.param_names(23) = {'csigl'};
M_.param_names_tex(23) = {'csigl'};
M_.param_names_long(23) = {'csigl'};
M_.param_names(24) = {'clandaw'};
M_.param_names_tex(24) = {'clandaw'};
M_.param_names_long(24) = {'clandaw'};
M_.param_names(25) = {'crdpi'};
M_.param_names_tex(25) = {'crdpi'};
M_.param_names_long(25) = {'crdpi'};
M_.param_names(26) = {'crpi'};
M_.param_names_tex(26) = {'crpi'};
M_.param_names_long(26) = {'crpi'};
M_.param_names(27) = {'crdy'};
M_.param_names_tex(27) = {'crdy'};
M_.param_names_long(27) = {'crdy'};
M_.param_names(28) = {'cry'};
M_.param_names_tex(28) = {'cry'};
M_.param_names_long(28) = {'cry'};
M_.param_names(29) = {'crr'};
M_.param_names_tex(29) = {'crr'};
M_.param_names_long(29) = {'crr'};
M_.param_names(30) = {'crhoa'};
M_.param_names_tex(30) = {'crhoa'};
M_.param_names_long(30) = {'crhoa'};
M_.param_names(31) = {'crhob'};
M_.param_names_tex(31) = {'crhob'};
M_.param_names_long(31) = {'crhob'};
M_.param_names(32) = {'crhog'};
M_.param_names_tex(32) = {'crhog'};
M_.param_names_long(32) = {'crhog'};
M_.param_names(33) = {'crhoqs'};
M_.param_names_tex(33) = {'crhoqs'};
M_.param_names_long(33) = {'crhoqs'};
M_.param_names(34) = {'crhoms'};
M_.param_names_tex(34) = {'crhoms'};
M_.param_names_long(34) = {'crhoms'};
M_.param_names(35) = {'crhopinf'};
M_.param_names_tex(35) = {'crhopinf'};
M_.param_names_long(35) = {'crhopinf'};
M_.param_names(36) = {'crhow'};
M_.param_names_tex(36) = {'crhow'};
M_.param_names_long(36) = {'crhow'};
M_.param_names(37) = {'ctrend'};
M_.param_names_tex(37) = {'ctrend'};
M_.param_names_long(37) = {'ctrend'};
M_.param_names(38) = {'cg'};
M_.param_names_tex(38) = {'cg'};
M_.param_names_long(38) = {'cg'};
M_.param_names(39) = {'cgamma'};
M_.param_names_tex(39) = {'cgamma'};
M_.param_names_long(39) = {'cgamma'};
M_.param_names(40) = {'clandap'};
M_.param_names_tex(40) = {'clandap'};
M_.param_names_long(40) = {'clandap'};
M_.param_names(41) = {'cbetabar'};
M_.param_names_tex(41) = {'cbetabar'};
M_.param_names_long(41) = {'cbetabar'};
M_.param_names(42) = {'cr'};
M_.param_names_tex(42) = {'cr'};
M_.param_names_long(42) = {'cr'};
M_.param_names(43) = {'cpie'};
M_.param_names_tex(43) = {'cpie'};
M_.param_names_long(43) = {'cpie'};
M_.param_names(44) = {'crk'};
M_.param_names_tex(44) = {'crk'};
M_.param_names_long(44) = {'crk'};
M_.param_names(45) = {'cw'};
M_.param_names_tex(45) = {'cw'};
M_.param_names_long(45) = {'cw'};
M_.param_names(46) = {'cikbar'};
M_.param_names_tex(46) = {'cikbar'};
M_.param_names_long(46) = {'cikbar'};
M_.param_names(47) = {'cik'};
M_.param_names_tex(47) = {'cik'};
M_.param_names_long(47) = {'cik'};
M_.param_names(48) = {'clk'};
M_.param_names_tex(48) = {'clk'};
M_.param_names_long(48) = {'clk'};
M_.param_names(49) = {'cky'};
M_.param_names_tex(49) = {'cky'};
M_.param_names_long(49) = {'cky'};
M_.param_names(50) = {'ciy'};
M_.param_names_tex(50) = {'ciy'};
M_.param_names_long(50) = {'ciy'};
M_.param_names(51) = {'ccy'};
M_.param_names_tex(51) = {'ccy'};
M_.param_names_long(51) = {'ccy'};
M_.param_names(52) = {'crkky'};
M_.param_names_tex(52) = {'crkky'};
M_.param_names_long(52) = {'crkky'};
M_.param_names(53) = {'cwhlc'};
M_.param_names_tex(53) = {'cwhlc'};
M_.param_names_long(53) = {'cwhlc'};
M_.param_names(54) = {'cwly'};
M_.param_names_tex(54) = {'cwly'};
M_.param_names_long(54) = {'cwly'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 7;
M_.endo_nbr = 43;
M_.param_nbr = 54;
M_.orig_endo_nbr = 41;
M_.aux_vars(1).endo_index = 42;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 29;
M_.aux_vars(1).orig_lead_lag = -1;
M_.aux_vars(1).orig_expr = 'pinf(-1)';
M_.aux_vars(2).endo_index = 43;
M_.aux_vars(2).type = 1;
M_.aux_vars(2).orig_index = 29;
M_.aux_vars(2).orig_lead_lag = -2;
M_.aux_vars(2).orig_expr = 'AUX_ENDO_LAG_28_1(-1)';
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(7, 7);
M_.Correlation_matrix = eye(7, 7);
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
M_.orig_eq_nbr = 41;
M_.eq_nbr = 43;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 3;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 3;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 3;
M_.lead_lag_incidence = [
 0 23 0;
 0 24 0;
 0 25 0;
 0 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 1 30 0;
 2 31 0;
 0 32 0;
 0 33 66;
 0 34 0;
 0 35 67;
 3 36 68;
 4 37 69;
 5 38 0;
 0 39 70;
 0 40 0;
 0 41 0;
 0 42 0;
 0 43 0;
 0 44 71;
 0 45 0;
 0 46 72;
 6 47 73;
 7 48 74;
 8 49 0;
 0 50 75;
 9 51 76;
 10 52 77;
 11 53 0;
 12 54 0;
 13 55 0;
 14 56 0;
 15 57 0;
 16 58 0;
 17 59 0;
 18 60 0;
 19 61 0;
 20 62 0;
 0 63 0;
 21 64 0;
 22 65 0;]';
M_.nstatic = 15;
M_.nfwrd   = 6;
M_.npred   = 16;
M_.nboth   = 6;
M_.nsfwrd   = 12;
M_.nspred   = 22;
M_.ndynamic   = 28;
M_.dynamic_tmp_nbr = [14; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'a' ;
  2 , 'name' , 'zcapf' ;
  3 , 'name' , 'rkf' ;
  4 , 'name' , 'kf' ;
  5 , 'name' , 'invef' ;
  6 , 'name' , 'pkf' ;
  7 , 'name' , 'cf' ;
  8 , 'name' , 'yf' ;
  9 , 'name' , '9' ;
  10 , 'name' , 'wf' ;
  11 , 'name' , 'kpf' ;
  12 , 'name' , 'mc' ;
  13 , 'name' , 'zcap' ;
  14 , 'name' , 'rk' ;
  15 , 'name' , 'k' ;
  16 , 'name' , 'inve' ;
  17 , 'name' , 'pk' ;
  18 , 'name' , 'c' ;
  19 , 'name' , 'y' ;
  20 , 'name' , '20' ;
  21 , 'name' , 'pinf' ;
  22 , 'name' , 'w' ;
  23 , 'name' , 'r' ;
  24 , 'name' , '24' ;
  25 , 'name' , 'b' ;
  26 , 'name' , 'g' ;
  27 , 'name' , 'qs' ;
  28 , 'name' , 'ms' ;
  29 , 'name' , 'spinf' ;
  30 , 'name' , 'epinfma' ;
  31 , 'name' , 'sw' ;
  32 , 'name' , 'ewma' ;
  33 , 'name' , 'kp' ;
  34 , 'name' , 'dy' ;
  35 , 'name' , 'dc' ;
  36 , 'name' , 'dinve' ;
  37 , 'name' , 'dw' ;
  38 , 'name' , 'pinfobs' ;
  39 , 'name' , 'pinf4' ;
  40 , 'name' , 'robs' ;
  41 , 'name' , 'labobs' ;
};
M_.mapping.labobs.eqidx = [41 ];
M_.mapping.robs.eqidx = [40 ];
M_.mapping.pinfobs.eqidx = [38 ];
M_.mapping.dy.eqidx = [34 ];
M_.mapping.dc.eqidx = [35 ];
M_.mapping.dinve.eqidx = [36 ];
M_.mapping.dw.eqidx = [37 ];
M_.mapping.ewma.eqidx = [31 32 ];
M_.mapping.epinfma.eqidx = [29 30 ];
M_.mapping.zcapf.eqidx = [2 4 8 ];
M_.mapping.rkf.eqidx = [1 2 3 6 ];
M_.mapping.kf.eqidx = [3 4 9 ];
M_.mapping.pkf.eqidx = [5 6 ];
M_.mapping.cf.eqidx = [7 8 10 ];
M_.mapping.invef.eqidx = [5 8 11 ];
M_.mapping.yf.eqidx = [8 9 23 ];
M_.mapping.labf.eqidx = [3 7 9 10 ];
M_.mapping.wf.eqidx = [1 3 10 ];
M_.mapping.rrf.eqidx = [6 7 ];
M_.mapping.mc.eqidx = [12 21 ];
M_.mapping.zcap.eqidx = [13 15 19 ];
M_.mapping.rk.eqidx = [12 13 14 17 ];
M_.mapping.k.eqidx = [14 15 20 ];
M_.mapping.pk.eqidx = [16 17 ];
M_.mapping.c.eqidx = [18 19 22 35 ];
M_.mapping.inve.eqidx = [16 19 33 36 ];
M_.mapping.y.eqidx = [19 20 23 34 ];
M_.mapping.lab.eqidx = [14 18 20 22 41 ];
M_.mapping.pinf.eqidx = [17 18 21 22 23 38 39 ];
M_.mapping.w.eqidx = [12 14 22 37 ];
M_.mapping.r.eqidx = [17 18 23 40 ];
M_.mapping.a.eqidx = [1 9 12 20 24 ];
M_.mapping.b.eqidx = [6 7 17 18 25 ];
M_.mapping.g.eqidx = [8 19 26 ];
M_.mapping.qs.eqidx = [5 11 16 27 33 ];
M_.mapping.ms.eqidx = [23 28 ];
M_.mapping.spinf.eqidx = [21 29 ];
M_.mapping.sw.eqidx = [22 31 ];
M_.mapping.kpf.eqidx = [4 11 ];
M_.mapping.kp.eqidx = [15 33 ];
M_.mapping.pinf4.eqidx = [39 ];
M_.mapping.ea.eqidx = [24 26 ];
M_.mapping.eb.eqidx = [25 ];
M_.mapping.eqs.eqidx = [27 ];
M_.mapping.eg.eqidx = [26 ];
M_.mapping.em.eqidx = [28 ];
M_.mapping.epinf.eqidx = [30 ];
M_.mapping.ew.eqidx = [32 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [8 9 14 15 16 25 26 27 29 30 31 32 33 34 35 36 37 38 39 40 42 43 ];
M_.exo_names_orig_ord = [1:7];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(43, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(7, 1);
M_.params = NaN(54, 1);
M_.endo_trends = struct('deflator', cell(43, 1), 'log_deflator', cell(43, 1), 'growth_factor', cell(43, 1), 'log_growth_factor', cell(43, 1));
M_.NNZDerivatives = [169; 0; -1; ];
M_.static_tmp_nbr = [11; 2; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(13) = .025;
ctou = M_.params(13);
M_.params(24) = 1.5;
clandaw = M_.params(24);
M_.params(38) = 0.18;
cg = M_.params(38);
M_.params(3) = 10;
curvp = M_.params(3);
M_.params(1) = 10;
curvw = M_.params(1);
M_.params(37) = 0.4312;
ctrend = M_.params(37);
M_.params(39) = 1+M_.params(37)/100;
cgamma = M_.params(39);
M_.params(6) = 0.1657;
constebeta = M_.params(6);
M_.params(11) = 100/(100+M_.params(6));
cbeta = M_.params(11);
M_.params(5) = 0.7869;
constepinf = M_.params(5);
M_.params(43) = 1+M_.params(5)/100;
cpie = M_.params(43);
M_.params(4) = 0.5509;
constelab = M_.params(4);
M_.params(9) = 0.1901;
calfa = M_.params(9);
M_.params(14) = 1.3808;
csigma = M_.params(14);
M_.params(18) = 1.6064;
cfc = M_.params(18);
M_.params(2) = 0.5187;
cgy = M_.params(2);
M_.params(12) = 5.7606;
csadjcost = M_.params(12);
M_.params(15) = 0.7133;
chabb = M_.params(15);
M_.params(20) = 0.7061;
cprobw = M_.params(20);
M_.params(23) = 1.8383;
csigl = M_.params(23);
M_.params(22) = 0.6523;
cprobp = M_.params(22);
M_.params(19) = 0.5845;
cindw = M_.params(19);
M_.params(21) = 0.2432;
cindp = M_.params(21);
M_.params(10) = 0.5462;
czcap = M_.params(10);
M_.params(26) = 2.0443;
crpi = M_.params(26);
M_.params(29) = 0.8103;
crr = M_.params(29);
M_.params(28) = 0.0882;
cry = M_.params(28);
M_.params(27) = 0.2247;
crdy = M_.params(27);
M_.params(30) = 0.9577;
crhoa = M_.params(30);
M_.params(31) = 0.2194;
crhob = M_.params(31);
M_.params(32) = 0.9767;
crhog = M_.params(32);
M_.params(33) = 0.7113;
crhoqs = M_.params(33);
M_.params(34) = 0.1479;
crhoms = M_.params(34);
M_.params(35) = 0.8895;
crhopinf = M_.params(35);
M_.params(36) = 0.9688;
crhow = M_.params(36);
M_.params(8) = 0.7010;
cmap = M_.params(8);
M_.params(7) = 0.8503;
cmaw = M_.params(7);
M_.params(40) = M_.params(18);
clandap = M_.params(40);
M_.params(41) = M_.params(11)*M_.params(39)^(-M_.params(14));
cbetabar = M_.params(41);
M_.params(42) = M_.params(43)/(M_.params(11)*M_.params(39)^(-M_.params(14)));
cr = M_.params(42);
M_.params(44) = M_.params(11)^(-1)*M_.params(39)^M_.params(14)-(1-M_.params(13));
crk = M_.params(44);
M_.params(45) = (M_.params(9)^M_.params(9)*(1-M_.params(9))^(1-M_.params(9))/(M_.params(40)*M_.params(44)^M_.params(9)))^(1/(1-M_.params(9)));
cw = M_.params(45);
M_.params(46) = 1-(1-M_.params(13))/M_.params(39);
cikbar = M_.params(46);
M_.params(47) = M_.params(39)*(1-(1-M_.params(13))/M_.params(39));
cik = M_.params(47);
M_.params(48) = (1-M_.params(9))/M_.params(9)*M_.params(44)/M_.params(45);
clk = M_.params(48);
M_.params(49) = M_.params(18)*M_.params(48)^(M_.params(9)-1);
cky = M_.params(49);
M_.params(50) = M_.params(47)*M_.params(49);
ciy = M_.params(50);
M_.params(51) = 1-M_.params(38)-M_.params(47)*M_.params(49);
ccy = M_.params(51);
M_.params(52) = M_.params(44)*M_.params(49);
crkky = M_.params(52);
M_.params(53) = M_.params(49)*M_.params(44)*(1-M_.params(9))*1/M_.params(24)/M_.params(9)/M_.params(51);
cwhlc = M_.params(53);
M_.params(54) = 1-M_.params(44)*M_.params(49);
cwly = M_.params(54);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.4582)^2;
M_.Sigma_e(2, 2) = (0.2400)^2;
M_.Sigma_e(3, 3) = (0.4526)^2;
M_.Sigma_e(4, 4) = (0.5291)^2;
M_.Sigma_e(5, 5) = (0.2449)^2;
M_.Sigma_e(6, 6) = (0.1410)^2;
M_.Sigma_e(7, 7) = (0.2446)^2;
options_.irf = 20;
options_.nograph = true;
options_.noprint = true;
options_.order = 2;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
labobs_index=find(strcmp(M_.endo_names,'labobs'));
robs_index=find(strcmp(M_.endo_names,'robs'));
pinfobs_index=find(strcmp(M_.endo_names,'pinfobs'));
dy_index=find(strcmp(M_.endo_names,'dy'));
dc_index=find(strcmp(M_.endo_names,'dc'));
dinve_index=find(strcmp(M_.endo_names,'dinve'));
dw_index=find(strcmp(M_.endo_names,'dw'));
std_L=(oo_.var(labobs_index,labobs_index))^(1/2);
std_R=(oo_.var(robs_index,robs_index))^(1/2);
std_pi=(oo_.var(pinfobs_index,pinfobs_index))^(1/2);
std_Y=(oo_.var(dy_index,dy_index))^(1/2);
std_C=(oo_.var(dc_index,dc_index))^(1/2);
std_I=(oo_.var(dinve_index,dinve_index))^(1/2);
std_w=(oo_.var(dw_index,dw_index))^(1/2);
targets=[2.908 0.830 0.615 0.856 0.685 2.256 0.564];
actual_dynare=[std_L std_R std_pi std_Y std_C std_I std_w];
create_matrix_quadratic_from_dynare
ALPHA_ZS_dynare=[zeros(n,nstatic) oo_.dr.ghx zeros(n,nfwrd)];
ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
matrix_quadratic.X=ALPHA_ZS_dynare;
BETA_ZS_dynare=oo_.dr.ghu; BETA_ZS_dynare=BETA_ZS_dynare(oo_.dr.inv_order_var,:);
matrix_quadratic.P=ALPHA_ZS_dynare;
matrix_quadratic.Q=BETA_ZS_dynare;
[errors_dynare] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_dynare=matrix_quadratic;
save('dynare_results','targets','actual_dynare','matrix_quadratic_dynare','errors_dynare')
options_.dr_cycle_reduction = true;
options_.dr_cycle_reduction_tol = 1e-16;
options_.irf = 0;
options_.noprint = true;
options_.order = 1;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
std_L=(oo_.var(labobs_index,labobs_index))^(1/2);
std_R=(oo_.var(robs_index,robs_index))^(1/2);
std_pi=(oo_.var(pinfobs_index,pinfobs_index))^(1/2);
std_Y=(oo_.var(dy_index,dy_index))^(1/2);
std_C=(oo_.var(dc_index,dc_index))^(1/2);
std_I=(oo_.var(dinve_index,dinve_index))^(1/2);
std_w=(oo_.var(dw_index,dw_index))^(1/2);
actual_dynare_cr=[std_L std_R std_pi std_Y std_C std_I std_w];
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
save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Smets_Wouters_2007_rep_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
