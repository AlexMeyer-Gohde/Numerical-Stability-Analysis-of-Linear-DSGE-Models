function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 13);

T = jermann98_rep_def_R_rr_compare.static_resid_tt(T, y, x, params);

T(9) = exp(params(5)+y(12))*(-(exp(y(1))*exp(y(3))))/(exp(y(1))*exp(y(1)));
T(10) = getPowerDeriv(T(6),1-params(9),1);
T(11) = getPowerDeriv(T(6),(-params(9)),1);
T(12) = getPowerDeriv(exp(y(2))-exp(y(2))*params(2)*exp((-(params(5)+y(12)))),(-params(1)),1);
T(13) = getPowerDeriv(exp(y(2))*exp(params(5)+y(12))-exp(y(2))*params(2),(-params(1)),1);

end
