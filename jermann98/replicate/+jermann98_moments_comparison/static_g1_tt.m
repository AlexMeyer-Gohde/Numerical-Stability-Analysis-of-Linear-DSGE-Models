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

assert(length(T) >= 7);

T = jermann98_moments_comparison.static_resid_tt(T, y, x, params);

T(5) = exp(params(6))*(-(exp(y(1))*exp(y(3))))/(exp(y(1))*exp(y(1)));
T(6) = getPowerDeriv(T(3),1-params(10),1);
T(7) = getPowerDeriv(T(3),(-params(10)),1);

end
