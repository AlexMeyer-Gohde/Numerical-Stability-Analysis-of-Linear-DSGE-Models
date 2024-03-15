function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
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

assert(length(T) >= 8);

T(1) = params(3)*exp(params(5))^(params(1)-1);
T(2) = exp(y(1))^params(4);
T(3) = exp(y(13))*exp((params(5)+y(12))*(-params(4)))*T(2);
T(4) = exp(y(5))^(1-params(4));
T(5) = (exp(params(5))*(1-(1-params(10))*exp((-params(5)))))^params(9);
T(6) = exp(params(5)+y(12))*exp(y(3))/exp(y(1));
T(7) = exp(params(5)+y(12))*params(4)*exp(y(4))/exp(y(1));
T(8) = T(7)+(1+params(12)*y(10)-params(10))/(1+params(12)*y(11)-params(12))-params(12)*T(6);

end
