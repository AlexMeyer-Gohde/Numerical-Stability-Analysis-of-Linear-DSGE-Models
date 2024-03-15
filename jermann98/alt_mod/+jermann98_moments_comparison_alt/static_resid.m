function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = jermann98_moments_comparison_alt.static_resid_tt(T, y, x, params);
end
residual = zeros(13, 1);
lhs = y(8);
rhs = (exp(y(2))-exp(y(2))*params(2)*exp((-params(6))))^(-params(1))-params(2)*params(3)*(exp(y(2))*exp(params(6))-exp(y(2))*params(2))^(-params(1));
residual(1) = lhs - rhs;
lhs = exp(y(13));
rhs = params(3)*exp(params(6)*(-params(1)));
residual(2) = lhs - rhs;
lhs = exp(y(1));
rhs = exp(y(1))*exp((-params(6)))*(1-params(13)+params(16)*y(9))+(1-params(16))*exp(y(3));
residual(3) = lhs - rhs;
lhs = exp(y(4));
rhs = T(2);
residual(4) = lhs - rhs;
lhs = exp(y(4));
rhs = exp(y(2))+exp(y(3));
residual(5) = lhs - rhs;
lhs = y(9);
rhs = params(11)/(1-params(10))*T(3)^(1-params(10))+params(12);
residual(6) = lhs - rhs;
lhs = y(10);
rhs = params(11)*T(3)^(-params(10));
residual(7) = lhs - rhs;
lhs = exp(y(12));
rhs = (1+params(16)*y(10)-params(16))*T(4);
residual(8) = lhs - rhs;
lhs = exp(y(13)+y(12));
rhs = 1;
residual(9) = lhs - rhs;
lhs = y(5);
rhs = params(6);
residual(10) = lhs - rhs;
lhs = y(6);
rhs = params(6);
residual(11) = lhs - rhs;
lhs = y(7);
rhs = params(6);
residual(12) = lhs - rhs;
lhs = y(11);
rhs = y(11)*params(7)+params(9)*x(1);
residual(13) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
