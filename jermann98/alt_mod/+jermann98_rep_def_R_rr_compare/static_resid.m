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
    T = jermann98_rep_def_R_rr_compare.static_resid_tt(T, y, x, params);
end
residual = zeros(17, 1);
lhs = y(9);
rhs = (exp(y(2))-exp(y(2))*params(2)*exp((-(params(5)+y(12)))))^(-params(1))-params(2)*T(1)*(exp(y(2))*exp(params(5)+y(12))-exp(y(2))*params(2))^(-params(1));
residual(1) = lhs - rhs;
lhs = exp(y(17));
rhs = T(1)*exp((params(5)+y(12))*(-params(1)));
residual(2) = lhs - rhs;
lhs = exp(y(1));
rhs = exp(y(1))*(1-params(10)+params(12)*y(10))*exp((-params(5))-y(12))+(1-params(12))*exp(y(3));
residual(3) = lhs - rhs;
lhs = exp(y(4));
rhs = T(3)*T(4);
residual(4) = lhs - rhs;
residual(5) = 1-exp(y(5));
lhs = exp(y(4));
rhs = exp(y(2))+exp(y(3));
residual(6) = lhs - rhs;
lhs = y(10);
rhs = T(5)/(1-params(9))*T(6)^(1-params(9))+exp(params(5))*(1-(1-params(10))*exp((-params(5))))*(-params(9))/(1-params(9));
residual(7) = lhs - rhs;
lhs = y(11);
rhs = T(5)*T(6)^(-params(9));
residual(8) = lhs - rhs;
lhs = exp(y(14));
rhs = exp(y(4))*(1-params(4))/exp(y(5));
residual(9) = lhs - rhs;
lhs = exp(y(15));
rhs = exp(y(4))-exp(y(5))*exp(y(14))-exp(y(3));
residual(10) = lhs - rhs;
lhs = exp(y(16));
rhs = (1+params(12)*y(11)-params(12))*T(8);
residual(11) = lhs - rhs;
lhs = exp((params(5)+y(12))*(-params(1)))*T(1)*exp(y(16));
rhs = 1;
residual(12) = lhs - rhs;
lhs = y(6);
rhs = params(5)+y(12);
residual(13) = lhs - rhs;
lhs = y(7);
rhs = params(5)+y(12);
residual(14) = lhs - rhs;
lhs = y(8);
rhs = params(5)+y(12);
residual(15) = lhs - rhs;
lhs = y(12);
rhs = params(7)*x(2);
residual(16) = lhs - rhs;
lhs = y(13);
rhs = y(13)*params(6)+params(8)*x(1);
residual(17) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
