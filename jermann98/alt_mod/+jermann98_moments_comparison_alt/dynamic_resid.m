function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = jermann98_moments_comparison_alt.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(13, 1);
lhs = y(15);
rhs = (exp(y(9))-params(2)*exp(y(2))*exp((-params(6))))^(-params(1))-params(2)*params(3)*(exp(y(21))*exp(params(6))-exp(y(9))*params(2))^(-params(1));
residual(1) = lhs - rhs;
lhs = exp(y(20));
rhs = params(3)*exp(params(6)*(-params(1)))*y(15)/y(5);
residual(2) = lhs - rhs;
lhs = exp(y(8));
rhs = exp((-params(6)))*(1-params(13)+params(16)*y(16))*exp(y(1))+(1-params(16))*exp(y(10));
residual(3) = lhs - rhs;
lhs = exp(y(11));
rhs = T(2);
residual(4) = lhs - rhs;
lhs = exp(y(11));
rhs = exp(y(9))+exp(y(10));
residual(5) = lhs - rhs;
lhs = y(16);
rhs = params(11)/(1-params(10))*T(3)^(1-params(10))+params(12);
residual(6) = lhs - rhs;
lhs = y(17);
rhs = params(11)*T(3)^(-params(10));
residual(7) = lhs - rhs;
lhs = exp(y(19));
rhs = T(4)*(1+params(16)*y(6)-params(16));
residual(8) = lhs - rhs;
lhs = exp(y(22)+y(23));
rhs = 1;
residual(9) = lhs - rhs;
lhs = y(12);
rhs = params(6)+y(9)-y(2);
residual(10) = lhs - rhs;
lhs = y(13);
rhs = params(6)+y(11)-y(4);
residual(11) = lhs - rhs;
lhs = y(14);
rhs = params(6)+y(10)-y(3);
residual(12) = lhs - rhs;
lhs = y(18);
rhs = params(7)*y(7)+params(9)*x(it_, 1);
residual(13) = lhs - rhs;

end
