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
    T = jermann98_rep_def_R_rr_compare.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(17, 1);
lhs = y(16);
rhs = (exp(y(9))-params(2)*exp(y(2))*exp((-(params(5)+y(19)))))^(-params(1))-params(2)*T(1)*(exp(y(25))*exp(params(5)+y(27))-exp(y(9))*params(2))^(-params(1));
residual(1) = lhs - rhs;
lhs = exp(y(24));
rhs = T(1)*exp((params(5)+y(19))*(-params(1)))*y(16)/y(5);
residual(2) = lhs - rhs;
lhs = exp(y(8));
rhs = (1-params(10)+params(12)*y(17))*exp((-params(5))-y(19))*exp(y(1))+(1-params(12))*exp(y(10));
residual(3) = lhs - rhs;
lhs = exp(y(11));
rhs = T(3)*T(4);
residual(4) = lhs - rhs;
residual(5) = 1-exp(y(12));
lhs = exp(y(11));
rhs = exp(y(9))+exp(y(10));
residual(6) = lhs - rhs;
lhs = y(17);
rhs = T(5)/(1-params(9))*T(6)^(1-params(9))+exp(params(5))*(1-(1-params(10))*exp((-params(5))))*(-params(9))/(1-params(9));
residual(7) = lhs - rhs;
lhs = y(18);
rhs = T(5)*T(6)^(-params(9));
residual(8) = lhs - rhs;
lhs = exp(y(21));
rhs = exp(y(11))*(1-params(4))/exp(y(12));
residual(9) = lhs - rhs;
lhs = exp(y(22));
rhs = exp(y(11))-exp(y(12))*exp(y(21))-exp(y(10));
residual(10) = lhs - rhs;
lhs = exp(y(23));
rhs = T(8)*(1+params(12)*y(6)-params(12));
residual(11) = lhs - rhs;
lhs = T(9);
rhs = 1;
residual(12) = lhs - rhs;
lhs = y(13);
rhs = y(19)+params(5)+y(9)-y(2);
residual(13) = lhs - rhs;
lhs = y(14);
rhs = y(19)+params(5)+y(11)-y(4);
residual(14) = lhs - rhs;
lhs = y(15);
rhs = y(19)+params(5)+y(10)-y(3);
residual(15) = lhs - rhs;
lhs = y(19);
rhs = params(7)*x(it_, 2);
residual(16) = lhs - rhs;
lhs = y(20);
rhs = params(6)*y(7)+params(8)*x(it_, 1);
residual(17) = lhs - rhs;

end
