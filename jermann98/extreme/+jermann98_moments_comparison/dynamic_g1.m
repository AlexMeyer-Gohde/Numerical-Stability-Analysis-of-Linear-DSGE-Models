function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = jermann98_moments_comparison.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(13, 24);
g1(1,2)=(-((-(params(2)*exp(y(2))*exp((-params(6)))))*T(8)));
g1(1,9)=(-(exp(y(9))*T(8)-params(2)*params(3)*(-(exp(y(9))*params(2)))*T(9)));
g1(1,21)=params(2)*params(3)*exp(y(21))*exp(params(6))*T(9);
g1(1,15)=1;
g1(2,5)=(-(params(3)*exp(params(6)*(-params(1)))*(-y(15))/(y(5)*y(5))));
g1(2,15)=(-(params(3)*exp(params(6)*(-params(1)))*1/y(5)));
g1(2,20)=exp(y(20));
g1(3,1)=(-(exp((-params(6)))*(1-params(13)+params(16)*y(16))*exp(y(1))));
g1(3,8)=exp(y(8));
g1(3,10)=(-((1-params(16))*exp(y(10))));
g1(3,16)=(-(exp(y(1))*exp((-params(6)))*params(16)));
g1(4,1)=(-(T(1)*exp(y(18))*exp(params(6)*(-params(5)))*exp(y(1))*getPowerDeriv(exp(y(1)),params(5),1)));
g1(4,11)=exp(y(11));
g1(4,18)=(-T(2));
g1(5,9)=(-exp(y(9)));
g1(5,10)=(-exp(y(10)));
g1(5,11)=exp(y(11));
g1(6,1)=(-(params(11)/(1-params(10))*T(5)*T(6)));
g1(6,10)=(-(params(11)/(1-params(10))*T(3)*T(6)));
g1(6,16)=1;
g1(7,1)=(-(params(11)*T(5)*T(7)));
g1(7,10)=(-(params(11)*T(3)*T(7)));
g1(7,17)=1;
g1(8,1)=(-((1+params(16)*y(6)-params(16))*(exp(params(6))*params(5)*(-(exp(y(1))*exp(y(11))))/(exp(y(1))*exp(y(1)))-params(16)*T(5))));
g1(8,10)=(-((1+params(16)*y(6)-params(16))*(-(params(16)*T(3)))));
g1(8,11)=(-(exp(params(6))*params(5)*exp(y(11))/exp(y(1))*(1+params(16)*y(6)-params(16))));
g1(8,16)=(-((1+params(16)*y(6)-params(16))*params(16)/(1+params(16)*y(17)-params(16))));
g1(8,6)=(-(params(16)*T(4)));
g1(8,17)=(-((1+params(16)*y(6)-params(16))*(-(params(16)*(1+params(16)*y(16)-params(13))))/((1+params(16)*y(17)-params(16))*(1+params(16)*y(17)-params(16)))));
g1(8,19)=exp(y(19));
g1(9,15)=exp(params(6)*(-params(1)))*params(3)*exp(y(23))*(-y(22))/(y(15)*y(15));
g1(9,22)=exp(params(6)*(-params(1)))*params(3)*exp(y(23))*1/y(15);
g1(9,23)=exp(params(6)*(-params(1)))*params(3)*exp(y(23))*y(22)/y(15);
g1(10,2)=1;
g1(10,9)=(-1);
g1(10,12)=1;
g1(11,4)=1;
g1(11,11)=(-1);
g1(11,13)=1;
g1(12,3)=1;
g1(12,10)=(-1);
g1(12,14)=1;
g1(13,7)=(-params(7));
g1(13,18)=1;
g1(13,24)=(-params(9));

end
