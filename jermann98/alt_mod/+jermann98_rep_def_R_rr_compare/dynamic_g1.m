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
    T = jermann98_rep_def_R_rr_compare.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(17, 30);
g1(1,2)=(-((-(params(2)*exp(y(2))*exp((-(params(5)+y(19))))))*T(13)));
g1(1,9)=(-(exp(y(9))*T(13)-params(2)*T(1)*(-(exp(y(9))*params(2)))*T(14)));
g1(1,25)=params(2)*T(1)*exp(y(25))*exp(params(5)+y(27))*T(14);
g1(1,16)=1;
g1(1,19)=(-(T(13)*(-(params(2)*exp(y(2))*(-exp((-(params(5)+y(19)))))))));
g1(1,27)=params(2)*T(1)*exp(y(25))*exp(params(5)+y(27))*T(14);
g1(2,5)=(-(T(1)*exp((params(5)+y(19))*(-params(1)))*(-y(16))/(y(5)*y(5))));
g1(2,16)=(-(T(1)*exp((params(5)+y(19))*(-params(1)))*1/y(5)));
g1(2,19)=(-(y(16)/y(5)*T(1)*(-params(1))*exp((params(5)+y(19))*(-params(1)))));
g1(2,24)=exp(y(24));
g1(3,1)=(-((1-params(10)+params(12)*y(17))*exp((-params(5))-y(19))*exp(y(1))));
g1(3,8)=exp(y(8));
g1(3,10)=(-((1-params(12))*exp(y(10))));
g1(3,17)=(-(exp(y(1))*params(12)*exp((-params(5))-y(19))));
g1(3,19)=(-(exp(y(1))*(1-params(10)+params(12)*y(17))*(-exp((-params(5))-y(19)))));
g1(4,1)=(-(T(4)*exp(y(20))*exp((params(5)+y(19))*(-params(4)))*exp(y(1))*getPowerDeriv(exp(y(1)),params(4),1)));
g1(4,11)=exp(y(11));
g1(4,12)=(-(T(3)*exp(y(12))*getPowerDeriv(exp(y(12)),1-params(4),1)));
g1(4,19)=(-(T(4)*T(2)*exp(y(20))*(-params(4))*exp((params(5)+y(19))*(-params(4)))));
g1(4,20)=(-(T(3)*T(4)));
g1(5,12)=(-exp(y(12)));
g1(6,9)=(-exp(y(9)));
g1(6,10)=(-exp(y(10)));
g1(6,11)=exp(y(11));
g1(7,1)=(-(T(5)/(1-params(9))*T(10)*T(11)));
g1(7,10)=(-(T(5)/(1-params(9))*T(6)*T(11)));
g1(7,17)=1;
g1(7,19)=(-(T(5)/(1-params(9))*T(6)*T(11)));
g1(8,1)=(-(T(5)*T(10)*T(12)));
g1(8,10)=(-(T(5)*T(6)*T(12)));
g1(8,18)=1;
g1(8,19)=(-(T(5)*T(6)*T(12)));
g1(9,11)=(-(exp(y(11))*(1-params(4))/exp(y(12))));
g1(9,12)=(-((-(exp(y(12))*exp(y(11))*(1-params(4))))/(exp(y(12))*exp(y(12)))));
g1(9,21)=exp(y(21));
g1(10,10)=exp(y(10));
g1(10,11)=(-exp(y(11)));
g1(10,12)=exp(y(12))*exp(y(21));
g1(10,21)=exp(y(12))*exp(y(21));
g1(10,22)=exp(y(22));
g1(11,1)=(-((1+params(12)*y(6)-params(12))*(params(4)*exp(params(5)+y(19))*(-(exp(y(1))*exp(y(11))))/(exp(y(1))*exp(y(1)))-params(12)*T(10))));
g1(11,10)=(-((1+params(12)*y(6)-params(12))*(-(params(12)*T(6)))));
g1(11,11)=(-(T(7)*(1+params(12)*y(6)-params(12))));
g1(11,17)=(-((1+params(12)*y(6)-params(12))*params(12)/(1+params(12)*y(18)-params(12))));
g1(11,6)=(-(params(12)*T(8)));
g1(11,18)=(-((1+params(12)*y(6)-params(12))*(-(params(12)*(1+params(12)*y(17)-params(10))))/((1+params(12)*y(18)-params(12))*(1+params(12)*y(18)-params(12)))));
g1(11,19)=(-((1+params(12)*y(6)-params(12))*(T(7)-params(12)*T(6))));
g1(11,23)=exp(y(23));
g1(12,16)=T(1)*exp(y(28))*exp((-params(1))*(params(5)+y(27)))*(-y(26))/(y(16)*y(16));
g1(12,26)=T(1)*exp(y(28))*exp((-params(1))*(params(5)+y(27)))*1/y(16);
g1(12,27)=y(26)/y(16)*T(1)*exp(y(28))*(-params(1))*exp((-params(1))*(params(5)+y(27)));
g1(12,28)=T(9);
g1(13,2)=1;
g1(13,9)=(-1);
g1(13,13)=1;
g1(13,19)=(-1);
g1(14,4)=1;
g1(14,11)=(-1);
g1(14,14)=1;
g1(14,19)=(-1);
g1(15,3)=1;
g1(15,10)=(-1);
g1(15,15)=1;
g1(15,19)=(-1);
g1(16,19)=1;
g1(16,30)=(-params(7));
g1(17,7)=(-params(6));
g1(17,20)=1;
g1(17,29)=(-params(8));

end
