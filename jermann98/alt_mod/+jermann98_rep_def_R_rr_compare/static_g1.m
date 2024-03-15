function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = jermann98_rep_def_R_rr_compare.static_g1_tt(T, y, x, params);
end
g1 = zeros(17, 17);
g1(1,2)=(-((exp(y(2))-exp(y(2))*params(2)*exp((-(params(5)+y(12)))))*T(12)-params(2)*T(1)*(exp(y(2))*exp(params(5)+y(12))-exp(y(2))*params(2))*T(13)));
g1(1,9)=1;
g1(1,12)=(-(T(12)*(-(exp(y(2))*params(2)*(-exp((-(params(5)+y(12)))))))-params(2)*T(1)*exp(y(2))*exp(params(5)+y(12))*T(13)));
g1(2,12)=(-(T(1)*(-params(1))*exp((params(5)+y(12))*(-params(1)))));
g1(2,17)=exp(y(17));
g1(3,1)=exp(y(1))-exp(y(1))*(1-params(10)+params(12)*y(10))*exp((-params(5))-y(12));
g1(3,3)=(-((1-params(12))*exp(y(3))));
g1(3,10)=(-(exp(y(1))*params(12)*exp((-params(5))-y(12))));
g1(3,12)=(-(exp(y(1))*(1-params(10)+params(12)*y(10))*(-exp((-params(5))-y(12)))));
g1(4,1)=(-(T(4)*exp(y(13))*exp((params(5)+y(12))*(-params(4)))*exp(y(1))*getPowerDeriv(exp(y(1)),params(4),1)));
g1(4,4)=exp(y(4));
g1(4,5)=(-(T(3)*exp(y(5))*getPowerDeriv(exp(y(5)),1-params(4),1)));
g1(4,12)=(-(T(4)*T(2)*exp(y(13))*(-params(4))*exp((params(5)+y(12))*(-params(4)))));
g1(4,13)=(-(T(3)*T(4)));
g1(5,5)=(-exp(y(5)));
g1(6,2)=(-exp(y(2)));
g1(6,3)=(-exp(y(3)));
g1(6,4)=exp(y(4));
g1(7,1)=(-(T(5)/(1-params(9))*T(9)*T(10)));
g1(7,3)=(-(T(5)/(1-params(9))*T(6)*T(10)));
g1(7,10)=1;
g1(7,12)=(-(T(5)/(1-params(9))*T(6)*T(10)));
g1(8,1)=(-(T(5)*T(9)*T(11)));
g1(8,3)=(-(T(5)*T(6)*T(11)));
g1(8,11)=1;
g1(8,12)=(-(T(5)*T(6)*T(11)));
g1(9,4)=(-(exp(y(4))*(1-params(4))/exp(y(5))));
g1(9,5)=(-((-(exp(y(5))*exp(y(4))*(1-params(4))))/(exp(y(5))*exp(y(5)))));
g1(9,14)=exp(y(14));
g1(10,3)=exp(y(3));
g1(10,4)=(-exp(y(4)));
g1(10,5)=exp(y(5))*exp(y(14));
g1(10,14)=exp(y(5))*exp(y(14));
g1(10,15)=exp(y(15));
g1(11,1)=(-((1+params(12)*y(11)-params(12))*(exp(params(5)+y(12))*params(4)*(-(exp(y(1))*exp(y(4))))/(exp(y(1))*exp(y(1)))-params(12)*T(9))));
g1(11,3)=(-((1+params(12)*y(11)-params(12))*(-(params(12)*T(6)))));
g1(11,4)=(-(T(7)*(1+params(12)*y(11)-params(12))));
g1(11,10)=(-params(12));
g1(11,11)=(-(params(12)*T(8)+(1+params(12)*y(11)-params(12))*(-(params(12)*(1+params(12)*y(10)-params(10))))/((1+params(12)*y(11)-params(12))*(1+params(12)*y(11)-params(12)))));
g1(11,12)=(-((1+params(12)*y(11)-params(12))*(T(7)-params(12)*T(6))));
g1(11,16)=exp(y(16));
g1(12,12)=T(1)*exp(y(16))*(-params(1))*exp((params(5)+y(12))*(-params(1)));
g1(12,16)=exp((params(5)+y(12))*(-params(1)))*T(1)*exp(y(16));
g1(13,6)=1;
g1(13,12)=(-1);
g1(14,7)=1;
g1(14,12)=(-1);
g1(15,8)=1;
g1(15,12)=(-1);
g1(16,12)=1;
g1(17,13)=1-params(6);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
