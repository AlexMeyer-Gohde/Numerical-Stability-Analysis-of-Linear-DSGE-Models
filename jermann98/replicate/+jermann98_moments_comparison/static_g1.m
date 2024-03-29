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
    T = jermann98_moments_comparison.static_g1_tt(T, y, x, params);
end
g1 = zeros(13, 13);
g1(1,2)=(-((exp(y(2))-exp(y(2))*params(2)*exp((-params(6))))*getPowerDeriv(exp(y(2))-exp(y(2))*params(2)*exp((-params(6))),(-params(1)),1)-params(2)*params(3)*(exp(y(2))*exp(params(6))-exp(y(2))*params(2))*getPowerDeriv(exp(y(2))*exp(params(6))-exp(y(2))*params(2),(-params(1)),1)));
g1(1,8)=1;
g1(2,13)=exp(y(13));
g1(3,1)=exp(y(1))-exp(y(1))*exp((-params(6)))*(1-params(13)+params(16)*y(9));
g1(3,3)=(-((1-params(16))*exp(y(3))));
g1(3,9)=(-(exp(y(1))*exp((-params(6)))*params(16)));
g1(4,1)=(-(T(1)*exp(y(11))*exp(params(6)*(-params(5)))*exp(y(1))*getPowerDeriv(exp(y(1)),params(5),1)));
g1(4,4)=exp(y(4));
g1(4,11)=(-T(2));
g1(5,2)=(-exp(y(2)));
g1(5,3)=(-exp(y(3)));
g1(5,4)=exp(y(4));
g1(6,1)=(-(params(11)/(1-params(10))*T(5)*T(6)));
g1(6,3)=(-(params(11)/(1-params(10))*T(3)*T(6)));
g1(6,9)=1;
g1(7,1)=(-(params(11)*T(5)*T(7)));
g1(7,3)=(-(params(11)*T(3)*T(7)));
g1(7,10)=1;
g1(8,1)=(-((1+params(16)*y(10)-params(16))*(exp(params(6))*params(5)*(-(exp(y(1))*exp(y(4))))/(exp(y(1))*exp(y(1)))-params(16)*T(5))));
g1(8,3)=(-((1+params(16)*y(10)-params(16))*(-(params(16)*T(3)))));
g1(8,4)=(-(exp(params(6))*params(5)*exp(y(4))/exp(y(1))*(1+params(16)*y(10)-params(16))));
g1(8,9)=(-params(16));
g1(8,10)=(-(params(16)*T(4)+(1+params(16)*y(10)-params(16))*(-(params(16)*(1+params(16)*y(9)-params(13))))/((1+params(16)*y(10)-params(16))*(1+params(16)*y(10)-params(16)))));
g1(8,12)=exp(y(12));
g1(9,12)=exp(params(6)*(-params(1)))*params(3)*exp(y(12));
g1(10,5)=1;
g1(11,6)=1;
g1(12,7)=1;
g1(13,11)=1-params(7);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
