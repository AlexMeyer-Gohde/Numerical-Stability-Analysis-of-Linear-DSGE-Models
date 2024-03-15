function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 4);

T(1) = params(14)^(1-params(5));
T(2) = exp(y(18))*exp(params(6)*(-params(5)))*exp(y(1))^params(5)*T(1);
T(3) = exp(params(6))*exp(y(10))/exp(y(1));
T(4) = exp(params(6))*params(5)*exp(y(11))/exp(y(1))+(1+params(16)*y(16)-params(13))/(1+params(16)*y(17)-params(16))-params(16)*T(3);

end
