function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 14);

T = jermann98_rep_def_R_rr_compare.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(10) = exp(params(5)+y(19))*(-(exp(y(1))*exp(y(10))))/(exp(y(1))*exp(y(1)));
T(11) = getPowerDeriv(T(6),1-params(9),1);
T(12) = getPowerDeriv(T(6),(-params(9)),1);
T(13) = getPowerDeriv(exp(y(9))-params(2)*exp(y(2))*exp((-(params(5)+y(19)))),(-params(1)),1);
T(14) = getPowerDeriv(exp(y(25))*exp(params(5)+y(27))-exp(y(9))*params(2),(-params(1)),1);

end
