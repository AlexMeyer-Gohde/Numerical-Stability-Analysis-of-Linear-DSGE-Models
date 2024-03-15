function [residual, g1] = dynamic_resid_g1(T, y, x, params, steady_state, it_, T_flag)
% function [residual, g1] = dynamic_resid_g1(T, y, x, params, steady_state, it_, T_flag)
%
% Wrapper function automatically created by Dynare
%

    if T_flag
        T = jermann98_moments_comparison.dynamic_g1_tt(T, y, x, params, steady_state, it_);
    end
    residual = jermann98_moments_comparison.dynamic_resid(T, y, x, params, steady_state, it_, false);
    g1       = jermann98_moments_comparison.dynamic_g1(T, y, x, params, steady_state, it_, false);

end