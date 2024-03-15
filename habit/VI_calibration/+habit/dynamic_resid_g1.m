function [residual, g1] = dynamic_resid_g1(T, y, x, params, steady_state, it_, T_flag)
% function [residual, g1] = dynamic_resid_g1(T, y, x, params, steady_state, it_, T_flag)
%
% Wrapper function automatically created by Dynare
%

    if T_flag
        T = habit.dynamic_g1_tt(T, y, x, params, steady_state, it_);
    end
    residual = habit.dynamic_resid(T, y, x, params, steady_state, it_, false);
    g1       = habit.dynamic_g1(T, y, x, params, steady_state, it_, false);

end
