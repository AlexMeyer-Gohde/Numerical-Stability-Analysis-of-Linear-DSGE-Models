BETA_ZS_dynare=oo_.dr.ghu(oo_.dr.inv_order_var,:);%BETA_ZS_dynare=BETA_ZS_dynare(1:2,:);
ALPHA_ZS_dynare=[zeros(n,nstatic) oo_.dr.ghx zeros(n,nfwrd)];
ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);%ALPHA_ZS_dynare=ALPHA_ZS_dynare(1:2,1:2);
matrix_quadratic.P=ALPHA_ZS_dynare;
matrix_quadratic.Q=BETA_ZS_dynare;