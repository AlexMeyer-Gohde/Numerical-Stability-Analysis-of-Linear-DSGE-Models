function [C,H]=bp(Bhat,Chat,Ahat,D1,D2,Gamma,N,eps3)
% Set Up the System
% Chat * x[t] = Ahat * x[t-1] + Bhat * E(x[t+1]|I[t])
%               + D1 * w[t] + D2 * E(w[t+1]|I[t])
% w[t] = Gamma * w[t-1] + v[t]

% Transform System to Canonical Form:
% x[t] = A * x[t-1] + B * E(x[t+1]|I[t]) 
%        + inv(Chat) * D1 * w[t] + inv(Chat) * D2 * E(w[t+1]|I[t])
B = inv(Chat)*Bhat;
A = inv(Chat)*Ahat;

[dim1,dim2] = size(A);
%N = 100; % Initial Forecasting Horizon (See Binder and Pesaran, 1995)

% Carry Out Recursions
Q = eye(dim1);
R = Chat\(D1*Gamma^N+D2*Gamma^(N+1));
j = 1;
while j <= N
   R = B/Q*R+Chat\(D1*Gamma^(N-j)+D2*Gamma^(N+1-j));
   Q = eye(dim1)-B/Q*A;
   j = j+1;
end

crit3 = 1;
%eps3 = 10^(-6);
QN = Q;
RN = R;

while crit3 > eps3
   Q = QN;
   R = RN;
   QN = eye(dim1);
   RN = Chat\(D1*Gamma^(N+1)+D2*Gamma^(N+2));
   j = 1;
   while j <= (N+1)
      RN = B/QN*RN+Chat\(D1*Gamma^(N+1-j)+D2*Gamma^(N+2-j));
      QN = eye(dim1)-B/QN*A;
      j = j+1;
   end
   crit3 = max(max(abs(QN\RN-Q\R)));
   N = N+1;
end

C = QN\A;
H = QN\R;

% Display Results
%disp('The decision rule is: '),
%disp('x[t] = C*x[t-1] + H*w[t], '),

