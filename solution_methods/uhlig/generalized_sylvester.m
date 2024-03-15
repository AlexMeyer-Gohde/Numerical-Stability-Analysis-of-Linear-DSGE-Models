function [Y]=generalized_sylvester(A,B,C,D) 
%solves A * Y + B * Y * C + D=0
n=length(A);
Y=zeros(n,n);                   % create Y-matrix
TY=zeros(n,n);

[H,T,Ws,Z]=hess(A,B);       % equation (3.3) with Hessenberg-triangular H instead of S
[U,R]=schur(C);           % equation (3.3)
F=-Ws*D*U;                  % equation (3.4)
for k=1:n                   % loop through columns of Y
    lhs=R(k,k)*T;
    lhs=lhs+H;
    rhs=TY(:,1:k-1)*R(1:k-1,k);
    rhs=F(:,k)-rhs;
    Y(:,k)=lhs\rhs;
    TY(:,k)=T*Y(:,k);
end
%  max(max(abs(Y_amg-Y)))
Y=Z*Y*U';
