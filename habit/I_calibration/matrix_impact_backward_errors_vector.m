function [results] = matrix_impact_backward_errors_vector(inputs)
%Matrix quadratic backward errors and conditioning number follow Higham and Kim (2001)
%0=A*X^2+B*X+C
%X is the n by n solvent
%A, B, and C are n by n matrices
%Alexander Meyer-Gohde
%31/03/2021
A=inputs.A;
B=inputs.B;
C=inputs.C;
X=inputs.X;
D=inputs.D;
W=inputs.W;
F=A*X+B;
results=NaN(3,3);
alpha=norm(A,'fro');
beta=norm(B,'fro');
gamma=norm(C,'fro');
delta=norm(D,'fro');
xi=norm(X,'fro');
phi=norm(F,'fro');
R=F*W+D;
n=size(X,1);
ne=size(W,2);
V=kron(eye(n),A*X+B)+kron(X',A);
QAV=kron(W',A)/V;
QT=kron(W',eye(n))-QAV;


%results.lower_backward_error=norm(R,'fro')/((alpha^2*norm(X^2,'fro')^2+beta^2*norm(X,'fro')^2+length(X)*gamma^2)^(1/2));
results(1,1)=norm(R,'fro')/((phi^2*norm(W,'fro')^2+delta^2*ne)^(1/2));%norm(R,'fro')/((phi^2*svds(W,1,'largest')^2+delta^2)^(1/2)); %
results(1,2)=norm(R,'fro')/((xi^2*svds(W,1,'largest')^2*svds(A,1,'largest')^2+alpha^2*svds(X*W,1,'largest')^2 +beta^2*svds(W,1,'largest')^2+delta^2)^(1/2));
results(1,3)=norm(R,'fro')/((alpha^2*svds(QT*kron((X')^2,eye(n)),1,'largest')^2 +beta^2*svds(QT*kron(X',eye(n)),1,'largest')^2+gamma^2*svds(QAV,1,'largest')^2+delta^2)^(1/2));

if length(A)<30
try
results(2,1)=norm(R,'fro')/((phi^2*svds(W,1,'smallest')^2+delta^2)^(1/2)); 
results(2,2)=norm(R,'fro')/((xi^2*svds(W,1,'smallest')^2*svds(A,1,'smallest')^2+alpha^2*svds(X*W,1,'smallest')^2 +beta^2*svds(W,1,'smallest')^2+delta^2)^(1/2));
results(2,3)=norm(R,'fro')/((alpha^2*svds(QT*kron(X',eye(n)),1,'smallest')^2 +beta^2*svds(QT,1,'smallest')^2+gamma^2*svds(QAV,1,'smallest')^2+delta^2)^(1/2));
  

%H=[kron(alpha*(X^2)',eye(size(X))) kron(beta*(X)',eye(size(X))) gamma*eye(size(X).^2)];
%results.upper_backward_error=norm(-R(:),2)/svds(H,1,'smallest');
%results(2)=norm(-R(:),2)/svds(H,1,'smallest');
%results(2)=norm(R(:),2)/((alpha^2*svds(X^2,1,'smallest')^2+beta^2*svds(X,1,'smallest')^2+gamma^2)^(1/2));
catch
    results(2,1:3)=NaN;
end

try
%P=kron(eye(size(X)),A*X)+kron(X',A)+kron(eye(size(X)),B);
%results.conditioning=norm(P\H,2)/norm(X,'fro');
results(3,1)=norm(kron(eye(ne),F)\([phi*kron(W',eye(n)) delta*eye(n*ne)]),2)/norm(W,'fro');
results(3,2)=norm(kron(eye(ne),F)\([xi*kron(W',A) alpha*kron(W'*X',eye(n)) beta*kron(W',eye(n)) delta*eye(n*ne)]),2)/norm(W,'fro');
results(3,3)=norm(kron(eye(ne),F)\([alpha*QT*kron(X',eye(n)) beta*QT -gamma*QAV delta*eye(n*ne)]),2)/norm(W,'fro');
catch
    results(3,1:3)=NaN;
end
end
end
