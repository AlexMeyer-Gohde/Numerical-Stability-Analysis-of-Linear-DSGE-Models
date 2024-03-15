function [results] = matrix_quadratic_backward_errors(inputs)
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

alpha=norm(A,'fro');
beta=norm(B,'fro');
gamma=norm(C,'fro');

R=A*X^2+B*X+C;
results.lower_backward_error=norm(R,'fro')/((alpha^2*norm(X^2,'fro')^2+beta^2*norm(X,'fro')^2+length(X)*gamma^2)^(1/2));

H=[kron(alpha*(X^2)',eye(size(X))) kron(beta*(X)',eye(size(X))) gamma*eye(2*size(X))];
results.upper_backward_error=norm(-R(:),2)/svds(H,1,'smallest');

P=kron(eye(size(X)),A*X)+kron(X',A)+kron(eye(size(X)),B);
results.conditioning=norm(P\H,2)/norm(X,'fro');
end
