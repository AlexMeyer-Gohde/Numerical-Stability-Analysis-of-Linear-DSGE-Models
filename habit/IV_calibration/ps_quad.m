function y = ps(X, m, tol, rl, marksize)
%PS     Dot plot of a pseudospectrum.
%       PS(A, M, TOL, RL) plots an approximation to a pseudospectrum
%       of the square matrix A, using M random perturbations of size TOL.
%       M defaults to a SIZE(A)-dependent value and TOL to 1E-3.
%       RL defines the type of perturbation:
%         RL =  0 (default): absolute complex perturbations of 2-norm TOL.
%         RL =  1:           absolute real perturbations of 2-norm TOL.
%         RL = -1:           componentwise real perturbations of size TOL.
%       The eigenvalues of A are plotted as crosses `x'.
%       PS(A, M, TOL, RL, MARKSIZE) uses the specified marker size instead
%       of a size that depends on the figure size, the matrix order, and M.
%       If MARKSIZE < 0, the plot is suppressed and the plot data is returned
%       as an output argument.
%       PS(A, 0) plots just the eigenvalues of A.

%       For a given TOL, the pseudospectrum of A is the set of
%       pseudo-eigenvalues of A, that is, the set
%       { e : e is an eigenvalue of A+E, for some E with NORM(E) <= TOL }.
%
%       References:
%       L. N. Trefethen, Computation of pseudospectra, Acta Numerica,
%          8:247-295, 1999.
%       L. N. Trefethen, Spectra and pseudospectra, in The Graduate
%          Student's Guide to Numerical Analysis '98, M. Ainsworth,
%          J. Levesley, and M. Marletta, eds., Springer-Verlag, Berlin,
%          1999, pp. 217-250.

%if diff(size(A)), error('Matrix must be square.'), end
A=X.A;
B=X.B;
C=X.C;

n = length(A);

if nargin < 5, marksize = 0; end
if nargin < 4, rl = 0; end
if nargin < 3, tol = 1e-3; end
if nargin < 2 | isempty(m), m = 5*max(1, round( 25*exp(-0.047*n) )); end

if m == 0
   e = polyeig(C,B,A);   ax = cpltaxes(e);
   plot(real(e), imag(e), 'x')
   axis(ax), axis('square')
   return
end

x = zeros(m*2*n,1);
i = sqrt(-1);

AA=[ zeros(size(A)) -A;  eye(size(A)) zeros(size(A))];
BB=[C B;   zeros(size(A)) eye(size(A))];
z = zeros(m*2*n,1);

for j = 1:m
   if rl == -1     % Componentwise.
      dA = -ones(n) + 2*rand(n);   % Uniform random numbers on [-1,1].
      dA = tol * A .* dA;
            dB = -ones(n) + 2*rand(n);   % Uniform random numbers on [-1,1].
      dB = tol * B .* dB;
            dC = -ones(n) + 2*rand(n);   % Uniform random numbers on [-1,1].
      dC = tol * C .* dC;
            dAA = -ones(2*n) + 2*rand(2*n);   % Uniform random numbers on [-1,1].
      dAA = tol * AA .* dAA;
            dBB = -ones(2*n) + 2*rand(2*n);   % Uniform random numbers on [-1,1].
      dBB = tol * BB .* dBB;
   else
      if rl == 0   % Complex absolute.
         dA = randn(n) + i*randn(n);
         dB = randn(n) + i*randn(n);
         dC = randn(n) + i*randn(n);         
         dAA = randn(2*n) + i*randn(2*n);
         dBB = randn(2*n) + i*randn(2*n);
      else         % Real absolute.
         dA = randn(n);
         dB = randn(n);
         dC = randn(n);
         dAA = randn(2*n);
         dBB = randn(2*n);
      end
      dA = tol/norm(dA)*dA;
      dB = tol/norm(dB)*dB;
      dC = tol/norm(dC)*dC;
      dAA = tol/norm(dAA)*dAA;
      dBB = tol/norm(dBB)*dBB;
   end
   x((j-1)*2*n+1:j*2*n) =  polyeig(C+dC,B+dB,A+dA);
   [s,t,~,~] = qz(AA+dAA,BB+dBB);
   z((j-1)*2*n+1:j*2*n) =  diag(s).\diag(t);
end

x=x(abs(x)<10);
z=z(abs(z)<10);

if marksize >= 0

   ax = cpltaxes(x);
   h = ezplot(@(x,y)x.^2+y.^2-1,[-1, 1, -1,1]);
   hold on
   set(h,'Color','k', 'LineWidth', 2); xlabel('');ylabel('');title('');
   h = plot(real(z),imag(z),'r.');
   axis(ax), axis('square')
      
   h = plot(real(x),imag(x),'b.');
   hold off

   % Next block adapted from SPY.M.
   if marksize == 0
      %units = get(gca,'units');
      %set(gca,'units','points');
      pos = get(gca,'position');
      nps = 2.4*sqrt(n*m);  % Factor based on number of pseudo-ei'vals plotted.
      myguess = 5*round(3*min(pos(3:4))/nps);
      marksize = max(1,myguess);
      %set(gca,'units',units);
   end

   hold on
   e = polyeig(C,B,A);
   plot(real(e),imag(e),'ko','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
   %set(h,'markersize',2*marksize);
   hold off

else

  y = x;

end

function x = cpltaxes(z)
%CPLTAXES   Determine suitable AXIS for plot of complex vector.
%           X = CPLTAXES(Z), where Z is a complex vector,
%           determines a 4-vector X such that AXIS(X) sets axes for a plot
%           of Z that has axes of equal length and leaves a reasonable amount
%           of space around the edge of the plot.

%           Called by FV, GERSH, PS and PSCONT.

% Set x and y axis ranges so both have the same length.

xmin = min(real(z)); xmax = max(real(z));
ymin = min(imag(z)); ymax = max(imag(z));

% Fix for rare case of `trivial data'.
if xmin == xmax, xmin = xmin - 1/2; xmax = xmax + 1/2; end
if ymin == ymax, ymin = ymin - 1/2; ymax = ymax + 1/2; end

if xmax-xmin >= ymax-ymin
   ymid = (ymin + ymax)/2;
   ymin =  ymid - (xmax-xmin)/2; ymax = ymid + (xmax-xmin)/2;
else
   xmid = (xmin + xmax)/2;
   xmin = xmid - (ymax-ymin)/2; xmax = xmid + (ymax-ymin)/2;
end
axis('square')

% Scale ranges by 1+2*alpha to give extra space around edges of plot.

alpha = 0.1;
x(1) = xmin - alpha*(xmax-xmin);
x(2) = xmax + alpha*(xmax-xmin);
x(3) = ymin - alpha*(ymax-ymin);
x(4) = ymax + alpha*(ymax-ymin);

if x(1) == x(2), x(2) = x(2) + 0.1; end
if x(3) == x(4), x(4) = x(3) + 0.1; end
