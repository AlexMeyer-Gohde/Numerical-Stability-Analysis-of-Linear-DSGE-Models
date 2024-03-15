VARNAMES = ['consumption',
            'capital    ',
            'technology '];
        
       h=Y(1);
beta=Y(2);
delta=Y(3);
alpha=Y(4);
sigma=Y(5);
rho=Y(6);
omega=Y(7);
K_aminus1=(1/beta-1+delta)/alpha;
CK=K_aminus1-delta;
A=[-sigma/(1-h) 0; 0 0];
B=[sigma*(1+h)/(1-h) alpha*beta*K_aminus1*(alpha-1);-CK -1];
C=[-sigma*h/(1-h) 0;0 (1-delta+alpha*K_aminus1)];
D=[beta*alpha*K_aminus1*rho; K_aminus1];

AA = [];
BB = [];
CC=[];
DD=[];
FF=A;
GG =B;
HH=C;
JJ=[];
KK=[];
LL=zeros(2,1);
MM = D;
NN = [rho];
Sigma = [ omega  ];

% Setting the options:
QQ=zeros(2,1);
  
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 3; % Index of output among the variables selected for HP filter
IMP_SELECT = [1:5,7];


% Starting the calculations:
DO_IMPRESP=0;
DO_SIMUL=0;
DO_MOMENTS=0;

do_it;

rp_uhlig=4*100*sigma/(1-h)*QQ(1)*alpha*K_aminus1*beta*(omega)^2
AA=[PP(1,1)-1  PP(1,2) QQ(1)];
BB=[PP QQ;0 0 rho ];
DD=[0;0;1];
log_c_std_uhlig=(kron(AA,AA)*((eye(9,9)-kron(BB,BB))\kron(DD,DD))*omega^2)^(1/2)