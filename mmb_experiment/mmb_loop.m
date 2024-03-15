%Description: ....
%....
%Alexander Meyer-Gohde, Johanna Saecker

clear
%load First_Run_AMG
%loop_k_start=loop_k;
loop_k_start=1;

addpath C:\dynare\5.1\matlab
addpath('..\algorithm\')
addpath('..\mmb_replication\')
addpath('..\solution_methods\')
%YourPath=''
%YourPath='C:\Users\saecker.ITS\PowerFolders\Newton_Saecker (Alexander Meyer-Gohde)\code\workspace\2021_10_15_js';
YourPath=pwd;
cd (YourPath)
addpath(pwd)

%fileID = fopen('mmb_test.txt','r');
fileID = fopen('mmb_names.txt','r');
mmbline = fgetl(fileID);        
%geht hier vielleicht auch besser (ff. 5 Zeilen von mathworks-website)
mmb_vec = cell(0,1);            
while ischar(mmbline)           
    mmb_vec{end+1,1} = mmbline; 
    mmbline = fgetl(fileID);    
end    
loop_n=size(mmb_vec,1);

%AMG_JS_Results=zeros(5,7,loop_n);

for loop_k=loop_k_start:loop_n
    %k=1;    %for testing
    %k=1;    %for testing
    [loop_k loop_k/loop_n]
%change directory to folder path in MMB
cd([YourPath '..\..\mmb_replication\mmb-rep-master_names\' mmb_vec{loop_k} '\' mmb_vec{loop_k} '_rep'])

%run dynare
%dynare ([mmb_vec{k} '_rep']) 
eval(['dynare ', mmb_vec{loop_k}, '_rep noclearall nograph nostrict'])
%%%%% current problem: dynare_to_matrix_quadratic needs to be located in
%%%%% mmb-rep-folders (e.g. BRA_SAMBA08_rep)
%AMG_JS_Results(1,:,loop_k)=[M_.nstatic, M_.nfwrd, M_.npred, M_.nboth, M_.nsfwrd, M_.nspred, M_.ndynamic];

create_matrix_quadratic_from_dynare;

BETA_ZS_dynare=oo_.dr.ghu(oo_.dr.inv_order_var,:);%BETA_ZS_dynare=BETA_ZS_dynare(1:2,:);
ALPHA_ZS_dynare=[zeros(n,nstatic) oo_.dr.ghx zeros(n,nfwrd)];
ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);%ALPHA_ZS_dynare=ALPHA_ZS_dynare(1:2,1:2);
matrix_quadratic.P=ALPHA_ZS_dynare;
matrix_quadratic.Q=BETA_ZS_dynare;
if n>100; [errors]=dsge_backward_errors_condition_sparse_minimal(matrix_quadratic); else [errors]=dsge_backward_errors_condition_full(matrix_quadratic);end
% eigval=oo_.dr.eigval;
% l_1=eigval(abs(eigval)>1);
% l_2=eigval(abs(eigval)<1);
% chordal_dist=abs(l_1-l_2')./((abs(l_1).^2+1).^(1/2).*(abs(l_2').^2+1).^(1/2));
% eig_seperation_dynare=min(min(chordal_dist));
num_results(:,:,4,loop_k)=errors;

[ALPHA_ZS_aim,BETA_ZS_aim,sorted_eigenvalues,existence_uniqueness]=...
    sp_solve(matrix_quadratic.A,matrix_quadratic.B,matrix_quadratic.C,zeros(M_.endo_nbr,M_.exo_nbr),matrix_quadratic.D,zeros(M_.exo_nbr,M_.exo_nbr),M_.endo_nbr,M_.endo_nbr,M_.exo_nbr,1);
matrix_quadratic.P=ALPHA_ZS_aim;
matrix_quadratic.Q=BETA_ZS_aim;
if n>100; [errors]=dsge_backward_errors_condition_sparse_minimal(matrix_quadratic); else [errors]=dsge_backward_errors_condition_full(matrix_quadratic);end
num_results(:,:,5,loop_k)=errors;

try
[BETA_ZS_solab,ll,ff,ALPHA_ZS_solab,eigval_solab]  = solab([ zeros(M_.endo_nbr,M_.endo_nbr) -matrix_quadratic.A;  eye(M_.endo_nbr,M_.endo_nbr) zeros(M_.endo_nbr,M_.endo_nbr)],[matrix_quadratic.C matrix_quadratic.B;   zeros(M_.endo_nbr,M_.endo_nbr) eye(M_.endo_nbr,M_.endo_nbr)],zeros(M_.exo_nbr,M_.exo_nbr),[matrix_quadratic.D; zeros(M_.endo_nbr,M_.exo_nbr)],M_.endo_nbr);
matrix_quadratic.P=ALPHA_ZS_solab;
matrix_quadratic.Q=BETA_ZS_solab;
if n>100; [errors]=dsge_backward_errors_condition_sparse_minimal(matrix_quadratic); else [errors]=dsge_backward_errors_condition_full(matrix_quadratic);end
num_results(:,:,1,loop_k)=errors;
end

try
addpath([YourPath '\..\solution_methods\gensys'])
g0=[eye(size(matrix_quadratic.A)), matrix_quadratic.B;zeros(size(matrix_quadratic.A)) matrix_quadratic.A];
g1=[zeros(size(matrix_quadratic.A)), -matrix_quadratic.C; eye(size(matrix_quadratic.A)) zeros(size(matrix_quadratic.A))];
c=zeros(2*size(matrix_quadratic.A,1),1);
psi=[-matrix_quadratic.D; zeros(size(matrix_quadratic.D))];
%pi=[matrix_quadratic.A; zeros(size(matrix_quadratic.A))];
pi=[zeros(size(matrix_quadratic.A)); matrix_quadratic.A];


% g0=[matrix_quadratic.A, matrix_quadratic.B;zeros(size(matrix_quadratic.A)) eye(size(matrix_quadratic.A))];
% %g1=[zeros(size(matrix_quadratic.A)), -matrix_quadratic.C; eye(size(matrix_quadratic.A)) zeros(size(matrix_quadratic.A))];
% g1=[zeros(size(matrix_quadratic.A)), -matrix_quadratic.C; eye(size(matrix_quadratic.A)) zeros(size(matrix_quadratic.A))];
% c=zeros(2*size(matrix_quadratic.A,1),1);
% psi=[-matrix_quadratic.D; zeros(size(matrix_quadratic.D))];
% %pi=[matrix_quadratic.A; zeros(size(matrix_quadratic.A))];
% pi=[zeros(size(matrix_quadratic.A)); eye(size(matrix_quadratic.A))];

[G1,Cgensys,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi);
ALPHA_ZS_gensys=G1(size(matrix_quadratic.A,1)+1:end,size(matrix_quadratic.A,1)+1:end);
BETA_ZS_gensys=impact(size(matrix_quadratic.A,1)+1:end,:);
matrix_quadratic.P=ALPHA_ZS_gensys;
matrix_quadratic.Q=BETA_ZS_gensys;
if n>100; [errors]=dsge_backward_errors_condition_sparse_minimal(matrix_quadratic); else [errors]=dsge_backward_errors_condition_full(matrix_quadratic);end
num_results(:,:,2,loop_k)=errors;
rmpath([YourPath '\..\solution_methods\gensys'])
end

if n<100
try
addpath([YourPath '\..\solution_methods\uhlig'])

AA=[]; BB=[]; CC=[]; DD=[];
FF=matrix_quadratic.A;
GG=matrix_quadratic.B;
HH=matrix_quadratic.C;
JJ=[]; KK=[];
MM=matrix_quadratic.D;
LL=zeros(size(MM));
NN=zeros(size(matrix_quadratic.D,2),size(matrix_quadratic.D,2));

[l_equ,m_states] = size(FF);
[l_equ,n_endog ] = size(FF);
[l_equ,k_exog  ] = size(MM);
warnings=[];
options;
solve;
ALPHA_ZS_uhlig=PP;
BETA_ZS_uhlig=QQ;
matrix_quadratic.P=ALPHA_ZS_uhlig;
matrix_quadratic.Q=BETA_ZS_uhlig;
if n>100; [errors]=dsge_backward_errors_condition_sparse_minimal(matrix_quadratic); else [errors]=dsge_backward_errors_condition_full(matrix_quadratic);end
num_results(:,:,3,loop_k)=errors;
rmpath([YourPath '\..\solution_methods\uhlig'])
end
end

cd([YourPath])

%run different Newton methods
%newton_solvent_2(matrix_quadratic);
%newton_solvent_3(matrix_quadratic);


sprintf('Iteration %d of %d',loop_k,loop_n)
%save certain results somewhere
clearvars -except loop_k loop_n num_results impact_results YourPath mmb_vec
save First_Run_AMG
end

