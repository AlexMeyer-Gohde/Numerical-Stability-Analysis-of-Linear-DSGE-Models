% load('model_info_run')
% model_info(8,:)=sum(model_info(1:4,:),1);endo_var=model_info(end,:);
% model_sizes=[sum(endo_var<20) sum(endo_var<50)-sum(endo_var<20) sum(endo_var<100)-sum(endo_var<50) sum(endo_var<500)-sum(endo_var<100) sum(endo_var<5000)-sum(endo_var<500)]')')

num_results_combined=NaN(15,5,7,99);
load('First_Run_AMG')
num_results_combined(:,:,1:5,:)=num_results;
load('Fourth_Run_AMG_JS_only_cr.mat')
num_results_combined(:,:,7,:)=num_results(:,:,7,:);
num_results=num_results_combined;
num_results(num_results==Inf)=NaN;
num_results(num_results==-Inf)=NaN;
save('Combined_Run_22_11_13.mat')

load('Combined_Run_22_11_13.mat')
number_of_models=99;
n_bins=20;

% impact_results_count_min=impact_results;
% impact_results_count_min(:,8,:)=ones(9,1,99);
% impact_results_count_min(impact_results_count_min=NaN)=1/0;
% 
% impact_results_count_max=impact_results;
% impact_results_count_max(:,8,:)=-ones(9,1,99);
% impact_results_count_min(impact_results_count_min=NaN)=-1/0;

num_results_solab=squeeze(num_results(:,:,1,1:number_of_models));
num_results_gensys=squeeze(num_results(:,:,2,1:number_of_models));
num_results_uhlig=squeeze(num_results(:,:,3,1:number_of_models));
num_results_dynare=squeeze(num_results(:,:,4,1:number_of_models));
num_results_aim=squeeze(num_results(:,:,5,1:number_of_models));
num_results_dynare_cr=squeeze(num_results(:,:,7,1:number_of_models));

% impact_results_solab=squeeze(impact_results(:,1,1:number_of_models));
% impact_results_gensys=squeeze(impact_results(:,2,1:number_of_models));
% impact_results_uhlig=squeeze(impact_results(:,3,1:number_of_models));
% impact_results_dynare=squeeze(impact_results(:,4,1:number_of_models));
% impact_results_aim=squeeze(impact_results(:,5,1:number_of_models));
% impact_results_dynare_cr=squeeze(impact_results(:,7,1:number_of_models));

figure
hold on
temp=log10(num_results_solab(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'c','LineWidth',2)
temp=log10(num_results_gensys(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'k','LineWidth',2)
temp=log10(num_results_uhlig(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'g','LineWidth',2)
temp=log10(num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'b','LineWidth',2)
temp=log10(num_results_aim(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'r','LineWidth',2)
temp=log10(num_results_dynare_cr(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'m','LineWidth',2)
legend('Klein (2000)','Sims (2001)','Uhlig (1999)','Dynare','Anderson (2010)','Dynare CR','AutoUpdate','off')
% [y,x]=ksdensity(log10(num_results_solab(15,:))); plot(x,y,'c','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_gensys(15,:))); plot(x,y,'k','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_uhlig(15,:))); plot(x,y,'g','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_dynare(15,:))); plot(x,y,'b','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_aim(15,:))); plot(x,y,'r','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_dynare_cr(15,:))); plot(x,y,'m','LineWidth',2)
ylabel('density')
xlabel('Forward Error, Log10')
%xlim([-16.5 -12.5])
%ylim([-1 3])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
%yline(0);
xline(log10(eps),'-.','Machine Precision','LineWidth',2);
hold off

figure
hold on
temp=log10(num_results_solab(14,:)./num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'c','LineWidth',2)
temp=log10(num_results_gensys(14,:)./num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'k','LineWidth',2)
temp=log10(num_results_uhlig(14,:)./num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'g','LineWidth',2)
%temp=log10(num_results_dynare(14,:)./num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'b','LineWidth',2)
temp=log10(num_results_aim(14,:)./num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'r','LineWidth',2)
temp=log10(num_results_dynare_cr(14,:)./num_results_dynare(14,:));temp=temp(isfinite(temp)==1);[y,x]=ksdensity(temp); plot(x,y,'m','LineWidth',2)
legend('Klein (2000)','Sims (2001)','Uhlig (1999)','Anderson (2010)','Dynare CR','AutoUpdate','off')
% [y,x]=ksdensity(log10(num_results_solab(15,:))); plot(x,y,'c','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_gensys(15,:))); plot(x,y,'k','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_uhlig(15,:))); plot(x,y,'g','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_dynare(15,:))); plot(x,y,'b','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_aim(15,:))); plot(x,y,'r','LineWidth',2)
% [y,x]=ksdensity(log10(num_results_dynare_cr(15,:))); plot(x,y,'m','LineWidth',2)
ylabel('density')
xlabel('Forward Error, Log10 relative to Dynare')
xlim([-2 2])
%ylim([-1 3])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
%yline(0);
xline(0,'-.','LineWidth',2);
hold off