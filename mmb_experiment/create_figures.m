load('First_Run_AMG.mat')
num_results_main=num_results;
load('Fourth_Run_AMG_JS_only_cr_comb.mat')
% num_results_cr=num_results;
% num_results=num_results_cr;
% num_results(:,:,:,99)=zeros([ 15     5     7]);
% num_results(:,:,7,99)=NaN([ 15     5 ]);
num_results(:,:,1:5,:)=num_results_main;
num_results(num_results==0)=NaN;

P_select=[2 4 5 9 7 10 12 14 15];
Q_select_row=[2 4 5 2 4 5 7 10 10 10 12 12 12 14 14 15 15];
Q_select_col=[2 2 2 3 3 3 2 2  3  4  2  3  4   3  4  3  4];
Q_select=sub2ind([15 5],Q_select_row,Q_select_col);
maxs=squeeze(max(num_results(P_select,1,:,:),[],4));
mins=squeeze(min(num_results(P_select,1,:,:),[],4));
results_P=[maxs(1:3,:);mins(4:5,:);maxs(6:end,:)];
results_P_export=compose('%1.3g',results_P);

maxs=squeeze(max(num_results(P_select,5,:,:),[],4));
mins=squeeze(min(num_results(P_select,5,:,:),[],4));
results_PQ=[maxs(1:3,:);mins(4:5,:);maxs(6:end,:)];
results_PQ_export=compose('%1.3g',results_PQ);


maxs=squeeze(max(num_results(:,:,:,:),[],4));
mins=squeeze(min(num_results(:,:,:,:),[],4));
results_Q_min=zeros(length(Q_select_col),size(maxs,3));
results_Q_max=zeros(length(Q_select_col),size(maxs,3));
for j=1:size(maxs,3)
    min_t=mins(:,:,j);
    max_t=maxs(:,:,j);
results_Q_min(:,j)=min_t(Q_select);
results_Q_max(:,j)=max_t(Q_select);
end
results_Q=[results_Q_max(1:6,:);results_Q_min(7,:);results_Q_max(8:end,:)];
results_Q_export=compose('%1.3g',results_Q);

newcolors=[0 0 0;
    0 0 1
1 0 0
0.929000000000000	0.694000000000000	0.125000000000000
0.494000000000000	0.184000000000000	0.556000000000000
0.466000000000000	0.674000000000000	0.188000000000000
0.301000000000000	0.745000000000000	0.933000000000000
0.933000000000000	0.000000000000000	0.933000000000000];


figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(4,1,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Backward Error Bound, P,  Log10')
xlim([-18 -12])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(10,1,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
it=0;
for j=[1:5,7]
    it=it+1;
    obj=squeeze(num_results(12,1,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y,'-.','Color',newcolors(it,:))
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Condition Numbers, P,  Log10')
xlim([-1 15])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
%xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
legend('show');

hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(14,1,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
it=0;
for j=[1:5,7]
    it=it+1;
    obj=squeeze(num_results(15,1,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y,'-.','Color',newcolors(it,:))
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bounds, P  Log10')
xlim([-17 -4])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');

hold off


figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5]
    obj=squeeze(num_results(14,1,j,:))./squeeze(num_results(14,1,7,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1, P,  Log10 rel. to Dynare CR')
xlim([-4 4])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(0,'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off














figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(4,3,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Backward Error Bound, Q 2,  Log10')
xlim([-18 -14])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(10,3,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
it=0;
for j=[1:5,7]
    it=it+1;
    obj=squeeze(num_results(12,3,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y,'-.','Color',newcolors(it,:))
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Condition Numbers, Q 2,  Log10')
xlim([-1 12])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
%xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
legend('show');

hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(14,3,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
it=0;
for j=[1:5,7]
    it=it+1;
    obj=squeeze(num_results(15,3,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y,'-.','Color',newcolors(it,:))
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bounds, Q 2,  Log10')
xlim([-17 -9])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');

hold off


figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:2 4:5 7]
    obj=squeeze(num_results(14,3,j,:))./squeeze(num_results(14,3,3,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
legend('Klein (2000)','Sims (2001)', 'Dynare QZ', 'Anderson (2010)','Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1, Q 2,  Log10 rel. to Uhlig (1999)')
xlim([-2 4])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(0,'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off











figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(4,5,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Backward Error Bound, P Q,  Log10')
xlim([-18 -10])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(10,5,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
it=0;
for j=[1:5,7]
    it=it+1;
    obj=squeeze(num_results(12,5,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y,'-.','Color',newcolors(it,:))
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Condition Numbers, P Q,  Log10')
xlim([-1 15])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
%xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
legend('show');

hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5,7]
    obj=squeeze(num_results(14,5,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
it=0;
for j=[1:5,7]
    it=it+1;
    obj=squeeze(num_results(15,5,j,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y,'-.','Color',newcolors(it,:))
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)', 'Dynare CR','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bounds, P Q,  Log10')
xlim([-17 -4])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');

hold off


figure
set(gcf,'DefaultAxesColorOrder',newcolors)
ax = gca; 
hold on
for j=[1:5]
    obj=squeeze(num_results(14,5,j,:))./squeeze(num_results(14,5,7,:));
    [y,x]=ksdensity(log10(obj(isfinite(obj)))); 
      plot(x,y)
 
end
legend('Klein (2000)','Sims (2001)', 'Uhlig (1999)', 'Dynare QZ', 'Anderson (2010)','AutoUpdate','off')
%legend('QZ','Newton','Baseline','MBI','LS','N-B', 'N-B Col','N-B 1/3','N-B LS', 'N-B Col LS','N-B 1/3 LS', 'N-B Opt', 'N-B Opt LS','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1, P Q,  Log10 rel. to Dynare CR')
xlim([-4 4])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(0,'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off