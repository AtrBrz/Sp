% Script Spline models of Regression analysis
dati=xlsread('DailyDelhiClimateTest');
%Time
x=dati(:,1);
%Mean_Temperature
y=dati(:,2);
%Dimension of data
n=length(x);
figure 
%Scatterplot
plot (x,y,'ko')
xlabel('Time (Days)')
ylabel('Mean Temperature')
title('Daily Climate')

%Different methods of smoothing spline with lambda=1
tiledlayout(3,1)
nexttile 
figure
plot(x,y,'ko')
hold on
%Lagrangian parameter
S=1e-5;
%Weights_i=1 i=1,...,114
w=ones(n,1);
%H.R Algorithm
[a,b,c,d,p,iter] = my_reinsch(x,y,S,w);
for i = 1:n-1
xx = linspace(x(i),x(i+1),n);
yy = a(i)*(xx - x(i)).^(3) + b(i)*(xx - x(i)).^(2) + c(i)*(xx - x(i)) + d(i);
t1=plot(xx,yy,"b");
end
title(['Smoothing spline R-H con l = ',num2str(p)])
hold off

%Spaps Algorithm
nexttile 
plot(x,y,'ko')
hold on
[sp,values,rho_spaps] = spaps(x,y,S,w);
for i = 1:n-1
xx = linspace(x(i),x(i+1),n);
yy_spaps = spval(sp,xx);
t2=plot(xx,yy_spaps,"r");
end
title(['Smoothing spline Spaps con \rho = ',num2str(rho_spaps)])
hold off

%Csaps Algorithm
nexttile 
plot(x,y,'ko')
hold on
smooth_factor_csaps = rho_spaps/(1+rho_spaps);
pp = csaps(x,y,smooth_factor_csaps,[],w);
for i = 1:n-1
xx = linspace(x(i),x(i+1),n);
yy_csaps=ppval(pp,xx);
t3=plot(xx,yy_csaps,"g");
end
title(['Smoothing spline Cpaps con p = ',num2str(smooth_factor_csaps)])
hold off

%Final plot with all smoothin methods
figure
plot(x,y,'ko')
hold on
%Smoothing spline
S=1e-5;
w=ones(n,1);
[a,b,c,d,p,iter] = my_reinsch(x,y,S,w);
for i = 1:n-1
xx = linspace(x(i),x(i+1),n);
yy = a(i)*(xx - x(i)).^(3) + b(i)*(xx - x(i)).^(2) + c(i)*(xx - x(i)) + d(i);
plot(xx,yy,"b");
end

[sp,values,rho_spaps] = spaps(x,y,S,w);
for i = 1:n-1
xx = linspace(x(i),x(i+1),n);
yy_spaps = spval(sp,xx);
plot(xx,yy_spaps,"r");
end

smooth_factor_csaps = rho_spaps/(1+rho_spaps);
pp = csaps(x,y,smooth_factor_csaps,[],w);
for i = 1:n-1
xx = linspace(x(i),x(i+1),n);
yy_csaps=ppval(pp,xx);
plot(xx,yy_csaps,"g");
end
hold off
title(['Smoothing spline with different methods'])
legend('Data','R-H Algorithm','Spaps Algorithm','Csaps Algorithm', ...
    'Location','northwest')


%Regression P-spline
%Different number of K<<n
figure
number_breakpoints=[10,20,30,40,50,60];
t=subplot(2,3,length(number_breakpoints));
for i=1:length(number_breakpoints)
[z,lambda]= psplineBspline_mod(x,y,number_breakpoints(i));
modelfitsep(i,:)=z'; %TABLE TO COMPARE RESULTS
subplot(2,3,i)
plot(x,y,'ko')
hold on
plot(x,z,'-.')
title(['Lambda=',num2str(lambda),'  Br-Points=',num2str(number_breakpoints(i))])
axis tight
grid on
end

%Compare results
%Mean value and standard deviation to compute confidence interval at 95%
mn = mean(modelfitsep,1); 
se = std(modelfitsep);
ff=figure(300)
for i=1:length(number_breakpoints)
numnodes=number_breakpoints(i);
hold on
yvec=modelfitsep(i,:)';
plot(xx,yvec)
end
hold on
y2vec = mn'-1.96*(se')./sqrt(length(number_breakpoints));
y1vec = mn';
y3vec = mn'+1.96*(se')./sqrt(length(number_breakpoints));
plot_ci(xx,[y1vec y2vec y3vec], 'PatchColor', 'r', 'PatchAlpha', 0.2, ...
'MainLineWidth', 1.2, 'MainLineStyle', '-.', 'MainLineColor', 'b', ...
'LineWidth', 1.2, 'LineStyle','--', 'LineColor', 'k');
grid on
title(['Computed values and 95% confidence band with lambda=',num2str(lambda)])
axis tight
grid on
legend('K=10','K=20','K=30','K=40','K=50','K=60','Location','northwest')
