% make plots for the simulation experiment for the paper 
% "On robust recovery of signals from indirect observations" 
% by Y. Bekri, A. Juditsky and Arkadi Nemirovski 

close all
load out256_32.mat
s=out.s;
errx_h=out.errorsx(:,1);
errx_i=out.errorsx(:,2);
errx_hi=out.errorsx(:,3);

errn_h=out.errorsn(:,1);
errn_i=out.errorsn(:,2);
errn_hi=out.errorsn(:,3);
% 
optin_h=out.H.Opti;
optin_i=out.I.Opti;
optin_hi=out.HI.Opti;

risk_h=out.H.risk;
risk_i=out.I.risk;
risk_hi=out.HI.risk;

boxplot([errx_h,errx_i,errx_hi])
ymi=0.9*min(min([errx_h,errx_i,errx_hi])); yma=1.1*max([risk_h,risk_i,risk_hi]);
ylim([ymi yma])
ax=gca; ax.YAxis.Scale ="log";
ax.TickLabelInterpreter = 'latex';
ax.FontSize=14;
ax.XTickLabel={'$\widehat{x}_{HG}$','$\widehat{x}_{IG}$','$\widehat{x}_{HIG}$'};

hold on
line([0.75,1.25],[risk_h,risk_h],'color','red','linewidth',2)
line([1.75,2.25],[risk_i,risk_i],'color','red','linewidth',2)
line([2.75,3.25],[risk_hi,risk_hi],'color','red','linewidth',2)
grid on

riskn_h=optin_h*sqrt(2*s);
riskn_i=optin_i*sqrt(2*s);
riskn_hi=optin_hi*sqrt(2*s);
figure
boxplot([errn_h,errn_i,errn_hi])
ymi=0.9*min(min([errn_h,errn_i,errn_hi])); yma=1.1*max([riskn_h,riskn_i,riskn_hi]);
ylim([ymi yma])
ax=gca; ax.YAxis.Scale ="log";
ax.TickLabelInterpreter = 'latex';
ax.FontSize=14;
ax.XTickLabel={'$\widehat{\nu}_{H}$','$\widehat{\nu}_{I}$','$\widehat{\nu}_{HI}$'};

hold on
line([0.75,1.25],[riskn_h,riskn_h],'color','red','linewidth',2)
line([1.75,2.25],[riskn_i,riskn_i],'color','red','linewidth',2)
line([2.75,3.25],[riskn_hi,riskn_hi],'color','red','linewidth',2)
grid on

p=out.p; 
tt=2*pi/p*[1:p];
kk=24;
figure
emax=max(max([out.xhat0(:,kk)-out.xx(:,kk),out.xhat1(:,kk)-out.xx(:,kk),out.xhat2(:,kk)-out.xx(:,kk)]));
emin=min(min([out.xhat0(:,kk)-out.xx(:,kk),out.xhat1(:,kk)-out.xx(:,kk),out.xhat2(:,kk)-out.xx(:,kk)]));

plot(tt,out.xhat0(:,kk)-out.xx(:,kk),'--', tt,out.xhat1(:,kk)-out.xx(:,kk),'-.', tt,out.xhat2(:,kk)-out.xx(:,kk),'-','linewidth',1)
axis([1/p,2*pi,1.1*emin,1.1*emax])
grid on
legend('$\widehat x_{HG}-x_*$', '$\widehat x_{IG}-x_*$','$\widehat x_{HIG}-x_*$','FontSize', 12,'Interpreter','latex')
legend('boxoff')

figure
emax=max(max([out.xhat0(:,kk),out.xhat1(:,kk),out.xhat2(:,kk),out.xx(:,kk)]));
emin=min(min([out.xhat0(:,kk),out.xhat1(:,kk),out.xhat2(:,kk),out.xx(:,kk)]));

plot(tt,out.xhat0(:,kk),'--', tt,out.xhat1(:,kk),'-.', tt,out.xhat2(:,kk),'-','linewidth',1)
hold on
plot(tt,out.xx(:,kk),':','linewidth',2)
axis([1/p,2*pi,1.1*emin,1.1*emax])
grid on
legend('$\widehat x_{HG}$', '$\widehat x_{IG}$','$\widehat x_{HIG}$', '$x_*$','FontSize', 12,'Interpreter','latex')
legend('boxoff')
disp(['errx H=',num2str(errx_h(kk)),'; errx I=',num2str(errx_i(kk)), '; errx HI=',num2str(errx_hi(kk))])



