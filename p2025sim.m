% Simulation experiment for the paper "On robust recovery of signals from
% indirect observations" by Y. Bekri, A. Juditsky and Arkadi Nemirovski 
%
% set simulation parameters
nsim=100; % # of repeated simulations of observation \omega  
nuss=10;  % sigma of nonvanishing entries in nu
cntrl.sigma=0.1;
cntrl.eps=0.05;
cntrl.L=[4,1,4];
eeps=cntrl.eps;
sigm=cntrl.sigma;
s=8;
A=randn(256,32);
[m,p]=size(A);
for i=1:p
    A(:,i)=A(:,i)/norm(A(:,i),2);
end
n=m;
cvx_solver('mosek')

% compute optimized contrasts
[CtrH, cntrl]=GHget(A, s, cntrl);
optiH=CtrH.optv.Opti;
[CtrI, cntrl]=GIget(A, s, cntrl);
optiI=CtrI.optv.Opti;

CGH=[CtrH.H,CtrH.G];
CGI=[CtrI.H,CtrI.G];

% compute corresponding error bounds
[tt1,errbH] = Gevalp(CGH, A, s,optiH, cntrl);
[tt2,errbI] = Gevalp(CGI, A, s,optiI, cntrl);

% evaluate bound for the "composite" contrast
CGHI=[CtrH.H,CtrH.G,CtrI.H,CtrI.G];  
cntrl.eps=cntrl.eps*2*m/(size(CGHI,2));
[tt3,OptiHI] = hevalp(CGHI, A, s, cntrl);
optiHI=min([optiI,optiH,OptiHI]);
[tt4,errbHI] = Gevalp(CGHI, A, s,optiHI, cntrl);

% run simulations

%generate x signals
Rn = 2*binornd(1,0.5,p,nsim)-1;
Rn(1:2,:)=rand(2,nsim).*Rn(1:2,:);
dh=2*pi/p;
T=getT(p,dh,cntrl.L);
xx=T'\Rn;

%generate nu signals
ik=zeros(s,nsim);
tnu=zeros(n,nsim);
for i=1:nsim
    ik=randperm(n,s);
    tnu(ik,i)=nuss*randn(s,1);
end
% generate observation omega
omeg=A*xx+tnu+cntrl.sigma*randn(m,nsim);

% estimate with optimized H and G
Nc=size(CGH,2);
vkapp=sigm*sqrt(2)*erfcinv(eeps/Nc/2);

hn2=sqrt(sum(CGH.^2,1)); 
rhsh=vkapp*hn2;
xhat0=zeros(p,nsim);
nuhat0=zeros(n,nsim);

% compute estimates
for k=1:nsim
    cvx_begin quiet
    variables nnu(n,1) z(p,1)
    minimize(norm(nnu,1))
    subject to
    abs(T'*z)<=1;
    for i=1:Nc
        abs(CGH(:,i)'*(nnu+A*z-omeg(:,k))) <=rhsh(i);
    end
    cvx_end
    nuhat0(:,k)=nnu;
    xhat0(:,k)=z;
end


% estimate with contrasts I and G'
Nc=size(CGI,2);
vkapp=sigm*sqrt(2)*erfcinv(eeps/Nc/2);

hn2=sqrt(sum(CGI.^2,1)); 
rhsh=vkapp*hn2;
xhat1=zeros(p,nsim);
nuhat1=zeros(n,nsim);

% compute estimates
for k=1:nsim
    cvx_begin quiet
    variables nnu(n,1) z(p,1)
    minimize(norm(nnu,1))
    subject to
    abs(T'*z)<=1;
    for i=1:Nc
        abs(CGI(:,i)'*(nnu+A*z-omeg(:,k))) <=rhsh(i);
    end
    cvx_end
    nuhat1(:,k)=nnu;
    xhat1(:,k)=z;
end


% estimate with composite [H,I,G,G']
Nc=size(CGHI,2);
vkapp=sigm*sqrt(2)*erfcinv(eeps/Nc/2);

hn2=sqrt(sum(CGHI.^2,1)); 
rhsh=vkapp*hn2;
xhat2=zeros(p,nsim);
nuhat2=zeros(n,nsim);

% compute estimates
for k=1:nsim
    cvx_begin quiet
    variables nnu(n,1) z(p,1)
    minimize(norm(nnu,1))
    subject to
    abs(T'*z)<=1;
    for i=1:Nc
        abs(CGHI(:,i)'*(nnu+A*z-omeg(:,k))) <=rhsh(i);
    end
    cvx_end
    nuhat2(:,k)=nnu;
    xhat2(:,k)=z;
end
% dereference outputs
out.s=s;
out.xx=xx;
out.xhat0=xhat0;
out.xhat1=xhat1;
out.xhat2=xhat2;

out.errorsx(:,1)=sqrt(sum((xhat0-xx).^2,1))';
out.errorsn(:,1)=sqrt(sum((nuhat0-tnu).^2,1))';

out.errorsx(:,2)=sqrt(sum((xhat1-xx).^2,1))';
out.errorsn(:,2)=sqrt(sum((nuhat1-tnu).^2,1))';

out.errorsx(:,3)=sqrt(sum((xhat2-xx).^2,1))';
out.errorsn(:,3)=sqrt(sum((nuhat2-tnu).^2,1))';


out.H.Opti=optiH;
out.I.Opti=optiI;
out.HI.Opti=optiHI;
out.H.risk=errbH;
out.I.risk=errbI;
out.HI.risk=errbHI;
out.p=p;
out.m=m;
out.n=n;
out.size.GH=size(CGH,2);
out.size.GI=size(CGI,2);
out.size.GHI=size(CGHI,2);

save('out256_32.mat', 'out')

