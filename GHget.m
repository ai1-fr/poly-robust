function [Ctr, cntrl]=GHget(A, s, cntrl)
%GHGET: optimization of contrasts H and G
% experiment in the paper "On robust recovery of signals from indirect observations" 
% by Y. Bekri, A. Juditsky and Arkadi Nemirovski 

%Call:  [Ctr, cntrl]=GHget(A, s, cntrl)
% 
% Inputs: straightforward 
%   For the composition of the cntrl structure see p2025sim.m
% Outputs: structure with fields
%   H: H contrast component
%   G: G contrast component
%   optv: error bound of recovery of nuisance signal
%   errbnd: optimized error bound for recovery of x
%
% by A. Juditsky and A. Nemirovski, Jan. 2025 

% derefrence inputs
ni=nargin;

[m,p]=size(A);
n=m;
if ni<3, error('Estimate parameters [CNTRL] undefined'); end
sigm=cntrl.sigma;
eeps=cntrl.eps;
L=cntrl.L;
if isfield(cntrl,'div')
    hdenom=cntrl.div;
else 
    hdenom=1/4;
    cntrl.div=hdenom;
end

if isfield(cntrl,'gthr')
    Gthrsh=cntrl.gthr;
else 
    Gthrsh=1e-6;
    cntrl.gthr=Gthrsh;
end
Ne=2*n;
vkapp=sigm*sqrt(2)*erfcinv(eeps/Ne);
dh=2*pi/p;
% hdenom=1/4; % value of 1-2s\kapppa

% compute "standard" H contrasts
[tcpuH,ctr1] = hgetp(A, s, cntrl);
%evaluate esimation risks
[tcpu1,Opti] = hevalp(ctr1.H, A, s, cntrl);
optv.Opti=Opti;
optv.Opt1=2*s*Opti;
optv.Opt2=sqrt(2*s)*Opti;

qs=2*vkapp+optv.Opt2;
barA=A/qs;
hA=ctr1.H'*A/qs;

% compute Theta matrix
T=getT(p,dh,L);
tt=cputime;

cvx_begin sdp quiet
variable Theta(n,n) symmetric
variable V(p,p) symmetric
variables mmu(p,1) lam nnu(n,1)
minimize(lam+4*sum(mmu)+sum(nnu)+trace(Theta))
subject to
mmu >= zeros(p,1);
nnu >= zeros(n,1);
lam >= 0;
Theta >= 0;
T*diag(mmu)*T'+hA'*diag(nnu)*hA-V >=0;
[eye(p)*lam, eye(p)/2; eye(p)/2, barA'*Theta*barA+V] >=0;
cvx_end

erisk1=cvx_optval;
[W,D]=eig(Theta);
G=W(:,diag(D) >= Gthrsh);

% dereference outputs
Ctr.optv=optv;
Ctr.H=ctr1.H;
Ctr.G=G;
Ctr.errbnd=erisk1;
Ctr.cpuh=tcpuH;
Ctr.cpuopti=tcpu1;
Ctr.cpug=cputime-tt;

end %endof GHget

% pG=size(G,2);

% Q2i=(8*vkapp^2+4*s*Opti^2)*eye(m);
%Q=eye(m)/sqrt(8*vkapp^2+4*s*Opti^2);

% % create M_j matrices
% M0=8*vkapp^2*eye(m);
% Mj=zeros(m,m,m);
% for j=1:m
%     Mj(:,:,j)=M0;
%     Mj(j,j,j)=Mj(j,j,j)+8*s^2*Opti^2;
% end
