function  [tcpu,errbnd]= Gevalp(G, A, s, Opti, cntrl)
%GEVALP:bounding the risk of the estimate with contrast G
% experiment in the paper "On robust recovery of signals from indirect observations" 
% by Y. Bekri, A. Juditsky and Arkadi Nemirovski 
%Call:  [tcpu,errbnd]= Gevalp(G, A, s, Opti, cntrl)
% 
% Inputs: straightforward 
%   For the composition of the cntrl structure see p2025sim.m
% Outputs: 
%   tcpu: elapsed cpu time
%   errbnd: computed risk bound
%
% by A. Juditsky and A. Nemirovski, Jan. 2025 
tt=cputime;
ni=nargin;
[m,p]=size(A);
pG=size(G,2);
n=m;
if ni<5, error('Estimate parameters [CNTRL] undefined'); end
sigm=cntrl.sigma;
eeps=cntrl.eps;
L=cntrl.L;
Ne=2*n;
vkapp=sigm*sqrt(2)*erfcinv(eeps/Ne);

barA=G'*A;
tau=zeros(pG,1);
for k=1:pG
    tau(k)=sum_largest(G(:,k),2*s)*Opti +2*vkapp*norm(G(:,k),2);
end
dh=2*pi/p;
T=getT(p,dh,L);

cvx_begin sdp quiet

variables lam mmu(p,1) gam(pG,1)
variable V(p,p) symmetric

minimize(lam+4*sum(mmu)+gam'*(tau.^2))
subject to
mmu >= zeros(p,1);
lam >= 0;
gam >= zeros(pG,1);

T*diag(mmu)*T'+barA'*diag(gam)*barA-V >= 0;
[eye(p)*lam, eye(p)/2; eye(p)/2, V] >= 0;
cvx_end
errbnd=cvx_optval;
tcpu=cputime-tt;
end
