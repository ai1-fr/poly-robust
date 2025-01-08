function [tcpu,binf] = hevalp(H, A, s, cntrl)
%HEVALP: improved evaluation of the accuracy of l_1-recovery for given contrast H 
% experiment in the paper "On robust recovery of signals from indirect observations" 
% by Y. Bekri, A. Juditsky and Arkadi Nemirovski 
% Call:  [tcpu,binf] = hevalp(H, A, s, cntrl)
% 
% Inputs: straightforward 
%   For the composition of the cntrl structure see p2025sim.m
% Outputs: structure with fields
%   tcpu: elapsed cpu time
%   binf: bound for inf-norm of the error 
%
% by A. Juditsky and A. Nemirovski, Jan. 2025   

% no=nargout;
tcpu=cputime;
ni=nargin;
%  
if ni<4, error('Estimate parameters [CNTRL] undefined'); end
sigm=cntrl.sigma;
eeps=cntrl.eps;
L=cntrl.L;

nH2=sqrt(sum(H.^2,1));
nH2=nH2(:);
[m,p]=size(A);
n=m;
Ne=2*n;
vkapp=sigm*sqrt(2)*erfcinv(eeps/Ne);

T=getT(p, 2*pi/p, L);

% cvx_solver('mosek')
% compute l_infty bound
ilinf=zeros(n,1);
for i=1:n
    ei=zeros(n,1); ei(i)=1;
    cvx_begin quiet
    variables w(n,1) v(p,1) u
    maximize u
    subject to
        u<=ei'*w;
        norm(w,Inf)<=u;
        norm(w,1)<=2*s*u;
        abs(H'*w+H'*A*v)<=2*vkapp*nH2;
        abs(T'*v)<=ones(p,1);
    cvx_end
    ilinf(i)=u;
end
binf=max(ilinf);
tcpu=cputime-tcpu;
end 
%endof as_eval
