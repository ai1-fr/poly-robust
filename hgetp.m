function [tcpu,crs1,crs0] = hgetp(A, s, cntrl)
% HGETP: optimization of H-component of the contrast 
% experiment in the paper "On robust recovery of signals from indirect observations" 
% by Y. Bekri, A. Juditsky and Arkadi Nemirovski 
%Call: [tcpu,crs1,crs0] = hgetp(A, s, cntrl) 
% 
% Inputs: straightforward 
%   For the composition of the cntrl structure see p2025sim.m
% Outputs: 
%   tcpu: elapsed cpu time
%   crs1: structure of the optimized contrast H with fields
%       H:    optimized H contrast component
%       errh: l2-norms of columns of H
%       div:  the denominator 1-2s\kappa of the risk bound for l1-recovery
%   crs0: structure of the optimized favtorizing contrast H' with fields
%       H:    optimized H contrast component
%       errh: l2-norms of columns of H
%       div:  the denominator 1-2s\kappa of the risk bound for l1-recovery
%
% by A. Juditsky and A. Nemirovski, Jan. 2025 
%  
tt=cputime;
ni=nargin;
no=nargout;

if ni<3, error('Estimate parameters [CNTRL] undefined'); end
sigm=cntrl.sigma;
eeps=cntrl.eps;
L=cntrl.L;
if isfield(cntrl,'div')
    ddenom=cntrl.div;
else 
    ddenom=1/4;
    cntrl.div=ddenom;
end



% if ni<4, L=[1,1,1]; end
% if ni<5, ddenom=1/4; end
% vkapp=sigm*sqrt(2)*erfcinv(eeps/Ne);

%cvx_solver('mosek')
[m,p]=size(A);
dd=2*pi/p;
Ne=2*m;
n=m;
Ne=2*n;
vkapp=sigm*sqrt(2)*erfcinv(eeps/Ne);
T=getT(p, dd, L);

H=zeros(m,m);
maxkap=(1-ddenom)/2/s;
eh2=zeros(n,1);

% compute "standard contrast
for i=1:n
    ei=zeros(n,1); ei(i)=1;
    % contrasts
    cvx_begin quiet
    variables h(m,1) t xi(p,1)
%    
    minimize t
%
    subject to
    norm(h,2)*vkapp+norm(xi,1) <= t;
    norm(ei-h,Inf) <= maxkap;
    A'*h-T*xi == 0;
    cvx_end
    eh2(i)=t;
    H(:,i)=h;
end
crs1.H=H;
crs1.errh=eh2;
crs1.div=1-2*s*maxkap;

if no>2
    %compute "factorizing contrast"

    % check maxkap value first
    find_mk=0;
    while ~find_mk
        div_mk=0;
        for i=1:n
            ei=zeros(n,1); ei(i)=1;
            % contrasts
            cvx_begin quiet
            variables h(m,1)

            minimize norm(ei-h,Inf)

            subject to
            A'*h == 0;
            cvx_end
            if cvx_optval>maxkap
                maxkap=cvx_optval+1e-4;
                div_mk=1;
                break;
            end
        end
        if ~div_mk
            find_mk=1;
        end
    end
    for i=1:n
        ei=zeros(n,1); ei(i)=1;
        % contrasts
        cvx_begin quiet
        variables h(m,1) t

        minimize t

        subject to
        norm(h,2)*vkapp <= t;
        norm(ei-h,Inf) <= maxkap;
        A'*h == 0;
        cvx_end
        eh2(i)=t;
        H(:,i)=h;
    end
    crs0.H=H;
    crs0.errh=eh2;
    crs0.div=1-2*s*maxkap;
end
tcpu=cputime-tt;
end %endof hget
