function t=getT(p,h,L)
% constraint |t(i,:)*x|<=L

t=zeros(p,p);
t(1,1)=1/L(1);
t(1:2,2)=[-1;1]/h/L(2);
for i=3:p
    t(i-2:i,i)=[1;-2;1]/h^2/L(3);
end

end %endof getT