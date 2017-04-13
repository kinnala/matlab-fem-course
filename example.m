clear all
[p,e,t]=initmesh(decsg([3 4 0 1 1 0 0 0 1 1]'),'Hmax',0.1);
t=t(1:3,:);
a=LinearAssembler(p,t);
K=a.assembleBilinear(@(u,du,v,dv,x) du{1}.*dv{1}+du{2}.*dv{2});
f=a.assembleLinear(@(v,dv,x) 1.0*v);
N=size(K,2);
D=unique(e(1:2,:));
I=setdiff(1:N,D);
x=zeros(N,1);
x(I)=K(I,I)\f(I);
figure;
trisurf(t',p(1,:),p(2,:),x)