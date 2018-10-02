function [trange,x] = rk2(f,trange,xtmp)
N=length(trange)-1; %number of time steps

x=zeros(N,length(xtmp));
x(1,:)=xtmp';

for n=1:N
    t=trange(n);
    Dt=trange(n+1)-t;
    k1 =  Dt*feval(f,t,xtmp);
    k2 = Dt*feval(f,t+Dt,xtmp+k1);
    xtmp=xtmp+(k1+k2)/2;
    x(n+1,:)=xtmp';
end