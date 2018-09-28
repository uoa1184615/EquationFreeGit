function [trange,x] = rk4(f,trange,xtmp)
N=length(trange)-1; %number of time steps

x=zeros(N,length(xtmp));
x(1,:)=xtmp;

for n=1:N
    t=trange(n);
    Dt=trange(n+1)-t;
    k1 =  Dt*feval(f,t,xtmp);
    k2 = Dt*feval(f,t+Dt/2,xtmp+k1/2);
    k3 = Dt*feval(f,t+Dt/2,xtmp+k2/2);
    k4 = Dt*feval(f,t+Dt,xtmp+k3);
    xtmp=xtmp+(k1+2*k2+2*k3+k4)/6;
    x(n+1,:)=xtmp';
end