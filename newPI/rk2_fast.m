function [trange,x] = rk2_fast(f,trange,x)
N=length(trange)-1; %number of time steps
x=x';

for n=1:N
    t=trange(n);
    Dt=trange(n+1)-t;
    
    k1 =  Dt*feval(f,t,x)';
    k2 = Dt*feval(f,t+Dt,x+k1)';
    x=x+(k1+k2)/2;
end