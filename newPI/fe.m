function [trange,x] = fe(f,trange,IC)
N=length(trange)-1; %number of time steps

x=zeros(N,length(IC));
x(1,:)=IC;

for n=1:N
   Dt=trange(n+1)-trange(n); 
   x(n+1,:) = x(n,:) + Dt*feval(f,trange(n),x(n,:)')';
end