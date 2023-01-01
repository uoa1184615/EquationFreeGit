
%{ 
Testing some divided difference interpolation \cite{DividedDifferences,
PolynomialInterpolation} AJR, Dec 2022.   Patch scheme:  place end-patch
edges at the boundary: for micro-Dirichlet place micro-grid point at bdry;
for micro-Neumann place bdry halfway between; set one or other.  Maybe use
a number with 0 being Dirichlet and 0.5 being Neumann.  Allow to choose the
patches equi-spaced, or Chebyshev points, or prescribed centre locations.
Possibly allow centre locations to vary in time, via displacements, as an
optinal extra---unclear how to work in multi-D so maybe leave out for now.
%}
n=6
if 0, x=sort(pi*rand(n,1)), f=sin(x), theCase=1
else  x=sort(cos(linspace(0,pi,n)')), f=exp(-x)*0+(1:n)'==3, theCase=0
end
% divided differences in triangular table
F=nan(n);
F(:,1)=f;
for j=1:n-1
  i=1:n-j;
  F(i,j+1)=(F(i+1,j)-F(i,j))./(x(i+j)-x(i));
end
F=F
% or in the Matrix Form, that seems to hard to use
df=f;
Tf=diag(df);
for j=1:n-1
  i=1:n-j;
  df=diff(df)./(x(i+j)-x(i));
  Tf=Tf+diag(df,j);
end
Tf=Tf
% Horner's evaluation of polynomial at points in xs
if theCase, xs=linspace(0,pi);
else xs=linspace(min(x)-0.1,max(x)+0.1);
end
fs=repmat(F(1,n),size(xs));
for i=n-1:-1:1
  fs=F(1,i)+(xs-x(i)).*fs;
end
plot(x,f,'o',xs,fs)
xlabel('x'),ylabel('f')
% what if we have n points---near the collocation points say
% and we use pth order interpolation, near-centred
for z=linspace(-0.5,0.5,10)
xs=x+z*[diff(x);0.2];
%xs=x([2:n n]);
p=2 % need order p < n the number of points
% i = index of left end of interpolation
i=max(1,min(1:n,n-ceil(p/2))-floor(p/2))
fs=F(i,p+1);
for q=p:-1:1
  fs=F(i,q)+(xs-x(i+q-1)).*fs;
end
hold on, plot(xs,fs,'.'), hold off
end
