% A microscale discretisation of Burgers' PDE on a lattice x.  
% AJR 5 Apr 2019 -- Jun 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\paragraph{Example of Burgers PDE inside patches}
As a microscale discretisation of Burgers' \pde\ 
\(u_t=u_{xx}-30uu_x\), here code \(\dot u_{ij} 
=\frac1{\delta x^2} (u_{i+1,j}-2u_{i,j}+u_{i-1,j}) 
-30u_{ij} \frac1{2\delta x}(u_{i+1,j}-u_{i-1,j})\).
Here there is only one field variable, and one in the
ensemble, so for simpler coding of the PDE we squeeze them
out (with no need to reshape when via patchSmooth1()).
\begin{matlab}
%}
function ut=BurgersPDE(t,u,patches)
  u=squeeze(u);     % omit singleton dimensions
  dx=diff(patches.x(1:2));  % microscale spacing
  i=2:size(u,1)-1;  % interior points in patches
  ut=nan+u;         % preallocate output array
  ut(i,:)=diff(u,2)/dx^2 ...
    -30*u(i,:).*(u(i+1,:)-u(i-1,:))/(2*dx);
end
%{
\end{matlab}
%}
