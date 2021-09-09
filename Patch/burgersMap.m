% Microscale Euler step of the Burgers PDE on a lattice in
% x.  Used by BurgersExample.m AJR, 4 Apr 2019 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{burgersMap()}: discretise the PDE microscale}
\label{sec:burgersMap}
This function codes the microscale Euler integration map of
the lattice differential equations inside the patches.  Only
the patch-interior values are mapped (\verb|patchSys1()|
overrides the edge-values anyway).
\begin{matlab}
%}
function u = burgersMap(t,u,patches)
  u = squeeze(u);
  dx = diff(patches.x(2:3));   
  dt = dx^2/2;
  i = 2:size(u,1)-1;
  u(i,:) = u(i,:) +dt*( diff(u,2)/dx^2 ...
     -20*u(i,:).*(u(i+1,:)-u(i-1,:))/(2*dx) );
end
%{
\end{matlab}
%}
