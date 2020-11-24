% Microscale discretisation of a nonlinear diffusion PDE in
% 2D space (x,y) in 2D patches.
% AJR, 5 Apr 2019 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\paragraph{Example of nonlinear diffusion PDE inside patches}
As a microscale discretisation of \(u_t=\delsq(u^3)\), code
\(\dot u_{ijkl} =\frac1{\delta x^2} (u_{i+1,j,k,l}^3
-2u_{i,j,k,l}^3 +u_{i-1,j,k,l}^3) + \frac1{\delta y^2}
(u_{i,j+1,k,l}^3 -2u_{i,j,k,l}^3 +u_{i,j-1,k,l}^3)\).
\begin{matlab}
%}
function ut = nonDiffPDE(t,u,patches)
  if nargin<3, global patches, end
  u = squeeze(u); % reduce to 4D
  dx = diff(patches.x(1:2));  % microgrid spacing
  dy = diff(patches.y(1:2));
  i = 2:size(u,1)-1;  j = 2:size(u,2)-1;  % interior patch points
  ut = nan+u;  % preallocate output array
  ut(i,j,:,:) = diff(u(:,j,:,:).^3,2,1)/dx^2 ...
               +diff(u(i,:,:,:).^3,2,2)/dy^2;
end
%{
\end{matlab}
%}
