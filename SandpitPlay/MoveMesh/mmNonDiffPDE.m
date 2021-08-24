% Microscale discretisation of a nonlinear diffusion PDE in
% 2D space (x,y) on moving 2D patches.
% AJR, Aug 2021
%!TEX root = doc.tex
%{
\subsection{\texttt{mmNonDiffPDE():} nonlinear diffusion PDE inside moving patches}
As a microscale discretisation of \(u_t=\Vv\cdot\grad u+\delsq(u^3)\), code
\(\dot u_{ijkl} =\cdots+\frac1{\delta x^2} (u_{i+1,j,k,l}^3
-2u_{i,j,k,l}^3 +u_{i-1,j,k,l}^3) + \frac1{\delta y^2}
(u_{i,j+1,k,l}^3 -2u_{i,j,k,l}^3 +u_{i,j-1,k,l}^3)\).
\begin{matlab}
%}
function ut = mmNonDiffPDE(t,u,M,patches)
  if nargin<3, global patches, end
  u = squeeze(u); % reduce to 4D
  Vx = shiftdim(M.Vx,2); % omit two singleton dimens
  Vy = shiftdim(M.Vy,2); % omit two singleton dimens
  dx = diff(patches.x(1:2));  % microgrid spacing
  dy = diff(patches.y(1:2));
  i = 2:size(u,1)-1;  j = 2:size(u,2)-1;  % interior patch points
  ut = nan+u;  % preallocate output array
  ut(i,j,:,:) = ...
      +Vx.*(u(i+1,j,:,:)-u(i-1,j,:,:))/(2*dx) ...
      +Vy.*(u(i,j+1,:,:)-u(i,j-1,:,:))/(2*dy) ...
      +diff(u(:,j,:,:).^3,2,1)/dx^2 ...
      +diff(u(i,:,:,:).^3,2,2)/dy^2 ;
end
%{
\end{matlab}
%}
