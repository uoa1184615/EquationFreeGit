% Computes the time derivatives of heterogeneous waves
% in 3D on patches.  AJR, Aug--Sep 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroWave3()}: heterogeneous Waves}
\label{sec:heteroWave3}

This function codes the lattice heterogeneous waves
inside the patches.  The wave \pde\ is
\begin{equation*}
u_t=v,\quad v_t=\grad(C\divv u)
\end{equation*}
for diagonal matrix~\(C\) which has microscale variations.
For 8D input arrays~\verb|u|, \verb|x|, \verb|y|,
and~\verb|z| (via edge-value interpolation of
\verb|patchSmooth3|, \cref{sec:patchSmooth3}), computes the
time derivative at each point in the interior of a patch,
output in~\verb|ut|.  The three 3D array of heterogeneous
coefficients,~$c^x_{ijk}$, $c^y_{ijk}$ and~$c^z_{ijk}$, have
previously been stored in~\verb|patches.cs| (4D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = heteroWave3(t,u,patches)
  if nargin<3, global patches, end
%{
\end{matlab}
Microscale space-steps, and interior point indices.  
\begin{matlab}
%}
  dx = diff(patches.x(2:3));  % x micro-scale step
  dy = diff(patches.y(2:3));  % y micro-scale step
  dz = diff(patches.z(2:3));  % z micro-scale step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  k = 2:size(u,3)-1; % z interior points in a patch
%{
\end{matlab}
Reserve storage and then assign interior patch values to the
heterogeneous diffusion time derivatives. Using \verb|nan+u|
appears quicker than \verb|nan(size(u),patches.codist)|
\begin{matlab}
%}
  ut = nan+u;  % preallocate output array
  ut(i,j,k,1,:) = u(i,j,k,2,:);
  ut(i,j,k,2,:) ...
  =diff(patches.cs(:,j,k,1,:).*diff(u(:,j,k,1,:),1),1)/dx^2 ...
  +diff(patches.cs(i,:,k,2,:).*diff(u(i,:,k,1,:),1,2),1,2)/dy^2 ...
  +diff(patches.cs(i,j,:,3,:).*diff(u(i,j,:,1,:),1,3),1,3)/dz^2; 
end% function
%{
\end{matlab}
%}
