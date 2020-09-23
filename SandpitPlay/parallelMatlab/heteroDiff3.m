% Computes the time derivatives of heterogeneous diffusion
% in 3D on patches.  Adapted from 2D heterogeneous diffusion.
% JEB & AJR, May 2020 -- Sep 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff3()}: heterogeneous diffusion}
\label{sec:heteroDiff3}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 8D input array~\verb|u| (via edge-value interpolation of
\verb|patchSmooth3|, \cref{sec:patchSmooth3}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
three 3D array of diffusivities,~$c^x_{ijk}$, $c^y_{ijk}$
and~$c^z_{ijk}$, have previously been stored
in~\verb|patches.cs| (4D). 
\begin{matlab}
%}
function ut = heteroDiff3(t,u,patches)
  dx = diff(patches.x(2:3));  % x micro-scale step
  dy = diff(patches.y(2:3));  % y micro-scale step
  dz = diff(patches.z(2:3));  % z micro-scale step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  iz = 2:size(u,3)-1; % y interior points in a patch
  ut = nan+u; % reserve storage
%  Above is quicker and has much less communication than this
%  if ~patches.parallel, ut=nan(size(u));
%  else ut=nan(size(u),patches.codist); end
  ut(ix,iy,iz,:,:,:,:,:) ...
  = diff(patches.cs(:,iy,iz,1,:).*diff(u(:,iy,iz,:,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,iz,2,:).*diff(u(ix,:,iz,:,:,:,:,:),1,2),1,2)/dy^2 ...
   +diff(patches.cs(ix,iy,:,3,:).*diff(u(ix,iy,:,:,:,:,:,:),1,3),1,3)/dz^2; 
end% function
%{
\end{matlab}
%}
