% heteroDiff3() computes the time derivatives of
% heterogeneous diffusion in 3D on patches.  Adapted from 2D
% heterogeneous diffusion. JEB & AJR, May--Sep 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff3()}: heterogeneous diffusion}
\label{sec:heteroDiff3}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 8D input array~\verb|u| (via
edge-value interpolation of \verb|patchEdgeInt3|, such as by
\verb|patchSmooth3|, \cref{sec:patchSmooth3}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
three 3D array of diffusivities,~$c^x_{ijk}$, $c^y_{ijk}$
and~$c^z_{ijk}$, have previously been stored
in~\verb|patches.cs| (4+D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = heteroDiff3(t,u,patches)
  if nargin<3, global patches, end
%{
\end{matlab}
Microscale space-steps.  
Q: is using \verb|ix,iy,iz| slower than \verb|2:end-1|??
\begin{matlab}
%}
  dx = diff(patches.x(2:3));  % x micro-scale step
  dy = diff(patches.y(2:3));  % y micro-scale step
  dz = diff(patches.z(2:3));  % z micro-scale step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  iz = 2:size(u,3)-1; % y interior points in a patch
%{
\end{matlab}
Reserve storage and then assign interior patch values to the
heterogeneous diffusion time derivatives. Using \verb|nan+u|
appears quicker than \verb|nan(size(u),patches.codist)|
\begin{matlab}
%}
  ut = nan+u; % reserve storage
  ut(ix,iy,iz,:,:,:,:,:) ...
  = diff(patches.cs(:,iy,iz,1,:).*diff(u(:,iy,iz,:,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,iz,2,:).*diff(u(ix,:,iz,:,:,:,:,:),1,2),1,2)/dy^2 ...
   +diff(patches.cs(ix,iy,:,3,:).*diff(u(ix,iy,:,:,:,:,:,:),1,3),1,3)/dz^2; 
end% function
%{
\end{matlab}
%}
