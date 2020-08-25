% Computes the time derivatives of heterogeneous diffusion
% in 3D on patches.  Adapted from 2D heterogeneous diffusion.
% AJR, Aug 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff31()}: heterogeneous diffusion}
\label{sec:heteroDiff31}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 4D input arrays~\verb|u|, \verb|x| (via edge-value interpolation of
\verb|patchSmooth1|, \cref{sec:patchSmooth1}), computes the
time derivative at each
point in the interior of a patch, output in~\verb|ut|.  The
three 3D array of diffusivities,~$c^x_{ijk}$, $c^y_{ijk}$
and~$c^z_{ijk}$, have previously been stored
in~\verb|patches.cs| (2D). 
\begin{matlab}
%}
function ut = heteroDiff31(t,u,x)
  global patches
  dx = diff(x(2:3));  % x space step
  [nx,nV,nE,Nx] = size(u);
  ix = 2:nx-1; % x interior points in a patch
  ut = nan(size(u));  % preallocate output array
  ny = round(sqrt(nV)); % size of square cross-section
  c = reshape(patches.cs,nx-1,ny,ny,3,[]);
  u = reshape(u,[nx ny ny patches.nEnsem size(u,4)]);
  ut(ix,iy,iz,:,:,:) ...
  = diff(patches.cs(:,iy,iz,1,:).*diff(u(:,iy,iz,:,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,iz,2,:).*diff(u(ix,:,iz,:,:,:,:,:),1,2),1,2)/dx^2 ...
   +diff(patches.cs(ix,iy,:,3,:).*diff(u(ix,iy,:,:,:,:,:,:),1,3),1,3)/dx^2; 
end% function
%{
\end{matlab}
%}
