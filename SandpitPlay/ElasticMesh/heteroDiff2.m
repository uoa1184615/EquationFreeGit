% Computes the time derivatives of heterogeneous diffusion
% in 2D on patches.  Adapted from 1D heterogeneous diffusion.
% JEB & AJR, May 2020 -- Nov 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff2()}: heterogeneous diffusion}
\label{sec:heteroDiff2}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 6D input arrays~\verb|u|, \verb|x|,
and~\verb|y| (via edge-value interpolation of
\verb|patchSys2|, \cref{sec:patchSys2}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
two 2D array of diffusivities,~$c^x_{ij}$ and~$c^y_{ij}$,
have previously been stored in~\verb|patches.cs| (3D). 
\begin{matlab}
%}
function ut = heteroDiff2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
  ut(ix,iy,:,:,:,:) ...
  = diff(patches.cs(:,iy,1,:).*diff(u(:,iy,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,2,:).*diff(u(ix,:,:,:,:,:),1,2),1,2)/dy^2; 
end% function
%{
\end{matlab}
%}
