% Computes the time derivatives of heterogeneous diffusion
% in 2D on patches with macroscale BCs.  AJR, Sep 2023 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff20()}: heterogeneous diffusion}
\label{sec:heteroDiff20}

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
function ut = heteroDiff20(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
The macroscale boundary conditions are Dirichlet zero at the
extreme edges of left-right extreme patches, and Neumann zero 
at extreme edges of top-bottom extreme patches.
\begin{matlab}
%}
  u( 1 ,:,:,:, 1 ,:)=0; % left-edge of leftmost is zero
  u(end,:,:,:,end,:)=0; % right-edge of rightmost is zero
  u(:, 1 ,:,:,:, 1 )=u(:,  2  ,:,:,:, 1 ); % bottom-edge of bottommost
  u(:,end,:,:,:,end)=u(:,end-1,:,:,:,end); % top-edge of topmost
%{
\end{matlab}
Code the microscale diffusion.
\begin{matlab}
%}
  ut(ix,iy,:,:,:,:) ...
  = diff(patches.cs(:,iy,1,:).*diff(u(:,iy,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,2,:).*diff(u(ix,:,:,:,:,:),1,2),1,2)/dy^2; 
end% function
%{
\end{matlab}
%}
