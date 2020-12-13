% rotFilmMicro() computes the time derivatives of
% heterogeneous advection-diffusion in 2D along a 1D channel
% on 1D array patches. AJR, Dec 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{rotFilmMicro()}: 2D shallow water flow on a rotating heterogeneous substrate}
\label{sec:rotFilmMicro}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 6D input arrays~\verb|huv| (via edge-value interpolation of
\verb|patchSmooth2|, \cref{sec:patchSmooth2}), computes the
time derivatives~\eqref{eqs:spinddt} at each point in the
interior of a patch, output in~\verb|huvt|.  The heterogeneous
bed drag and diffusivities,~$b_{ij}$
and~$\nu_{ij}$, have previously been merged and
stored in the array~\verb|patches.cs| (2D\({}\times3\)). 
\begin{matlab}
%}
function huvt = rotFilmMicro(t,huv,p)
  [nx,ny,~,~,~,~]=size(huv); % micro-grid points in patches
  i = 2:nx-1;          % x interior points in a patch
  j = 2:ny-1;          % y interior points in a patch
  dx = diff(p.x(2:3)); % x space step
  dy = diff(p.y(2:3)); % y space step
  huvt = nan+huv;      % preallocate output array
%{
\end{matlab}
Set indices of fields in the arrays.
Need to store different diffusivity values for the \(x,y\)-directions as they are evaluated at different points in space.
\begin{matlab}
%}
  h=1; u=2; v=3; 
  b=1; nux=2; nuy=3;
%{
\end{matlab}
Evaluate conservation of mass \pde~\eqref{eq:spindhdt}:
\begin{matlab}
%}
  huvt(i,j,h,:,:,:) = ...
    - (huv(i+1,j,h,:,:,:).*huv(i+1,j,u,:,:,:) ...
      -huv(i-1,j,h,:,:,:).*huv(i-1,j,u,:,:,:))/(2*dx) ...
    - (huv(i,j+1,h,:,:,:).*huv(i,j+1,v,:,:,:) ...
      -huv(i,j-1,h,:,:,:).*huv(i,j-1,v,:,:,:))/(2*dy) ;
%{
\end{matlab}
Evaluate the \(x\)-direction momentum \pde~\eqref{eq:spindvdt}:
\begin{matlab}
%}
  huvt(i,j,u,:,:,:) = ...
    - p.cs(i,j,b).*huv(i,j,u,:,:,:) + p.f.*huv(i,j,v,:,:,:) ...
    - huv(i,j,u,:,:,:).*(huv(i+1,j,u,:,:,:)-huv(i-1,j,u,:,:,:))/(2*dx) ...
    - huv(i,j,v,:,:,:).*(huv(i,j+1,u,:,:,:)-huv(i,j-1,u,:,:,:))/(2*dy) ...
    - p.g*(huv(i+1,j,h,:,:,:)-huv(i-1,j,h,:,:,:))/(2*dx) ... 
    + diff(p.cs(:,j,nux).*diff(huv(:,j,u,:,:,:),[],1),[],1)/dx^2 ...
    + diff(p.cs(i,:,nuy).*diff(huv(i,:,u,:,:,:),[],2),[],2)/dy^2 ;
%{
\end{matlab}
Evaluate the \(y\)-direction momentum \pde~\eqref{eq:spindvdt}:
\begin{matlab}
%}
  huvt(i,j,v,:,:,:) = ...
    - p.cs(i,j,b).*huv(i,j,v,:,:,:) - p.f.*huv(i,j,u,:,:,:) ...
    - huv(i,j,u,:,:,:).*(huv(i+1,j,v,:,:,:)-huv(i-1,j,v,:,:,:))/(2*dx) ...
    - huv(i,j,v,:,:,:).*(huv(i,j+1,v,:,:,:)-huv(i,j-1,v,:,:,:))/(2*dy) ...
    - p.g*(huv(i,j+1,h,:,:,:)-huv(i,j-1,h,:,:,:))/(2*dy) ... 
    + diff(p.cs(:,j,nux).*diff(huv(:,j,v,:,:,:),[],1),[],1)/dx^2 ...
    + diff(p.cs(i,:,nuy).*diff(huv(i,:,v,:,:,:),[],2),[],2)/dy^2 ;
end% function
%{
\end{matlab}
%}
