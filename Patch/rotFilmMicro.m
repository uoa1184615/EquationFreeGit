% rotFilmMicro() computes the time derivatives of a 2D
% shallow water flow on a rotating heterogeneous substrate
% on 2D  patches in space.  AJR, Dec 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{rotFilmMicro()}: 2D shallow water flow
on a rotating heterogeneous substrate}
\label{sec:rotFilmMicro}

This function codes the heterogeneous shallow water
flow~\eqref{eqs:spinddt} inside 2D patches. The \pde{}s are
discretised on the multiscale lattice in terms of evolving
variables~$h_{ijIJ}$, $u_{ijIJ}$ and~$v_{ijIJ}$.  For 6D
input array~\verb|huv| (via edge-value interpolation of
\verb|patchEdgeInt2()|, \cref{sec:patchSys2}), computes
the time derivatives~\eqref{eqs:spinddt} at each point in
the interior of a patch, output in~\verb|huvt|.  The
heterogeneous bed drag and diffusivities,~$b_{ij}$
and~$\nu_{ij}$, have previously been merged and stored in
the array~\verb|patches.cs| (2D${}\times3$): herein
\verb|patches| is named~\verb|p|. 
\begin{matlab}
%}
function huvt = rotFilmMicro(t,huv,p)
  [nx,ny,~]=size(huv); % micro-grid points in patches
  i = 2:nx-1;          % x interior points in a patch
  j = 2:ny-1;          % y interior points in a patch
  dx = diff(p.x(2:3)); % x space step
  dy = diff(p.y(2:3)); % y space step
  huvt = nan+huv;      % preallocate output array
%{
\end{matlab}
Set indices of fields in the arrays. Need to store different
diffusivity values for the $x,y$-directions as they are
evaluated at different points in space.
\begin{matlab}
%}
  h=1; u=2; v=3; 
  b=1; nux=2; nuy=3;
%{
\end{matlab}
Use a staggered micro-grid so that $\verb|h(i,j)| =h_{ij}$,
$\verb|u(i,j)| =u_{i+1/2,j}$, and $\verb|v(i,j)|
=v_{i,j+1/2}$.  We need the following to interpolate some
quantities to other points on the staggered micro-grid. But
the first two statements fill-in two needed corner values
because they are not (currently) interpolated by
\verb|patchEdgeInt2()|.
\begin{matlab}
%}
huv(1,ny,u,:,:,:) = huv(2,ny,u,:,:,:)+huv(1,ny-1,u,:,:,:) ...
                   -huv(2,ny-1,u,:,:,:);
huv(nx,1,v,:,:,:) = huv(nx,2,v,:,:,:)+huv(nx-1,1,v,:,:,:) ...
                   -huv(nx-1,2,v,:,:,:);
v4u = (huv(i,j-1,v,:,:,:)+huv(i+1,j,v,:,:,:) ...
      +huv(i,j,v,:,:,:)+huv(i+1,j-1,v,:,:,:))/4;
u4v = (huv(i,j+1,u,:,:,:)+huv(i-1,j,u,:,:,:) ...
      +huv(i,j,u,:,:,:)+huv(i-1,j+1,u,:,:,:))/4;
h2u = (huv(2:nx,:,h,:,:,:)+huv(1:nx-1,:,h,:,:,:))/2;
h2v = (huv(:,2:ny,h,:,:,:)+huv(:,1:ny-1,h,:,:,:))/2;
%{
\end{matlab}
Evaluate conservation of mass \pde~\eqref{eq:spindhdt}
(needing averages of~$h$ at half-grid points):
\begin{matlab}
%}
  huvt(i,j,h,:,:,:) = ...
    - (h2u(i,j  ,:,:,:,:).*huv(i  ,j,u,:,:,:) ...
      -h2u(i-1,j,:,:,:,:).*huv(i-1,j,u,:,:,:) )/dx ...
    - (h2v(i,j  ,:,:,:,:).*huv(i,j  ,v,:,:,:) ...
      -h2v(i,j-1,:,:,:,:).*huv(i,j-1,v,:,:,:))/dy ;
%{
\end{matlab}
Evaluate the $x$-direction momentum
\pde~\eqref{eq:spindvdt} (needing to interpolate
component~$v$ to $u$-points):
\begin{matlab}
%}
  huvt(i,j,u,:,:,:) = ...
    - p.cs(i,j,b).*huv(i,j,u,:,:,:) + p.f.*v4u ...
    - huv(i,j,u,:,:,:).*(huv(i+1,j,u,:,:,:)-huv(i-1,j,u,:,:,:))/(2*dx) ...
    - v4u.*(huv(i,j+1,u,:,:,:)-huv(i,j-1,u,:,:,:))/(2*dy) ...
    - p.g*(huv(i+1,j,h,:,:,:)-huv(i,j,h,:,:,:))/dx ... 
    + diff(p.cs(:,j,nux).*diff(huv(:,j,u,:,:,:),[],1),[],1)/dx^2 ...
    + diff(p.cs(i,:,nuy).*diff(huv(i,:,u,:,:,:),[],2),[],2)/dy^2 ;
%{
\end{matlab}
Evaluate the $y$-direction momentum
\pde~\eqref{eq:spindvdt} (needing to interpolate
component~$u$ to $v$-points):
\begin{matlab}
%}
  huvt(i,j,v,:,:,:) = ...
    - p.cs(i,j,b).*huv(i,j,v,:,:,:) - p.f.*u4v ...
    - u4v.*(huv(i+1,j,v,:,:,:)-huv(i-1,j,v,:,:,:))/(2*dx) ...
    - huv(i,j,v,:,:,:).*(huv(i,j+1,v,:,:,:)-huv(i,j-1,v,:,:,:))/(2*dy) ...
    - p.g*(huv(i,j+1,h,:,:,:)-huv(i,j,h,:,:,:))/dy ... 
    + diff(p.cs(:,j,nux).*diff(huv(:,j,v,:,:,:),[],1),[],1)/dx^2 ...
    + diff(p.cs(i,:,nuy).*diff(huv(i,:,v,:,:,:),[],2),[],2)/dy^2 ;
end% function
%{
\end{matlab}
%}
