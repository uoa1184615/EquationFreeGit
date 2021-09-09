% chanDispMicro() computes the time derivatives of
% heterogeneous advection-diffusion in 2D along a long thin
% channel on 1D array patches. AJR, Nov 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{chanDispMicro()}: heterogeneous 2D
advection-diffusion in a long thin channel}
\label{sec:chanDispMicro}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 4D input arrays of
concentration~\verb|c| and spatial lattice~\verb|x| (via
edge-value interpolation of \verb|patchSys1|,
\cref{sec:patchSys1}), computes the time
derivative~\eqref{eq:ddeChanDisp} at each point in the
interior of a patch, output in~\verb|ct|.  The heterogeneous
advections and diffusivities,~$u_i(y_j)$
and~$\kappa_i(y_{j+1/2})$, have previously been merged and
stored in the one array~\verb|patches.cs| (2D). 
\begin{matlab}
%}
function ct = chanDispMicro(t,c,p)
  [nx,ny,~,~]=size(c); % micro-grid points in patches
  ix = 2:nx-1;         % x interior points in a patch
  dx = diff(p.x(2:3)); % x space step
  dy = 2/ny;           % y space step
  ct = nan+c;          % preallocate output array
  pcs = reshape(p.cs,nx-1,[],2);
%{
\end{matlab}
Compute the cross-channel flux using `ghost' nodes at
channel boundaries, so that the flux is zero at $y=\pm1$
either because the boundary values are replicated so the
differences are zero, or because the diffusivities in
\verb|cs| are zero at the channel boundaries.
\begin{matlab}
%}
  ydif = pcs(ix,1:2:end,2) ...
         .*(c(ix,[1:end end],:,:)-c(ix,[1 1:end],:,:))/dy;
%{
\end{matlab}
Now evaluate advection-diffusion time
derivative~\eqref{eq:ddeChanDisp}. Could use upwind
advection and no longitudinal diffusion, or, as here,
centred advection and diffusion.
\begin{matlab}
%}
  ct(ix,:,:,:) = (ydif(:,2:end,:,:)-ydif(:,1:end-1,:,:))/dy ...
  + diff(pcs(:,2:2:end,2).*diff(c))/dx^2 ...
  - p.Pe*pcs(ix,2:2:end,1).*(c(ix+1,:,:,:)-c(ix-1,:,:,:))/(2*dx); 
end% function
%{
\end{matlab}
%}
