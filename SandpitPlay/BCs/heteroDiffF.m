% Computes the time derivatives of forced heterogeneous
% diffusion in 1D on patches.  AJR, Apr 2019 -- 3 Jan 2023
%!TEX root = doc.tex
%{
\subsection{\texttt{heteroDiffF()}: forced heterogeneous diffusion}
\label{sec:heteroDiffF}

This function codes the lattice heterogeneous diffusion
inside the patches with forcing and with microscale boundary
conditions on the macroscale boundaries.  Computes the time
derivative at each point in the interior of a patch, output
in~\verb|ut|.  The column vector of diffusivities~\(a_i\)
has been stored in struct~\verb|patches.cs|, as has the
array of forcing coefficients.
\begin{matlab}
%}
function ut = heteroDiffF(t,u,patches)
  global microTimePeriod
  % macroscale Dirichlet BCs
  u( 1 ,:,:, 1 )=0; % left-edge of leftmost is zero
  u(end,:,:,end)=0; % right-edge of rightmost is zero
  % interior forced diffusion
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan+u;          % preallocate output array
  if microTimePeriod>0 % optional time fluctuations
     at = cos(2*pi*t/microTimePeriod)/30; 
  else at=0; end
  ut(i,:,:,:) = diff((patches.cs(:,1,:)+at).*diff(u))/dx^2 ...
      +patches.f2(i,:,:,:)*t^2+patches.f1(i,:,:,:)*t; 
end% function
%{
\end{matlab}
%}