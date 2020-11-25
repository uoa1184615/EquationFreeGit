% Computes the time derivatives of heterogeneous wave
% in 1D on patches.  Used by homoWaveEdgy1.m,
% AJR, 26 Nov 2019 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroWave()}: wave in heterogeneous
media with weak viscous damping}
\label{sec:heteroWave}

This function codes the lattice heterogeneous wave equation,
with weak viscosity, inside the patches.  For 3D input
array~\verb|u| (\(u_{ij} = \verb|u(i,1,j)|\) and \(v_{ij} =
\verb|u(i,2,j)|\)) and 2D array~\verb|x| (obtained in full
via edge-value interpolation of \verb|patchSmooth1|,
\cref{sec:patchSmooth1}), computes the time derivatives at
each point in the interior of a patch, output in~\verb|ut|:
\begin{equation*}
\D t{u_{ij}}=v_{ij}\,,\quad
\D t{v_{ij}}= \frac1{dx^2}\delta[c_{i-1/2}\delta u_{ij}]
+\frac{0.02}{dx^2}\delta^2 v_{ij}\,.
\end{equation*}
The column vector (or possibly array) of diffusion
coefficients~\(c_i\) have previously been stored in
struct~\verb|patches|.
\begin{matlab}
%}
function ut = heteroWave(t,u,patches)
  u = squeeze(u);
  dx = diff(patches.x(2:3));    % space step
  i = 2:size(u,1)-1;    % interior points in a patch
  ut = nan(size(u));    % preallocate output array
  ut(i,1,:) = u(i,2,:); % du/dt=v then dvdt=
  ut(i,2,:) = diff(patches.cs.*diff(u(:,1,:)))/dx^2 ...
        +0.02*diff(u(:,2,:),2)/dx^2; 
end% function
%{
\end{matlab}
%}
