% Codes a nonlinear water wave PDE on a staggered 1D grid
% inside patches in space.  Used by waterWaveExample.m
% AJR, 4 Apr 2019 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{waterWavePDE()}: water wave PDE}
\label{sec:waterWavePDE}

This function codes the staggered lattice equation inside
the patches for the nonlinear wave-like \pde\
system~\cref{eqs:patch:N}. Also, regularise the absolute
value appearing the the \pde{}s via the one-line
function~\verb|rabs()|.
\begin{matlab}
%}
function Ut = waterWavePDE(t,U,patches)
  rabs = @(u) sqrt(1e-4 + u.^2);
%{
\end{matlab}
As before, set the micro-grid spacing, reserve space for
time derivatives, and index the patch-interior points of the
micro-grid.
\begin{matlab}
%}
  dx = diff(patches.x(2:3));
  U = squeeze(U);
  Ut = nan(size(U));  ht = Ut;
  i = 2:size(U,1)-1;
%{
\end{matlab}
Need to estimate~\(h\) at all the \(u\)-points, so 
into~\verb|V| use averages, and linear extrapolation to
patch-edges.
\begin{matlab}
%}
  ii = i(2:end-1);
  V = Ut;
  V(ii,:) = (U(ii+1,:)+U(ii-1,:))/2;
  V(1:2,:) = 2*U(2:3,:)-V(3:4,:);
  V(end-1:end,:) = 2*U(end-2:end-1,:)-V(end-3:end-2,:);
%{
\end{matlab}
Then estimate \(\D x{(hu)}\) from~\(u\) and the
interpolated~\(h\) at the neighbouring micro-grid points.
\begin{matlab}
%}
  ht(i,:) = -(U(i+1,:).*V(i+1,:)-U(i-1,:).*V(i+1,:))/(2*dx);
%{
\end{matlab}
Correspondingly estimate the terms in the momentum \pde:
\(u\)-values in~\(\verb|U|_i\) and~\(\verb|V|_{i\pm1}\); and
\(h\)-values in~\(\verb|V|_i\) and~\(\verb|U|_{i\pm1}\).
\begin{matlab}
%}
  Ut(i,:) = -0.985*(U(i+1,:)-U(i-1,:))/(2*dx) ...
    -0.003*U(i,:).*rabs(U(i,:)./V(i,:)) ...
    -1.045*U(i,:).*(V(i+1,:)-V(i-1,:))/(2*dx) ...
    +0.26*rabs(V(i,:).*U(i,:)).*(V(i+1,:)-2*U(i,:)+V(i-1,:))/dx^2/2;
%{
\end{matlab}
where the mysterious division by two in the second
derivative is due to using the averaged values of~\(u\) in
the estimate:
\begin{eqnarray*}
u_{xx}&\approx&\frac1{4\delta^2}(u_{i-2}-2u_i+u_{i+2})
\\&=&\frac1{4\delta^2}(u_{i-2}+u_i-4u_i+u_i+u_{i+2})
\\&=&\frac1{2\delta^2}\left(\frac{u_{i-2}+u_i}2-2u_i+\frac{u_i+u_{i+2}}2\right)
\\&=&\frac1{2\delta^2}\left(\bar u_{i-1}-2u_i+\bar u_{i+1}\right).
\end{eqnarray*}
Then overwrite the unwanted~\(\dot u_{ij}\) with the
corresponding wanted~\(\dot h_{ij}\).
\begin{matlab}
%}
  Ut(patches.hPts) = ht(patches.hPts);
end
%{
\end{matlab}
%}
