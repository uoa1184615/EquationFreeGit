% Codes the ideal wave PDE on a staggered 1D grid inside
% patches in space.  Used by waterWaveExample.m
% AJR, 4 Apr 2019 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{idealWavePDE()}: ideal wave PDE}
\label{sec:idealWavePDE}

This function codes the staggered lattice equation inside
the patches for the ideal wave \pde\ system \(h_t=-u_x\) and
\(u_t=-h_x\). Here code for a staggered micro-grid,
index~\(i\), of staggered macroscale patches, index~\(j\):
the array
\begin{equation*}
U_{ij}=\begin{cases} u_{ij}&i+j\text{ even},\\
h_{ij}& i+j\text{ odd}.
\end{cases}
\end{equation*}
The output~\verb|Ut| contains the merged time derivatives of
the two staggered fields. So set the micro-grid spacing and
reserve space for time derivatives.
\begin{matlab}
%}
function Ut = idealWavePDE(t,U,patches)
  dx = diff(patches.x(2:3));
  U = squeeze(U);
  Ut = nan(size(U));  ht = Ut;
%{
\end{matlab}
Compute the \pde\ derivatives only at interior micro-grid
points of the patches.
\begin{matlab}
%}
  i = 2:size(U,1)-1;
%{
\end{matlab}
Here `wastefully' compute time derivatives for both \pde{}s
at all grid points---for simplicity---and then merge the
staggered results. Since \(\dot h_{ij} \approx -(u_{i+1,j}
-u_{i-1,j}) /(2\cdot dx) =-(U_{i+1,j} -U_{i-1,j}) /(2\cdot
dx)\) as adding\slash subtracting one from the index of a
\(h\)-value is the location of the neighbouring \(u\)-value
on the staggered micro-grid.
\begin{matlab}
%}
  ht(i,:) = -(U(i+1,:)-U(i-1,:))/(2*dx);
%{
\end{matlab}
Since \(\dot u_{ij} \approx -(h_{i+1,j} -h_{i-1,j}) /(2\cdot
dx) =-(U_{i+1,j} -U_{i-1,j}) /(2\cdot dx)\) as adding\slash
subtracting one from the index of a \(u\)-value is the
location of the neighbouring \(h\)-value on the staggered
micro-grid.
\begin{matlab}
%}
  Ut(i,:) = -(U(i+1,:)-U(i-1,:))/(2*dx);
%{
\end{matlab}
Then overwrite the unwanted~\(\dot u_{ij}\) with the
corresponding wanted~\(\dot h_{ij}\).
\begin{matlab}
%}
  Ut(patches.hPts) = ht(patches.hPts);
end
%{
\end{matlab}
%}
