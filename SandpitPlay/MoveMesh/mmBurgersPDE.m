% mmBurgersPDE() is a microscale discretisation of Burgers'
% PDE on multiscale patches of a lattice x, when the patches
% move according to some scheme.  AJR, Aug 2021
%!TEX root = doc.tex
%{
\subsection{\texttt{mmBurgersPDE()}: Burgers PDE inside a
moving mesh of patches} 
For the evolving scalar field~\(u(t,x)\), we code a
microscale discretisation of Burgers' \pde\ \(u_t = \epsilon
u_{xx} -uu_x\), for say \(\epsilon=0.02\)\,, when the
patches of microscale lattice move with various
velocities~\(V\).  
\begin{matlab}
%}
function ut = mmBurgersPDE(t,u,V,D,patches)
epsilon = 0.02;
%{
\end{matlab}

\paragraph{Generic input/output variables}
\begin{itemize}
\item \verb|t| (scalar) current time---not used here as the
\pde\ has no explicit time dependence (autonomous).
\item \verb|u| (\(n\times1\times1\times N\)) field values on
the patches of microscale lattice.
\item \verb|V| (\(1\times1\times1\times N\)) moving velocity of
the \(j\)th~patch.
\item \verb|D| (\(1\times1\times1\times N\)) displacement of
the \(j\)th~patch from the fixed spatial positions stored in
\verb|patches.x|---not used here as the \pde\ has no
explicit space dependence (homogeneous).
\item \verb|patches| struct of patch configuration
information.
\item \verb|ut| (\(n\times1\times1\times N\)) output
computed values of the time derivatives \(Du/Dt\) on the
patches of microscale lattice.
\end{itemize}

Here there is only one field variable, and one in the
ensemble, so for simpler coding of the \pde\ we squeeze them
out (no need to reshape when via \verb|mmPatchSmooth1|).
\begin{matlab}
%}
  u=squeeze(u);     % omit singleton dimensions
  V=shiftdim(V,2);  % omit two singleton dimens
%{
\end{matlab}

\paragraph{Burgers PDE}
In terms of the moving derivative \(Du/Dt:=u_t+Vu_x\) the
\pde\ becomes \(Du/Dt=\epsilon u_{xx}+(V-u)u_x\)\,. So code
for every patch that \(\dot u_{ij} =\frac\epsilon{h^2}
(u_{i+1,j}-2u_{i,j}+u_{i-1,j}) +(V_j-u_{ij})
\frac1{2h}(u_{i+1,j}-u_{i-1,j})\) at all interior lattice
points.
\begin{matlab}
%}
  dx=diff(patches.x(1:2)); % microscale spacing
  i=2:size(u,1)-1;  % interior points in patches
  ut=nan+u;         % preallocate output array
  ut(i,:) = epsilon*diff(u,2)/dx^2 ...
      +(V-u(i,:)).*(u(i+1,:)-u(i-1,:))/(2*dx);
end
%{
\end{matlab}
%}
