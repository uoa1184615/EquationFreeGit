% Simulate heterogeneous diffusion in 2D on patches as an
% example application of patches in space. Here compare to
% high contrast cases and discussion by Elise Fressart and
% Barbara Verfurth 2303.15151
% AJR, 23 Oct 2023
%{
\documentclass[10pt,a5paper]{article}
% fancyvrb does code listing, including line numbers
\usepackage{fancyvrb}
\newenvironment{matlab}%
    {\Verbatim[numbers=left,firstnumber=\the\inputlineno]}%
    {\endVerbatim}
% Also get fancyvrb to omit %{ and %} pairs, 
% although this requires they always be used
\makeatletter
  \edef\FancyVerbStartString{\@percentchar\@charrb} 
  \edef\FancyVerbStopString{\@percentchar\@charlb} 
\makeatother
\usepackage{ajr,mycleveref}
% These are recommended by Rob J Hyndman (2011)
% \footnote{\url{http://robjhyndman.com/researchtips/latex-floats/}}
\setcounter{topnumber}{2}
\setcounter{bottomnumber}{2}
\setcounter{totalnumber}{4}
\renewcommand{\topfraction}{0.85}
\renewcommand{\bottomfraction}{0.85}
\renewcommand{\textfraction}{0.15}
\renewcommand{\floatpagefraction}{0.7}

\title{Computational homogenisation of a 2D diffusion or
waves with high contrast inclusion}
\author{A.J. Roberts}
\date{\today}

\begin{document}

\maketitle

The issues raised by Elise Fressart and Barbara Verfurth
2303.15151 are the homogenisation-like modelling of wave
propagation through material with microscale high-contrast
`elasticity'. Here I address the analogous problem of
diffusion in high-contrast material as the issues are much
the same, the discussion is more rigorous, and cognate
results for waves are deduced by taking the square-root of
the eigenvalues in order to get frequencies of the waves.

One can change the parameters of this code, but as is it
solve diffusion in a 2D domain of~\([0,1]^2\) with
macroscale periodic boundary conditions. The heterogeneity
period~\(\epsilon\) of each cell is set to~\(1/9\). Within
each cell, it resolves the sub-period structure on a
\(16\times16\) microscale grid of spacing \(dx=0.0069\)
which is perfectly adequate for demonstrating the typical
behaviour. As done by Fressart and Verfurth, The
heterogeneity is that the diffusion (elasticity) coefficient
is one outside the centre square of each cell, and~\(a_0\)
inside the centre square---they reported the cases \(a_0 \in
\{1/2, 1/2^5, 1/2^{10}\}\).

To clearly differentiate, where possible, the difference
between the desired macroscale homogenisation and the
microscale sub-cell dynamics, I here invoke the patch
scheme, see \cref{figSim7}. \begin{figure}\centering
\caption{\label{figSim7}example of patch scheme with nine
patches on the unit square simulated to time~\(0.05\).  With
inclusions having low diffusivity \(a_0=1/2^7\), the
diffusion into and out of the inclusions takes quite a long
time.  For waves, the waves within each inclusion would
bounce around inside the inclusion and only slowly
leak/radiate outside.}
\includegraphics[width=\textwidth]{hiContrastSim7}
\end{figure}%
I choose to only resolve macroscale modes with wavelengths
longer than~\(0.5\) by choosing \(3\times3\) patches in the
domain. Each patch is one cell, hence of side
length~\(1/9\). Where discussed, the sub-patch dynamics are
essentially the same as the sub-cell dynamics. The patches
are coupled by spectral interpolation to ensure high
accuracy for any macroscale modes---whatever the `macroscale
modes' might be.


\paragraph{Example of homogeneous diffusion/waves} The
example of diffusion in homogeneous material, \(a_0=1\),
illustrates the distinction that the patch scheme makes
between macroscale and sub-cell modes.
\begin{figure}\centering
\caption{\label{figHomogeneous}eigenvalues of the
homogeneous diffusion (\(a_0=1\)).  It shows a spectral gap
from \(-79\) to \(-2501\) separating the macroscale of
interest from the sub-cell microscale.}
\includegraphics[width=\textwidth]{hiContrastHomogeneous}
\end{figure}%
We easily characterise the dynamics of the problem by
exploring the eigenvalues. \cref{figHomogeneous} plots all
the eigenvalues, as a function of their index on
quasi-log-log axes.
\begin{itemize}
\item The sole eigenvalue \(\lambda=0\) represents
conservation of stuff.
\item The next group of four at \(\lambda=-39\) represent
the four macroscale modes/waves with wavenumber
\((0,\pm2\pi)\) or~\((\pm2\pi,0)\).
\item The next group of four at \(\lambda=-79\) represent
the four macroscale modes/waves with wavenumber
\((\pm2\pi,\pm2\pi)\).
\item The remaining eigenvalues \(\lambda<-2500\) represent
high-wavenumber, microscale, sub-cell modes separated by a
spectral gap of ratio~\(31\).  
\end{itemize}
In wave problems, the spectral gap would be between slow,
macroscale waves of frequencies\({}<9\), and fast,
microscale, sub-cell, waves of frequencies\({}>50\).


\paragraph{High contrast erodes the spectral gap} Decreasing
the diffusivity inside the inclusion down to \(a_0=2^{-10}\)
changes the problem to one of high-contrast.
\cref{figHetero} plots the leading eigenvalues as a function
of~\(a_0\). \begin{figure}\centering
\caption{\label{figHetero}leading 19 eigenvalues of
heterogeneous diffusion as function of the inclusion's
diffusivity~\(a_0\).  The spectral gap for
large~\(a_0\approx 0.1\) closes as \(a_0\)~decreases
through~\(0.01\) as sub-inclusion modes become
long-lasting---equivalently as sub-inclusion waves become
slow.}
\includegraphics[width=\textwidth]{hiConstrastHetero19}
\end{figure}%
For such very small~\(a_0\), all the sub-inclusion modes
decay slowly so the become long-lasting modes, and so should
be considered as part of the `homogenised'
modelling---unless one can guarantee from initial conditions
(or otherwise) that they do not arise. For waves, for such
very small~\(a_0\), all the sub-inclusion modes become low
frequency waves and so similarly should be considered as
part of the `homogenised' modelling---unless the initial
conditions (or otherwise) ensure that they do not arise.



\section{\texttt{hiContrastDiff2}: computational
homogenisation of a 2D diffusion with high contrast
inclusion}
\label{sec:hiContrastDiff2}



First set heterogeneous diffusivities constant in each of
inclusion and exterior. 
\begin{matlab}
%}
a0 = 1/2^7
mPeriod = 16
xi=(0.5:mPeriod)/mPeriod;
incl = (abs(xi'-1/2)<1/4)&(abs(xi-1/2)<1/4);
cHetr = incl*a0+(~incl)*1;
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  Use macroscale periodic and
spectral interpolation.  In 2D we get only real eigenvalues
by using edgy interpolation.  
\begin{matlab}
%}
edgyInt = true; 
nPatch = 3
nSubP = mPeriod+2  
dx = 1/(mPeriod*nPatch) % this is for full domain
dx = dx/3 % use smaller periodicity separated by gaps
configPatches2(@heteroDiff2,[0 1],'periodic',nPatch ...
    ,0,dx,nSubP ,'EdgyInt',edgyInt ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Simulate}
Set initial conditions of a simulation (although what is FVs v0?).
\begin{matlab}
%}
global patches
sigma = 0.1
u0 = exp( -(patches.x-0.5).^2/sigma^2-(patches.y-0.5).^2/sigma^2 ); 
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
    [ts,us] = ode23(@patchSys2, 0.05*linspace(0,1).^2, u0(:));
%{
\end{matlab}

\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
disp('plot animation of solution field')
figure(1), clf, colormap(flipud(parula))
%{
\end{matlab}
Get spatial coordinates and pad them with NaNs to separate
patches.
\begin{matlab}
%}
x = squeeze(patches.x); y = squeeze(patches.y);
x(end+1,:)=nan;  y(end+1,:)=nan; % pad with nans
%{
\end{matlab}
For every time step draw the surface and pause for a short
display.
\begin{matlab}
%}
for i = 1:length(ts)
%{
\end{matlab}
Get the row vector of data,  form into the 6D array via the
interpolation to the edges, then pad with Nans between
patches, and reshape to suit the surf function.
\begin{matlab}
%}
  u = squeeze( mean( patchEdgeInt2(us(i,:)) ,4));
  u(end+1,:,:,:)=nan; u(:,end+1,:,:)=nan;
  u = reshape(permute(u,[1 3 2 4]), [numel(x) numel(y)]);
%{
\end{matlab}
If the initial time then draw the surface with labels,
otherwise just update the surface data.
\begin{matlab}
%}
  if i==1
       hsurf = surf(x(:),y(:),u'); view(60,40) 
       xlim([0 1]), ylim([0 1]),   caxis([0 1])
       xlabel('$x$'), ylabel('$y$'), zlabel('$u(x,y)$')
  else set(hsurf,'ZData', u');
  end
  legend(['time = ' num2str(ts(i),2)],'Location','north')
  pause(0.05)
%{
\end{matlab}
finish the animation loop and if-plot.
\begin{matlab}
%}
end%for over time
%{
\end{matlab}






\subsection{Compute Jacobian and its spectrum}
Let's explore the dynamics via the Jacobian. Find which
elements of the 6D array are interior micro-grid points and
hence correspond to dynamical variables.
\begin{matlab}
%}
    u0 = zeros([nSubP,nSubP,1,1,nPatch,nPatch]);
    u0([1 end],:,:) = nan;
    u0(:,[1 end],:) = nan;
    i = find(~isnan(u0));
%{
\end{matlab}
Construct the Jacobian of the scheme as the matrix of the
linear transformation, obtained by transforming the standard
unit vectors.
\begin{matlab}
%}
    Jac = nan(length(i));
    sizeJacobian = size(Jac)
    for j = 1:length(i)
      u = u0(:)+(i(j)==(1:numel(u0))');
      tmp = patchSys2(0,u);
      Jac(:,j) = tmp(i);
    end
%{
\end{matlab}
Test for symmetry, with error if we know it should be
symmetric.
\begin{matlab}
%}
    notSymmetric=norm(Jac-Jac')    
    assert(notSymmetric<1e-7,'failed symmetry')
%{
\end{matlab}
Find all the eigenvalues (as \verb|eigs| is unreliable).
\begin{matlab}
%}
    [evecs,evals] = eig((Jac+Jac')/2,'vector');
%{
\end{matlab}
Sort eigenvalues on their real-part with most positive
first, and most negative last.  List leading, and plot all.
\begin{matlab}
%}
    [~,k] = sort(-real(evals));
    evals=evals(k); evecs=evecs(:,k);
    leadingEvals=evals(1:2*nPatch^2+1)
    figure(2),clf
    plot(evals,'.')
    xlabel('index'),ylabel('eigenvalue $\lambda$')
    quasiLogAxes(1,10)
%{
\end{matlab}




\subsection{\texttt{heteroDiff2()}: heterogeneous diffusion}
\label{sec:heteroDiff2}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 6D input arrays~\verb|u|, \verb|x|,
and~\verb|y| (via edge-value interpolation of
\verb|patchSys2|), computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.  The
two 2D array of diffusivities,~$c^x_{ij}$ and~$c^y_{ij}$,
have previously been stored in~\verb|patches.cs| (3D). 
\begin{matlab}
%}
function ut = heteroDiff2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;        % preallocate output array
  ut(i,j,:,:,:,:) ...
  = diff(patches.cs(:,j,:).*diff(u(:,j,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(i,:,:).*diff(u(i,:,:,:,:,:),1,2),1,2)/dy^2; 
end% function
%{
\end{matlab}
\end{document}
%}