% Solve for steady state of multiscale heterogeneous diffusion
% in 2D on patches as an example application, varied from
% example of section 5.1 of Abdulle et al., (2020b).  AJR,
% 31 Jan 2023
%!TEX root = doc.tex
%{
\section{\texttt{abdulleDiffEquil2}: equilibrium of a 2D
multiscale heterogeneous diffusion via small patches}
\label{sec:abdulleDiffEquil2}

Here we find the steady state~\(u(x,y)\) to the
heterogeneous \pde\ \cite[inspired by][\S5.1]{Abdulle2020b}
\begin{equation*}
u_t=\divv[a(x,y)\grad u]+10,
\end{equation*}
on square domain \([0,1]^2\) with zero-Dirichlet BCs, for
coefficient `diffusion' matrix, varying with
period~\(\epsilon\) of (their~(45))
\begin{equation*}
a:=\frac{2+1.8\sin2\pi x/\epsilon}{2+1.8\cos2\pi y/\epsilon}
+\frac{2+\sin2\pi y/\epsilon}{2+1.8\cos2\pi x/\epsilon}.
\end{equation*}
\cref{fig:abdulleDiffEquil2} shows solutions have some 
nice microscale wiggles reflecting the heterogeneity.
\begin{figure}
\centering\caption{\label{fig:abdulleDiffEquil2}%
Equilibrium of the macroscale diffusion problem of Abdulle
with boundary conditions of Dirichlet zero-value except for
\(x=0\) which is Neumann (\cref{sec:abdulleDiffEquil2}). 
Here the patches have a Chebyshev-like spatial distribution.
The patch size is chosen large enough to see within.}
\includegraphics[scale=0.8]{Figs/abdulleDiffEquil2}
\end{figure}

Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches
%global OurCf2eps, OurCf2eps=true %option to save plot
%{
\end{matlab}


First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the spatial micro-grid
lattice.  Then \verb|configPatches2| replicates the
heterogeneity to fill each patch.  (These diffusion
coefficients should really recognise the half-grid-point
shifts, but let's not bother.)
\begin{matlab}
%}
mPeriod = 6
x = (0.5:mPeriod)'/mPeriod;  y=x';
a = (2+1.8*sin(2*pi*x))./(2+1.8*sin(2*pi*y)) ...
   +(2+    sin(2*pi*y))./(2+1.8*sin(2*pi*x));
diffusivityRange = [min(a(:)) max(a(:))]
%{
\end{matlab}
Set the periodicity~\(\epsilon\), here big enough so we can
see the patches, and other microscale parameters.
\begin{matlab}
%}
epsilon = 0.04 
dx = epsilon/mPeriod
nPeriodsPatch = 1 % any integer
nSubP = nPeriodsPatch*mPeriod+2 % when edgy int
%{
\end{matlab}



\paragraph{Patch configuration} Choose either Dirichlet
(default) or Neumann on the left boundary in coordination
with micro-code in \cref{sec:abdulleDiffForce2}
\begin{matlab}
%}
Dom.bcOffset = zeros(2);
if 1, Dom.bcOffset(1)=0.5; end% left Neumann
%{
\end{matlab}
Say use \(7\times7\) patches in \((0,1)^2\), fourth order
interpolation, and either `equispace' or `chebyshev':
\begin{matlab}
%}
nPatch = 7 
Dom.type='chebyshev';
configPatches2(@abdulleDiffForce2,[0 1],Dom ...
    ,nPatch ,4 ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',a );
%{
\end{matlab}




\paragraph{Solve for steady state} Set initial guess of
zero, with \verb|NaN| to indicate patch-edge values.
Index~\verb|i| are the indices of patch-interior points, and
the number of unknowns is then its length.
\begin{matlab}
%}
u0 = zeros(nSubP,nSubP,1,1,nPatch,nPatch);
u0([1 end],:,:) = nan;  u0(:,[1 end],:) = nan;
patches.i = find(~isnan(u0));
nVariables = numel(patches.i)
%{
\end{matlab}
Solve by iteration.  Use \verb|fsolve| for simplicity and
robustness (and using \verb|optimoptions| to omit trace
information), and give magnitudes.
\begin{matlab}
%}
tic;
uSoln = fsolve(@theRes2,u0(patches.i) ...
        ,optimoptions('fsolve','Display','off')); 
solnTime = toc
normResidual = norm(theRes2(uSoln))
normSoln = norm(uSoln)
%{
\end{matlab}
Store the solution vector into the patches, and interpolate,
but have not bothered to set boundary values so they stay
NaN from the interpolation.
\begin{matlab}
%}
u0(patches.i) = uSoln;
u0 = patchEdgeInt2(u0);
%{
\end{matlab}




\paragraph{Draw solution profile} Separate patches with
NaNs, then reshape arrays to suit 2D space surface plots.
\begin{matlab}
%}
figure(1), clf, colormap(0.8*hsv)
patches.x(end+1,:,:)=nan;  u0(end+1,:,:)=nan;  
patches.y(:,end+1,:)=nan;  u0(:,end+1,:)=nan;
u = reshape(permute(squeeze(u0),[1 3 2 4]) ...
    , [numel(patches.x) numel(patches.y)]);
%{
\end{matlab}
Draw the patch solution surface, with boundary-values
omitted as already~\verb|NaN| by not bothering to set them.
\begin{matlab}
%}
mesh(patches.x(:),patches.y(:),u'); view(60,55) 
xlabel('space $x$'), ylabel('space $y$'), zlabel('$u(x,y)$')
ifOurCf2eps(mfilename) %optionally save plot
%{
\end{matlab}






\subsection{\texttt{abdulleDiffForce2()}: microscale
discretisation inside patches of forced diffusion PDE}
\label{sec:abdulleDiffForce2}

This function codes the lattice heterogeneous diffusion of
the \pde\ inside the patches.  For 6D input arrays~\verb|u|,
\verb|x|, and~\verb|y|, computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.   
\begin{matlab}
%}
function ut = abdulleDiffForce2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
Set Dirichlet boundary value of zero around the square
domain, but also cater for zero Neumann condition on the
left boundary.
\begin{matlab}
%}
  u( 1 ,:,:,:, 1 ,:)=0; % left edge of left patches
  u(end,:,:,:,end,:)=0; % right edge of right patches
  u(:, 1 ,:,:,:, 1 )=0; % bottom edge of bottom patches
  u(:,end,:,:,:,end)=0; % top edge of top patches
  if 1, u(1,:,:,:,1,:)=u(2,:,:,:,1,:); end% left Neumann
%{
\end{matlab}
Compute the time derivatives via stored forcing and
coefficients.  Easier to code by conflating the last four
dimensions into the one~\verb|,:|.
\begin{matlab}
%}
  ut(i,j,:) = diff(patches.cs(:,j).*diff(u(:,j,:)))/dx^2 ...
      + diff(patches.cs(i,:).*diff(u(i,:,:),1,2),1,2)/dy^2 ...
      + 10; 
end%function abdulleDiffForce2
%{
\end{matlab}


\subsection{\texttt{theRes2()}: function to zero}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system, and returns the vector of patch-interior time
derivatives.
\begin{matlab}
%}
function f=theRes2(u)
  global patches i
  v=nan(size(patches.x+patches.y));
  v(patches.i)=u;
  f=patchSys2(0,v(:),patches);
  f=f(patches.i);
end%function theRes2
%{
\end{matlab}
%}
