% Solve for steady state of twoscale heterogeneous diffusion
% in 2D on patches as an example application, varied from
% example of section 5.1 of Abdulle et al., (2020b).  AJR,
% Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{abdulleDiffEquil2}: equilibrium of a 2D
twoscale heterogeneous diffusion via small patches}
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
+\frac{2+\sin\pi y/\epsilon}{2+1.8\cos2\pi x/\epsilon}.
\end{equation*}
The solution shows some nice little microscale wiggles.

Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches i
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
%{
\end{matlab}
Set the periodicity, via~\(\epsilon\), and other microscale
parameters.
\begin{matlab}
%}
nPeriodsPatch = 1 % any integer
epsilon = 2^(-4) % not tiny, so we can see patches
dx = epsilon/mPeriod
nSubP = nPeriodsPatch*mPeriod+2 % when edgy int
%{
\end{matlab}



\paragraph{Patch configuration}
Choose either Dirichlet (default) or Neumann on the left
boundary in coordination with micro-code in
\cref{sec:abdulleDiffForce2}
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
Dom.type='equispace';
configPatches2(@abdulleDiffForce2,[0 1],Dom ...
    ,nPatch ,4 ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',a );
%{
\end{matlab}




\paragraph{Solve for steady state}
Set initial guess of zero, with \verb|NaN| to indicate
patch-edge values.  Index~\verb|i| are the indices of
patch-interior points, and the number of unknowns is then
its length.
\begin{matlab}
%}
u0 = zeros(nSubP,nSubP,1,1,nPatch,nPatch);
u0([1 end],:,:) = nan;  u0(:,[1 end],:) = nan;
i = find(~isnan(u0));
nVariables = numel(i)
%{
\end{matlab}
Solve by iteration.  Use \verb|fsolve| for simplicity and
robustness (and using \verb|optimoptions| to omit trace
information).
\begin{matlab}
%}
tic;
uSoln = fsolve(@theRes,u0(i) ...
        ,optimoptions('fsolve','Display','off')); 
solnTime = toc
%{
\end{matlab}
Store the solution into the patches, and give magnitudes.
\begin{matlab}
%}
u0(i) = uSoln;
normSoln = norm(uSoln)
normResidual = norm(theRes(uSoln))
%{
\end{matlab}




\paragraph{Draw solution profile}
First reshape arrays to suit 2D space surface plots.
\begin{matlab}
%}
figure(1), clf, colormap(hsv)
x = squeeze(patches.x); y = squeeze(patches.y);
u = reshape(permute(squeeze(u0),[1 3 2 4]), [numel(x) numel(y)]);
%{
\end{matlab}
Draw the patch solution surface, with edge-values omitted as
already~\verb|NaN| by not bothering to interpolate them.
\begin{matlab}
%}
surf(x(:),y(:),u'); view(60,55) 
xlabel('x'), ylabel('y'), zlabel('u(x,y)')
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
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
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
  ut(ix,iy,:) = diff(patches.cs(:,iy).*diff(u(:,iy,:)))/dx^2 ...
      + diff(patches.cs(ix,:).*diff(u(ix,:,:),1,2),1,2)/dy^2 ...
      + 10; 
end%function abdulleDiffForce2
%{
\end{matlab}


\subsection{\texttt{theRes()}: function to zero}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system, and returns the vector of patch-interior time
derivatives.
\begin{matlab}
%}
function f=theRes(u)
  global patches i
  v=nan(size(patches.x+patches.y));
  v(i)=u;
  f=patchSys2(0,v(:),patches);
  f=f(i);
end%function theRes
%{
\end{matlab}

Fin.
%}
