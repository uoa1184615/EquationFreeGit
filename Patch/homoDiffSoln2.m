% Solve for steady state of heterogeneous diffusion in 2D on
% patches as an example application. The microscale is of
% known period so we interpolate next-to-edge values to get
% opposite edge values.  This version implements scenarios
% inspired by Biezemans et al. (2022) \S3, p.12.  AJR, Apr
% 2022
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffSoln2}: steady state of a 2D
heterogeneous diffusion via small patches}
\label{sec:homoDiffSoln2}

Here we find the steady state~\(u(x,y)\) to the heterogeneous \pde
\begin{equation*}
u_t=\divv[c(x,y)\grad u]-u+f,
\quad\text{for } f=100\sin(\pi x)\sin(\pi y).
\end{equation*}
The heterogeneous diffusion~\(c\) varies over two orders of
magnitude in small space distance~\(\epsilon\).  I
include~\(-u\) in the \pde\ to ensure a steady state with
periodic BCs.

\cite{Biezemans2022} discussed an example homogenisation in
2D with heterogeneity of period~\(\epsilon:=\pi/150\) in
both directions.  Ensure integer multiple of heterogeneity
periods in the domain, and initially use three times
bigger~\(\epsilon\).
\begin{matlab}
%}
epsilon = 1/round(50/pi)
%{
\end{matlab}
\cite{Biezemans2022} choose microscale mesh spacing
of~$1/1024$, so the number of micro-grid points in one
period would be~\(1024\epsilon\). But \emph{initially} use
less.
\begin{matlab}
%}
mPeriod = round(128*epsilon) %round(1024*epsilon)
%{
\end{matlab}
So the migro-grid spacing is exactly
\begin{matlab}
%}
dx = epsilon/mPeriod   
%{
\end{matlab}

\paragraph{Diffusivities}
Now form one period of the heterogeneity diffusivities.
\cite{Biezemans2022} used \(c=1 +100 \cos^2(\pi x/\epsilon)
\sin^2(\pi y/\epsilon) \).  Need to shift phases of the
diffusivity by half-micro-grid for diffusivities in each
direction to form two diffusivity matrices on the microscale
lattice.
\begin{matlab}
%}
cHetr=[];
v=(  1:mPeriod)/mPeriod;
h=(0.5:mPeriod)/mPeriod;
cHetr(:,:,1) = 1+100*cos(pi*h').^2*sin(pi*v).^2;
cHetr(:,:,2) = 1+100*cos(pi*v').^2*sin(pi*h).^2;
%{
\end{matlab}
Plot surfaces of the diffusivity.
\begin{matlab}
%}
figure(2),surf(h,v,cHetr(:,:,2))
hold on,  surf(v,h,cHetr(:,:,1))
hold off, alpha 0.5, drawnow
%{
\end{matlab}

\paragraph{Patch configuration}
As is common, \cite{Biezemans2022} implemented
zero-Dirichlet BCs on $(0,1)^2$.  Here these are
more-or-less encompassed by implementing periodic BCs on
$(-1,1)^2$. Initially use \(8\times8\) patches to have
\(4\times4\) patches in \((0,1)^2\), which then have patch
spacing~\(H\).
\begin{matlab}
%}
nPatch = [8 8]
H = 2./nPatch 
HepsilonRatio = H/epsilon  
%{
\end{matlab}
Best when each patch spans an integral number of periods
plus two grid steps.  The smallest patches are 
\begin{matlab}
%}
nSubP = [1 1]*mPeriod+2 
%{
\end{matlab}
Consequently, the ratio of space computed on, to the space
in the domain is the product of the following ratios in each
direction, namely about~8\% here.
\begin{matlab}
%}
ratio = ((nSubP-1)*dx)./H
%{
\end{matlab}
Specify spectral interpolation.  The edgy interpolation is
self-adjoint \cite[]{Bunder2020a} leading to a symmetric
matrix problem.
\begin{matlab}
%}
configPatches2(@hetDiffForce2,[-1 1 -1 1],nan,nPatch ...
    ,0,ratio,nSubP ,'EdgyInt',true  ...
    ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Solve for steady state}
Set initial guess of zero, with \verb|NaN| to indicate
patch-edge values.  Index~\verb|i| are the indices of
patch-interior points, and the number of unknowns is then
its length.
\begin{matlab}
%}
global patches i
u0 = zeros([nSubP,1,1,nPatch]);
u0([1 end],:,:) = nan; u0(:,[1 end],:) = nan;
i = find(~isnan(u0));
nVars = numel(i)
%{
\end{matlab}
Solve by iteration.  Could use \verb|fsolve| for nonlinear
problems, but for linear it is much faster to use
Conjugate-Gradient algorithm.  \verb|gmres| is competitive,
but appears to take twice as long.
% uSoln=gmres(@(u) rhsb-theRes(u),rhsb,[],1e-9,maxIt);
\begin{matlab}
%}
tic;
if 0, uSoln=fsolve(@theRes,u0(i)); 
%{
\end{matlab}
The above is for nonlinear \pde{}s.  For linear \pde{}s,
determine the \textsc{rhs} vector, and make a function that
computes the matrix vector product.
\begin{matlab}
%}
else 
    maxIt = ceil(nVars/10);
    rhsb = theRes(u0(i));
    uSoln = pcg(@(u) rhsb-theRes(u),rhsb,1e-9,maxIt);
end
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
Draw the patch solution surface in the positive quadrant,
with edge-values omitted as already~\verb|NaN| by not
bothering to interpolate them.
\begin{matlab}
%}
surf(x(:),y(:),u'); view(60,40) 
maxu = ceil(max(u(:))*10)/10;
axis([0 1 0 1 0 maxu]),   caxis([0 maxu])
xlabel('x'), ylabel('y'), zlabel('u(x,y)')
%{
\end{matlab}



\paragraph{\texttt{theRes()}: function to zero}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system, and returns the vector of patch-interior time
derivatives.
\begin{matlab}
%}
function f=theRes(u)
  global i patches
  v=nan(size(patches.x+patches.y));
  v(i)=u;
  f=patchSys2(0,v(:),patches);
  f=f(i);
end
%{
\end{matlab}


\paragraph{\texttt{hetDiffForce2()}: heterogeneous diffusion PDE}
This function, based upon \cref{sec:heteroDiff2}, codes the
lattice heterogeneous diffusion of the \pde\ inside the
patches.  For 6D input arrays~\verb|u|, \verb|x|,
and~\verb|y|, computes the time derivativeat each point in
the interior of a patch, output in~\verb|ut|.  The two 2D
array of diffusivities,~$c^x_{ij}$ and~$c^y_{ij}$, are
stored in~\verb|patches.cs| (3D). 
\begin{matlab}
%}
function ut = hetDiffForce2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
  fu = -u+100*sin(pi*patches.x).*sin(pi*patches.y);
  ut(ix,iy,:,:,:,:) ...
  = diff(patches.cs(:,iy,1).*diff(u(:,iy,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,2).*diff(u(ix,:,:,:,:,:),1,2),1,2)/dy^2 ...
   +fu(ix,iy,:,:,:,:); 
end% function
%{
\end{matlab}

Fin.
%}
