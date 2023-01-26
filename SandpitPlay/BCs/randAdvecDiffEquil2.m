% Solve for steady state of two-scale heterogeneous
% diffusion in 2D on patches as an example application, from
% section 6.2 of Bonizzoni, 2211.15221.  AJR, Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{randAdvecDiffEquil2}: equilibrium of a 2D
random heterogeneous advection-diffusion via small patches}
\label{sec:randAdvecDiffEquil2}

Here we find the steady state~\(u(x,y)\) of the
heterogeneous \pde\ (inspired by Bonizzoni et al.\footnote{
\protect \url{http://arxiv.org/abs/2211.15221}} \S6.2)
\begin{equation*}
u_t=\mu_1\delsq u -(\cos\mu_2,\sin\mu_2)\cdot\grad u -u +f\,,
\end{equation*}
on domain \([0,1]^2\) with Neumann BCs, for microscale
random diffusion and advection coefficients,
\(\mu_1\in[0.01,0.1]\) and \(\mu_2\in[0,2\pi)\), and for
forcing
\begin{equation*}
f:=\exp\left[-\frac{(x-\mu_3)^2+(x-\mu_4)^2}{\mu_5^2}\right],
\end{equation*}
smoothly varying in space for fixed \(\mu_3, \mu_4 \in
[0.25,0.75]\) and \(\mu_5 \in [0.1,0.25]\). The above system
is dominantly diffusive for lengths scales \(\ell<0.01 =
\min\mu_1\).

Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches
%{
\end{matlab}


First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the spatial lattice.  Then
\verb|configPatches2| replicates the heterogeneity to fill
each patch. 
\begin{matlab}
%}
mPeriod = 4
mu1 = 10.^(-1-rand(mPeriod))
mu2 = 2*pi*rand(mPeriod)
cs = cat(3,mu1,cos(mu2),sin(mu2));
meanDiffAdvec=squeeze(mean(mean(cs)))
%{
\end{matlab}
Set the periodicity,~\(\epsilon\), and other microscale
parameters.
\begin{matlab}
%}
nPeriodsPatch = 1 % any integer
epsilon = 2^(-4) % so we can see patches
dx = epsilon/mPeriod
nSubP = nPeriodsPatch*mPeriod+2 % for edgy int
%{
\end{matlab}



\paragraph{Patch configuration}
Say use \(7\times7\) patches in \((0,1)^2\), fourth order
interpolation, either `equispace' or `chebyshev', and the
offset for Neumann boundary conditions:
\begin{matlab}
%}
nPatch = 7 
Dom.type= 'equispace';
Dom.bcOffset = 0.5; 
configPatches2(@randAdvecDiffForce2,[0 1],Dom ...
    ,nPatch ,4 ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',cs );
%{
\end{matlab}
Compute the time-constant forcing, and store in struct
\verb|patches| for access by the microcode of
\cref{sec:randAdvecDiffForce2}.
\begin{matlab}
%}
mu = [ 0.25+0.5*rand(1,2) 0.1+0.15*rand ]
patches.fu = exp(-((patches.x-mu(1)).^2+(patches.y-mu(2)).^2)/mu(3)^2);
%{
\end{matlab}




\paragraph{Solve for steady state}
Set initial guess of zero, with \verb|NaN| to indicate
patch-edge values.  Index~\verb|i| are the indices of
patch-interior points, store in global patches for access by
\verb|theRes|, and the number of unknowns is then its number
of elements.
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
information).
\begin{matlab}
%}
tic;
uSoln = fsolve(@theRes,u0(patches.i) ...
        ,optimoptions('fsolve','Display','off')); 
solnTime = toc
%{
\end{matlab}
Store the solution into the patches, and give magnitudes.
\begin{matlab}
%}
u0(patches.i) = uSoln;
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






\subsection{\texttt{randAdvecDiffForce2()}: microscale
discretisation inside patches of forced diffusion PDE}
\label{sec:randAdvecDiffForce2}

This function codes the lattice heterogeneous diffusion of
the \pde\ inside the patches.  For 6D input arrays~\verb|u|,
\verb|x|, and~\verb|y|, computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.   
\begin{matlab}
%}
function ut = randAdvecDiffForce2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
Set Neumann boundary condition of zero derivative around the square
domain.
\begin{matlab}
%}
  u( 1 ,:,:,:, 1 ,:)=u(  2  ,:,:,:, 1 ,:); % left edge of left patches
  u(end,:,:,:,end,:)=u(end-1,:,:,:,end,:); % right edge of right patches
  u(:, 1 ,:,:,:, 1 )=u(:,  2  ,:,:,:, 1 ); % bottom edge of bottom patches
  u(:,end,:,:,:,end)=u(:,end-1,:,:,:,end); % top edge of top patches
%{
\end{matlab}
Compute the time derivatives via stored forcing and
coefficients.  Easier to code by conflating the last four
dimensions into the one~\verb|,:|. 
\begin{matlab}
%}
  ut(ix,iy,:) ...
  = patches.cs(ix,iy,1).*(diff(u(:,iy,:),2,1)/dx^2 ...
                         +diff(u(ix,:,:),2,2)/dy^2)...
   -patches.cs(ix,iy,2).*(u(ix+1,iy,:)-u(ix-1,iy,:))/(2*dx) ...
   -patches.cs(ix,iy,3).*(u(ix,iy+1,:)-u(ix,iy-1,:))/(2*dy) ...
   -u(ix,iy,:) +patches.fu(ix,iy,:); 
end%function randAdvecDiffForce2
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
  global patches
  v=nan(size(patches.x+patches.y));
  v(patches.i)=u;
  f=patchSys2(0,v(:),patches);
  f=f(patches.i);
end%function theRes
%{
\end{matlab}

Fin.
%}
