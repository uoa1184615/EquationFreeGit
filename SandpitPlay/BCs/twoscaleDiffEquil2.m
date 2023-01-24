% Solve for steady state of twoscale heterogeneous diffusion
% in 2D on patches as an example application, from section
% 5.3.1 of Freese, 2211.13731.  AJR, Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{twoscaleDiffEquil2}: equilibrium of a 2D
twoscale heterogeneous diffusion via small patches}
\label{sec:twoscaleDiffEquil2}

Here we find the steady state~\(u(x,y)\) to the
heterogeneous \pde\ (inspired by Freese et al.\footnote{
\protect \url{http://arxiv.org/abs/2211.13731}} \S5.3.1)
\begin{equation*}
u_t=A(x,y)\grad\grad u-f,
\end{equation*}
on domain \([-1,1]^2\) with Dirichlet BCs, for coefficient
`diffusion' matrix, varying with period~\(2\epsilon\) on the
microscale \(\epsilon=2^{-7}\), of
\begin{equation*}
A:=\begin{bmatrix} 2& a\\a & 2 \end{bmatrix}
\quad \text{with } a:=\sin(\pi x/\epsilon)\sin(\pi y/\epsilon),
\end{equation*}
and for forcing \(f:=(x+\cos3\pi x)y^3\).


Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches i
%{
\end{matlab}


First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the spatial lattice. Set the
phase of the heterogeneity so that each patch centre is a
point of symmetry of the diffusivity.  Then
\verb|configPatches2| replicates the heterogeneity to fill
each patch. 
\begin{matlab}
%}
mPeriod = 6
z = (0.5:mPeriod)'/mPeriod;  
A = sin(2*pi*z).*sin(2*pi*z');
%{
\end{matlab}
Set the periodicity, via~\(\epsilon\), and other microscale
parameters.
\begin{matlab}
%}
nPeriodsPatch = 1 % any integer
epsilon = 2^(-5) % so we can see patches
dx = (2*epsilon)/mPeriod
nSubP = nPeriodsPatch*mPeriod+2 % for edgy int
%{
\end{matlab}



\paragraph{Patch configuration}
Say use \(7\times7\) patches in \((-1,1)^2\), fourth order
interpolation, and either `equispace' or `chebyshev':
\begin{matlab}
%}
nPatch = 7 
configPatches2(@twoscaleDiffForce2,[-1 1],'equispace' ...
    ,nPatch ,4 ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',A );
%{
\end{matlab}
Compute the time-constant forcing, and store in struct
\verb|patches| for access by the microcode of
\cref{sec:twoscaleDiffForce2}.
\begin{matlab}
%}
patches.fu = 100*(patches.x+cos(3*pi*patches.x)).*patches.y.^3;
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






\subsection{\texttt{twoscaleDiffForce2()}: microscale
discretisation inside patches of forced diffusion PDE}
\label{sec:twoscaleDiffForce2}

This function codes the lattice heterogeneous diffusion of
the \pde\ inside the patches.  For 6D input arrays~\verb|u|,
\verb|x|, and~\verb|y|, computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.   
\begin{matlab}
%}
function ut = twoscaleDiffForce2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
Set Dirichlet boundary value of zero around the square
domain.
\begin{matlab}
%}
  u( 1 ,:,:,:, 1 ,:)=0; % left edge of left patches
  u(end,:,:,:,end,:)=0; % right edge of right patches
  u(:, 1 ,:,:,:, 1 )=0; % bottom edge of bottom patches
  u(:,end,:,:,:,end)=0; % top edge of top patches
%{
\end{matlab}
Compute the time derivatives via stored forcing and
coefficients.  Easier to code by conflating the last four
dimensions into the one~\verb|,:|. 
\begin{matlab}
%}
  ut(ix,iy,:) ...
  = 2*diff(u(:,iy,:),2,1)/dx^2 +2*diff(u(ix,:,:),2,2)/dy^2 ...
   +2*patches.cs(ix,iy).*( u(ix+1,iy+1,:) -u(ix-1,iy+1,:) ...
     -u(ix+1,iy-1,:) +u(ix-1,iy-1,:) )/(4*dx*dy) ...
   -patches.fu(ix,iy,:); 
end%function twoscaleDiffForce2
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
