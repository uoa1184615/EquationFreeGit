% Solve for steady state of monoscale heterogeneous
% diffusion in 2D on patches as an example application, from
% section 5.2 of Freese, 2211.13731.  AJR, Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{monoscaleDiffEquil2}: equilibrium of a 2D
monoscale heterogeneous diffusion via small patches}
\label{sec:monoscaleDiffEquil2}

Here we find the steady state~\(u(x,y)\) to the
heterogeneous \pde (inspired by Freese et al.\footnote{
\protect \url{http://arxiv.org/abs/2211.13731}} \S5.2)
\begin{equation*}
u_t=A(x,y)\grad\grad u-f,
\end{equation*}
on domain \([-1,1]^2\) with Dirichlet BCs, for coefficient
`diffusion' matrix
\begin{equation*}
A:=\begin{bmatrix} 2& a\\a & 2 \end{bmatrix}
\quad \text{with } a:=\sign(xy)
\text{ or }a:=\sin(\pi x)\sin(\pi y),
\end{equation*}
and for forcing~\(f(x,y)\) such that the exact equilibrium is 
\begin{equation*}
u=x\big(1-e^{1-|x|}\big)y\big(1-e^{1-|y|}\big).
\end{equation*}
But for simplicity, let's do \(u=x(1-x^2)y(1-y^2)\) for
which we code~\(f\) later---as determined by this computer
algebra.

\begin{verbatim}
on gcd; factor sin;
%let { df(sign(~x),~x)=>0
%    , df(abs(~x),~x)=>sign(x) 
%    , abs(~x)^2=>abs(x), sign(~x)^2=>1 };
%u:=x*(1-exp(1-abs(x)))*y*(1-exp(1-abs(y)));
u:=x*(1-x^2)*y*(1-y^2);
a:=sin(pi*x)*sin(pi*y);
f:=2*df(u,x,x)+2*a*df(u,x,y)+2*df(u,y,y);
\end{verbatim}


Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches i
%{
\end{matlab}


\paragraph{Patch configuration}
Initially use \(7\times7\) patches in the square \((-1,1)^2\). For
continuous forcing we may have small patches of any
reasonable microgrid spacing---here the microgrid error
dominates.
\begin{matlab}
%}
nPatch = 7
nSubP = 5 
dx = 0.03
%{
\end{matlab}
Specify some order of interpolation.  
\begin{matlab}
%}
configPatches2(@monoscaleDiffForce2,[-1 1 -1 1],'equispace' ...
    ,nPatch ,4 ,dx ,nSubP ,'EdgyInt',true );
%{
\end{matlab}
Compute the time-constant coefficient and time-constant
forcing, and store them in struct \verb|patches| for access
by the microcode of \cref{sec:monoscaleDiffForce2}.
\begin{matlab}
%}
  x=patches.x;  y=patches.y;
  patches.A = sin(pi*x).*sin(pi*y);
  patches.fu = ...
    +2*patches.A.*(9*x.^2.*y.^2-3*x.^2-3*y.^2+1) ...
    +12*x.*y.*(x.^2+y.^2-2);
%{
\end{matlab}
By construction, the \pde\ has analytic solution
\begin{matlab}
%}
uAnal = x.*(1-x.^2).*y.*(1-y.^2);
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
nVars = numel(i)
%{
\end{matlab}
Solve by iteration.  Use \verb|fsolve| for simplicity and
robustness (using \verb|optimoptions| to omit its trace
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
errors = uAnal(i)-uSoln;
normError = norm(errors)
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



\subsection{\texttt{monoscaleDiffForce2()}: microscale
discretisation inside patches of forced diffusion PDE}
\label{sec:monoscaleDiffForce2}

This function codes the lattice heterogeneous diffusion of
the \pde\ inside the patches.  For 6D input arrays~\verb|u|,
\verb|x|, and~\verb|y|, computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.   
\begin{matlab}
%}
function ut = monoscaleDiffForce2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
Set Dirichlet boundary value of zero around the square
domain, or code some function variation.
\begin{matlab}
%}
u( 1 ,:,:,:, 1 ,:)=0; % left edge of left patches
u( 1 ,:,:,:, 1 ,:)=(1+patches.y)/2; % or code function of y
u(end,:,:,:,end,:)=0; % right edge of right patches
u(:, 1 ,:,:,:, 1 )=0; % bottom edge of bottom patches
u(:,end,:,:,:,end)=0; % top edge of top patches
u(:,end,:,:,:,end)=1; % or code function of x
%{
\end{matlab}
Compute the time derivatives via stored forcing and
coefficients.  Easier to code by conflating the last four
dimensions into the one~\verb|,:|. 
\begin{matlab}
%}
  ut(ix,iy,:) ...
  = 2*diff(u(:,iy,:),2,1)/dx^2 +2*diff(u(ix,:,:),2,2)/dy^2 ...
   +2*patches.A(ix,iy,:).*( u(ix+1,iy+1,:) -u(ix-1,iy+1,:) ...
     -u(ix+1,iy-1,:) +u(ix-1,iy-1,:) )/(4*dx*dy) ...
   -patches.fu(ix,iy,:); 
end%function monoscaleDiffForce2
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
