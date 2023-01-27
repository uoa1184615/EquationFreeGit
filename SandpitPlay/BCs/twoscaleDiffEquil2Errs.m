% Explore errors inn the steady state of twoscale
% heterogeneous diffusion in 2D on patches as an example,
% inspired by section 5.3.1 of Freese et al., 2211.13731. 
% AJR, Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{twoscaleDiffEquil2Errs}: errors in
equilibria of a 2D twoscale heterogeneous diffusion via
small patches}
\label{sec:twoscaleDiffEquil2Errs}

Here we find the steady state~\(u(x,y)\) to the
heterogeneous \pde\ (inspired by Freese et al.\footnote{
\protect \url{http://arxiv.org/abs/2211.13731}} \S5.3.1)
\begin{equation*}
u_t=A(x,y)\grad\grad u+f,
\end{equation*}
on domain \([-1,1]^2\) with Dirichlet BCs, for coefficient
`diffusion' matrix, varying with some microscale
period~\(\epsilon\) (here \(\epsilon\approx 0.12,0.06\)), of
\begin{equation*}
A:=\begin{bmatrix} 2& a\\a & 2 \end{bmatrix}
\quad \text{with } a:=\sin(\pi x/\epsilon)\sin(\pi y/\epsilon),
\end{equation*}
and for forcing \(f:=10(x+y+\cos\pi x)\) (for which the
solution has magnitude up to one).\footnote{Freese et al.\
had forcing \(f:=(x+\cos3\pi x)y^3\), but here we want
smoother forcing so we get meaningful results on few patches
in up to a minute or two computation.  For the same reason
we do not invoke their smaller \(\epsilon\approx 0.01\).}

Here we explore the errors for increasing number~\(N\) of
patches (in both directions). Find mean-abs errors to be the
following acceptable decrease (for different orders of
interpolation):
\begin{equation*}
\def\ten#1{\cdot10^{-#1}}
\begin{array}{clll}
N&5&9&17%&33&65
\\\hline
\text{second-order} &0.069 &0.034 &0.012  \\
\text{fourth-order} &0.033 &0.0086 &0.00069 
\end{array}
\end{equation*}
Need more expensive larger~\(N\) to clarify the rate of decay.

Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches 
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
To use a hierarchy of patches with \verb|nPatch| of 5, 9,17,
\ldots, we need up to \(N\)~patches plus one~\verb|dx| to
fit into the domain interval.  Cater for up to some
full-domain simulation---can compute \(\verb|log2Nmax|=5\)
(\(\epsilon=0.06\)) in a couple of minutes:
\begin{matlab}
%}
log2Nmax = 4
nPatchMax=2^log2Nmax+1
%{
\end{matlab}
Set the periodicity~\(\epsilon\), and other microscale
parameters. 
\begin{matlab}
%}
nPeriodsPatch = 1 % any integer
nSubP = nPeriodsPatch*mPeriod+2 % for edgy int
epsilon = 2/(nPatchMax*nPeriodsPatch+1/mPeriod)
dx = epsilon/mPeriod
ordInt = 4
%{
\end{matlab}

\paragraph{For various numbers of patches}
Assume five patches is the coarsest patches.  Want place to
store common results for the solutions.  Assign \verb|Ps| to
be the indices of the common patches
\begin{matlab}
%}
us=[]; xs=[]; ys=[];
for log2N=log2Nmax:-1:2
    if log2N==log2Nmax
         Ps=1:2^(log2N-2):nPatchMax
    else Ps=(Ps+1)/2
    end
%{
\end{matlab}
Use patches in \((-1,1)^2\), fourth order
interpolation, and either `equispace' or `chebyshev':
\begin{matlab}
%}
nPatch = 2^log2N+1
configPatches2(@twoscaleDiffForce2,[-1 1],'equispace' ...
    ,nPatch ,ordInt ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',A );
%{
\end{matlab}
Compute the time-constant forcing, and store in struct
\verb|patches| for access by the microcode of
\cref{sec:twoscaleDiffForce2}.
\begin{matlab}
%}
if 1
  patches.fu = 10*(patches.x+cos(pi*patches.x)+patches.y);
else patches.fu = 8+0*patches.x+0*patches.y;
end
%{
\end{matlab}

\paragraph{Solve for steady state}
Set initial guess of either zero or a subsample of the next
finer solution, with \verb|NaN| to indicate patch-edge
values.  Index~\verb|i| are the indices of patch-interior
points, and the number of unknowns is then its length.
\begin{matlab}
%}
if log2N==log2Nmax
  u0 = zeros(nSubP,nSubP,1,1,nPatch,nPatch);
else u0 = u0(:,:,:,:,1:2:end,1:2:end);
end
u0([1 end],:,:) = nan;  u0(:,[1 end],:) = nan;
patches.i = find(~isnan(u0));
nVariables = numel(patches.i)
%{
\end{matlab}
Solve via iterative solver \verb|bicgstab| (or \verb|gmres|)
first, and if that fails then use \verb|fsolve| for
simplicity and robustness (and using \verb|optimoptions| to
omit trace information).
\begin{matlab}
%}
tic;
maxIt = ceil(nVariables/10);
rhsb = theRes(u0(patches.i));
[uSoln,flag] = bicgstab(@(u) rhsb-theRes(u),rhsb,1e-9,maxIt);
bicgTime = toc
if flag>0, disp('bicg failed, trying fsolve')
    tic;
    uSoln = fsolve(@theRes,u0(patches.i) ...
        ,optimoptions('fsolve','Display','off'));
    fsolveTime = toc
end%if flag
%{
\end{matlab}
Store the solution into the patches, and give
magnitudes---Inf norm is max(abs()).
\begin{matlab}
%}
normSoln = norm(uSoln,Inf)
normResidual = norm(theRes(uSoln),Inf)
u0(patches.i) = uSoln;
u0 = patchEdgeInt2(u0);
%{
\end{matlab}
Concatenate the solution on common patches into stores.
\begin{matlab}
%}
us=cat(5,us,squeeze(u0(:,:,:,:,Ps,Ps)));
xs=cat(3,xs,squeeze(patches.x(:,:,:,:,Ps,:)));
ys=cat(3,ys,squeeze(patches.y(:,:,:,:,:,Ps)));
%{
\end{matlab}

End loop.  Check grids were aligned, then compute errors
compared to the full-domain solution.
\begin{matlab}
%}
end%for log2N
assert(max(abs(reshape(diff(xs,1,3),[],1)))<1e-12,'x-coord failure')
assert(max(abs(reshape(diff(ys,1,3),[],1)))<1e-12,'y-coord failure')
errs = us-us(:,:,:,:,1);
meanAbsErrs = mean(abs(reshape(errs,[],size(us,5))))
ratioErrs = meanAbsErrs(2:end)./meanAbsErrs(1:end-1)
%{
\end{matlab}




\paragraph{Plot solution in common patches}
First reshape arrays to suit 2D space surface plots,
inserting nans to separate patches.
\begin{matlab}
%}
x = xs(:,:,1); y = ys(:,:,1); u=us;
x(end+1,:)=nan; y(end+1,:)=nan;
u(end+1,:,:)=nan; u(:,end+1,:)=nan;
u = reshape(permute(u,[1 3 2 4 5]),numel(x),numel(y),[]);
%{
\end{matlab}
Plot the patch solution surfaces, with colour offset between
surfaces (best if \(u\)-field has a range of one): blues are
the full-domain solution, reds the coarsest patches.
\begin{matlab}
%}
figure(1), clf, colormap(jet)
for p=1:size(u,3)
     mesh(x(:),y(:),u(:,:,p)',p-1+u(:,:,p)'); 
     hold on;
end, hold off
view(60,55), colorbar
xlabel('x'), ylabel('y'), zlabel('u(x,y)')
%{
\end{matlab}

\paragraph{Plot error surfaces}
Plot the error surfaces, with colour offset between surfaces
(best if \(u\)-field has a range of one): dark blue is the
full-domain zero error, reds the coarsest patches.
\begin{matlab}
%}
err=u-u(:,:,1);  
maxAbsErr=max(abs(err(:)))
figure(2), clf, colormap(jet)
for p=1:size(u,3)
     mesh(x(:),y(:),err(:,:,p)',p-1+err(:,:,p)'/maxAbsErr/2); 
     hold on;
end, hold off
view(60,55)
xlabel('x'), ylabel('y'), zlabel('errors in u(x,y)')
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
   +patches.fu(ix,iy,:); 
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
