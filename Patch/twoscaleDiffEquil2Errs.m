% Explore errors in the steady state of twoscale
% heterogeneous diffusion in 2D on patches as an example,
% inspired by section 5.3.1 of Freese et al., 2211.13731.
% AJR, 31 Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{twoscaleDiffEquil2Errs}: errors in
equilibria of a 2D twoscale heterogeneous diffusion via
small patches}
\label{sec:twoscaleDiffEquil2Errs}

\begin{figure}
\centering\caption{\label{fig:twoscaleDiffEquil2Errsus} For
various numbers of patches as indicated on the colorbar,
plot the equilibrium of the multiscale diffusion problem of
Freese with Dirichlet zero-value boundary conditions
(\cref{sec:twoscaleDiffEquil2Errs}).  We only compare
solutions only in these 25~common patches.}
\includegraphics[scale=0.8]{Figs/twoscaleDiffEquil2Errsus}
\end{figure}%
Here we find the steady state~\(u(x,y)\) to the
heterogeneous \pde\ (inspired by Freese et al.\footnote{
\protect \url{http://arxiv.org/abs/2211.13731}} \S5.3.1)
\begin{equation*}
u_t=A(x,y)\grad\grad u+f,
\end{equation*}
on domain \([-1,1]^2\) with Dirichlet BCs, for coefficient
`diffusion' matrix, varying with some microscale
period~\(\epsilon\) (here \(\epsilon\approx 0.24, 0.12,
0.06, 0.03\)), of
\begin{equation*}
A:=\begin{bmatrix} 2& a\\a & 2 \end{bmatrix}
\quad \text{with } a:=\sin(\pi x/\epsilon)\sin(\pi y/\epsilon),
\end{equation*}
and for forcing \(f:=10(x+y+\cos\pi x)\) (for which the
solution has magnitude up to one).\footnote{Freese et al.\
had forcing \(f:=(x+\cos3\pi x)y^3\), but here we want
smoother forcing so we get meaningful results in a minute or
two computation.\footnote{Except in the `usergiven' case,
for $N=65$, that is $152\,100$ unknowns, it takes an hour to
compute the Jacobian, then chokes.}  For the same reason we
do not invoke their smaller \(\epsilon\approx 0.01\).}

\begin{figure}
\centering\caption{\label{fig:twoscaleDiffEquil2Errs} For
various numbers of patches as indicated on the colorbar,
plot the equilibrium of the multiscale diffusion problem of
Freese with Dirichlet zero-value boundary conditions
(\cref{sec:twoscaleDiffEquil2Errs}).  We only compare
solutions only in these 25~common patches.}
\includegraphics[scale=0.8]{Figs/twoscaleDiffEquil2Errs}
\end{figure}%
Here we explore the errors for increasing number~\(N\) of
patches (in both directions). Find mean-abs errors to be the
following (for different orders of interpolation and patch
distribution):
\begin{equation*}
\begin{array}{rcccc}
N&5&9&17&33%&65
\\\hline
\text{equispace, 2nd-order} &6\E2 &3\E2 &1\E2 &3\E3
\\
\text{equispace, 4th-order} &3\E2 &8\E3 &7\E4 &7\E5
\\
\text{chebyshev, 4th-order} &1\E2 &2\E2 &6\E3 &2\E3
\\
\text{usergiven, 4th-order} &1\E2 &2\E2 &4\E3 &\text{n/a}
\\
\text{equispace, 6th-order} &3\E2 &1\E3 &1\E4 &2\E5
\\\hline
\end{array}
\end{equation*}


\paragraph{Script start} Clear, and initiate global patches.
Choose the type of patch distribution to be either
`equispace', `chebyshev', or `usergiven'. Also set order of
interpolation (fourth-order is good start).
\begin{matlab}
%}
clear all
global patches 
%global OurCf2eps, OurCf2eps=true %option to save plot
switch 1
    case 1, Dom.type = 'equispace'
    case 2, Dom.type = 'chebyshev'
    case 3, Dom.type = 'usergiven'
end% switch
ordInt = 4
%{
\end{matlab}


\paragraph{First configure the patch system} Establish the
microscale heterogeneity has micro-period \verb|mPeriod| on
the spatial lattice. Then \verb|configPatches2| replicates
the heterogeneity as needed to fill each patch. 
\begin{matlab}
%}
mPeriod = 6
z = (0.5:mPeriod)'/mPeriod;  
A = sin(2*pi*z).*sin(2*pi*z');
%{
\end{matlab}
To use a hierarchy of patches with \verb|nPatch| of~5, 9,
17, \ldots, we need up to \(N\)~patches plus one~\verb|dx|
to fit into the domain interval.  Cater for up to some
full-domain simulation---can compute \(\verb|log2Nmax|=5\)
(\(\epsilon=0.06\)) within minutes:
\begin{matlab}
%}
log2Nmax = 4 % >2 up to 6 OKish
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
%{
\end{matlab}

\paragraph{For various numbers of patches} Choose five
patches to be the coarsest number of patches.  Define
variables to store common results for the solutions from
differing patches.  Assign \verb|Ps| to be the indices of
the common patches: for equispace set to the five common
patches, but for `chebyshev' the only common ones are the
three centre and boundary-adjacent patches.
\begin{matlab}
%}
us=[]; xs=[]; ys=[]; nPs=[];
for log2N=log2Nmax:-1:2
    if log2N==log2Nmax
         Ps=linspace(1,nPatchMax ...
           ,5-2*all(Dom.type=='chebyshev'))
    else Ps=(Ps+1)/2
    end
%{
\end{matlab}
Set the number of patches in \((-1,1)\):
\begin{matlab}
%}
    nPatch = 2^log2N+1
%{
\end{matlab}
In the case of `usergiven', we set the standard Chebyshev
distribution of the patch-centres, which involves
overlapping of patches near the boundaries! (instead of the
coded chebyshev which has boundary layers of abutting
patches, and non-overlapping Chebyshev between the boundary
layers). 
\begin{matlab}
%}
    if all(Dom.type=='usergiven')
      halfWidth = dx*(nSubP-1)/2;
      X1 = -1+halfWidth;  X2 = 1-halfWidth;
      Dom.X = (X1+X2)/2-(X2-X1)/2*cos(linspace(0,pi,nPatch));
      Dom.Y = Dom.X;
    end
%{
\end{matlab}
Configure the patches:
\begin{matlab}
%}
    configPatches2(@twoscaleDiffForce2,[-1 1],Dom,nPatch  ...
        ,ordInt ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',A );
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

\paragraph{Solve for steady state} Set initial guess of
either zero or a subsample of the previous, next-finer,
solution. \verb|NaN| indicates patch-edge values. 
Index~\verb|i| are the indices of patch-interior points, and
the number of unknowns is then its length.
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
First try to solve via iterative solver \verb|bicgstab|, via
the generic patch system wrapper \verb|theRes|
(\cref{sec:theRes}).
\begin{matlab}
%}
    tic;
    maxIt = ceil(nVariables/10);
    rhsb = theRes( zeros(size(patches.i)) );
    [uSoln,flag] = bicgstab(@(u) rhsb-theRes(u),rhsb ...
                   ,1e-9,maxIt,[],[],u0(patches.i));
    bicgTime = toc
%{
\end{matlab}
However, the above often fails (and \verb|fsolve| sometimes
takes too long here), so then try a preconditioned version
of \verb|bicgstab|.  The preconditioner is derived from the
Jacobian which is expensive to find (four minutes for
\(N=33\), one hour for $N=65$), but we do so as follows. 
\begin{matlab}
%}
    if flag>0, disp('**** bicg failed, trying ILU preconditioner')
        disp(['Computing Jacobian: wait roughly ' ...
              num2str(nPatch^4/4500,2) ' secs'])
        tic  
        Jac=sparse(nVariables,nVariables);
        for j=1:nVariables
            Jac(:,j)=sparse( rhsb-theRes((1:nVariables)'==j) );
        end
        formJacTime=toc
%{
\end{matlab}
Compute an incomplete \(LU\)-factorization, and use it as
preconditioner to \verb|bicgstab|.
\begin{matlab}
%}
        tic
        [L,U] = ilu(Jac,struct('type','ilutp','droptol',1e-4));
        LUfillFactor = (nnz(L)+nnz(U))/nnz(Jac)
        [uSoln,flag] = bicgstab(@(u) rhsb-theRes(u),rhsb ...
                   ,1e-9,maxIt,L,U,u0(patches.i));
        precondSolveTime=toc
        assert(flag==0,'preconditioner fails bicgstab. Lower droptol?')
    end%if flag
%{
\end{matlab}
Store the solution into the patches, and give
magnitudes---Inf norm is max(abs()).
\begin{matlab}
%}
    normResidual = norm(theRes(uSoln),Inf)
    normSoln = norm(uSoln,Inf)
    u0(patches.i) = uSoln;
    u0 = patchEdgeInt2(u0);
    u0( 1 ,:,:,:, 1 ,:)=0; % left edge of left patches
    u0(end,:,:,:,end,:)=0; % right edge of right patches
    u0(:, 1 ,:,:,:, 1 )=0; % bottom edge of bottom patches
    u0(:,end,:,:,:,end)=0; % top edge of top patches
    assert(normResidual<1e-5,'poor--bad solution found')
%{
\end{matlab}
Concatenate the solution on common patches into stores.
\begin{matlab}
%}
    us=cat(5,us,squeeze(u0(:,:,:,:,Ps,Ps)));
    xs=cat(3,xs,squeeze(patches.x(:,:,:,:,Ps,:)));
    ys=cat(3,ys,squeeze(patches.y(:,:,:,:,:,Ps)));
    nPs = [nPs;nPatch];
%{
\end{matlab}

End loop.  Check micro-grids are aligned, then compute
errors compared to the full-domain solution (or the highest
resolution solution for the case of `usergiven').
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




\paragraph{Plot solution in common patches} First reshape
arrays to suit 2D space surface plots, inserting nans to
separate patches.
\begin{matlab}
%}
x = xs(:,:,1); y = ys(:,:,1); u=us;
x(end+1,:)=nan;   y(end+1,:)=nan;
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
     mesh(x(:),y(:),u(:,:,p)',p+u(:,:,p)'); 
     hold on;
end, hold off
view(60,55) 
colorbar('Ticks',1:size(u,3) ...
    ,'TickLabels',[num2str(nPs) ['x';'x';'x'] num2str(nPs)]);
xlabel('space x'), ylabel('space y'), zlabel('u(x,y)')
ifOurCf2eps([mfilename 'us'])%optionally save
%{
\end{matlab}



\paragraph{Plot error surfaces} Plot the error surfaces,
with colour offset between surfaces (best if \(u\)-field has
a range of one): dark blue is the full-domain zero error,
reds the coarsest patches.
\begin{matlab}
%}
err=u(:,:,1)-u;  
maxAbsErr=max(abs(err(:)));
figure(2), clf, colormap(jet)
for p=1:size(u,3)
     mesh(x(:),y(:),err(:,:,p)',p+err(:,:,p)'/maxAbsErr); 
     hold on;
end, hold off
view(60,55)
colorbar('Ticks',1:size(u,3) ...
    ,'TickLabels',[num2str(nPs) ['x';'x';'x'] num2str(nPs)]);
xlabel('space x'), ylabel('space y')
zlabel('errors in u(x,y)')
ifOurCf2eps(mfilename)%optionally save
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
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
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
  ut(i,j,:) ...
  = 2*diff(u(:,j,:),2,1)/dx^2 +2*diff(u(i,:,:),2,2)/dy^2 ...
   +2*patches.cs(i,j).*( u(i+1,j+1,:) -u(i-1,j+1,:) ...
     -u(i+1,j-1,:) +u(i-1,j-1,:) )/(4*dx*dy) ...
   +patches.fu(i,j,:); 
end%function twoscaleDiffForce2
%{
\end{matlab}
%}
