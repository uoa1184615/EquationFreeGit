% Test the accuracy of the patch scheme with Dirichlet
% boundaries. Compare an example with various numbers of
% patches, the equilibrium in forced heterogeneous diffusion
% in 1D on patches.  Adapted the example from the second
% example of Eckhardt (2210.04536, sec 6.2.1).  Implement
% Dirichlet BCs using new facilities in the patch toolbox,
% and try various distributions of patches.  
% AJR, 29 Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{EckhardtEquilibErrs}: explore errors in
equilibria of a 1D heterogeneous diffusion on small patches}
\label{sec:EckhardtEquilibErrs}


\cref{sec:EckhardtEquilib} finds the equilibrium, of the
forced heterogeneous system with a forcing corresponding to
that applied at time \(t=1\).  Computational efficiency
comes from only computing the microscale heterogeneity on
small spatially sparse patches. Here we explore the errors
as the number~\(N\) of patches increases, see \cref{fig:EckhardtEquilibErrsus,fig:EckhardtEquilibErrs}. 
\begin{figure}
\centering\begin{tabular}{@{}c@{\ }c@{}}
\parbox[t]{10em}{\caption{\label{fig:EckhardtEquilibErrsus}%
Equilibrium of the heterogeneous diffusion problem for
relatively large \(\epsilon=0.03\) so we can see the
patches.  The solution is obtained with various numbers of
patches, but  we only compare solutions in these five common
patches.
}} &
\def\extraAxisOptions{mark size=1pt}
\raisebox{-\height}{\input{../Patch/Figs/EckhardtEquilibErrsus}}
\end{tabular}
\end{figure}%
Find mean-abs errors to be the following:
\begin{equation*}
\begin{array}{ccccccc}
&\qquad N=&5&9&17&33&65\\\hline
\text{equispace}&\text{second-order} 
&8\E3 &1\E2 &1\E2 &4\E3 &9\E4 
\\
\text{equispace}&\text{fourth-order} 
&2\E3 &7\E4 &1\E4 &9\E6 &5\E7 
\\
\text{equispace}&\text{sixth-order} 
&2\E3 &2\E5 &4\E7 &1\E8 &2\E{10} 
\\\hline
\text{chebyshev}&\text{second-order} 
&4\E2 &6\E2 &3\E2 &2\E2 &2\E2 
\\
\text{chebyshev}&\text{fourth-order} 
&9\E4 &3\E3 &6\E4 &3\E4 &2\E4 
\\
\text{chebyshev}&\text{sixth-order} 
&9\E4 &3\E5 &1\E5 &4\E6 &1\E6 
\\\hline
\text{usergiven}&\text{second-order} 
&4\E2 &6\E2 &3\E2 &9\E3 &2\E3 
\\
\text{usergiven}&\text{fourth-order} 
&8\E4 &3\E3 &6\E4 &4\E5 &2\E6 
\\
\text{usergiven}&\text{sixth-order} 
&8\E4 &3\E5 &1\E5 &2\E7 &3\E9
\\\hline
\end{array}
\end{equation*}
For `chebyshev' this assessment of errors is a bit dodgy as
it is based only on the centre and boundary patches. The
`usergiven' distribution is for overlapping patches with
Chebyshev distribution of centres---a spatial `christmas
tree'\footnote{But the error assessment is with respect to
finest patch-grid, no longer with a full domain solution}.
Curiously, and with above caveats, here my `smart' chebyshev
is the worst, the overlapping Chebyshev is good, but
\emph{equispace appears usually the best.}

\begin{figure}
\centering\begin{tabular}{@{}c@{\ }c@{}}
\parbox[t]{9em}{\caption{\label{fig:EckhardtEquilibErrs}%
Errors in the equilibrium of the heterogeneous diffusion
problem for relatively large \(\epsilon=0.03\).  The
solution is obtained with various numbers of patches, but we
only plot the errors within these five common patches.  
}} &
\def\extraAxisOptions{mark size=1pt}
\raisebox{-\height}{\input{../Patch/Figs/EckhardtEquilibErrs}}
\end{tabular}
\end{figure}
The above errors are for simple sin forcing.   What if we
make not so simple with exp modification of the forcing? 
The errors shown below are very little different (despite
the magnitude of the solution being a little larger).
\begin{equation*}
\begin{array}{ccccccc}
&\qquad N=&5&9&17&33&65
\\\hline
%\text{equispace}&\text{second-order} 
%&8\E3 &1\E2 &1\E2 &4\E3 &9\E4 
%\\
\text{equispace}&\text{fourth-order} 
&4\E3 &7\E4 &1\E4 &8\E6 &5\E7 
\\
%\text{equispace}&\text{sixth-order} 
%&2\E3 &2\E5 &4\E7 &1\E8 &2\E{10} 
%\\\hline
%\text{chebyshev}&\text{second-order} 
%&4\E2 &6\E2 &3\E2 &2\E2 &2\E2 
%\\
\text{chebyshev}&\text{fourth-order} 
&7\E4 &2\E3 &5\E4 &3\E4 &1\E4 
\\
%\text{chebyshev}&\text{sixth-order} 
%&9\E4 &3\E5 &1\E5 &4\E6 &1\E6 
%\\\hline
%\text{usergiven}&\text{second-order} 
%&4\E2 &6\E2 &3\E2 &9\E3 &2\E3 
%\\
\text{usergiven}&\text{fourth-order} 
&2\E3 &3\E3 &5\E4 &4\E5 &2\E6 
%\\
%\text{usergiven}&\text{sixth-order} 
%&8\E4 &3\E5 &1\E5 &2\E7 &3\E9
\\\hline
\end{array}
\end{equation*}



Clear, and initiate global patches.  Choose the type of
patch distribution to be either 'equispace', 'chebyshev', or
`usergiven'. Also set order of interpolation (fourth-order
is good start).
\begin{matlab}
%}
clear all
global patches 
%global OurCf2eps, OurCf2eps=true %option to save plots
switch 1
    case 1, Dom.type = 'equispace'
    case 2, Dom.type = 'chebyshev'
    case 3, Dom.type = 'usergiven'
    end% switch
ordInt = 4
%{
\end{matlab}


\paragraph{First configure the patch system}
Establish the microscale heterogeneity has micro-period
\verb|mPeriod| on the lattice, and coefficients to match
Eckhardt2210.04536 \S6.2.1. 
\begin{matlab}
%}
mPeriod = 6
z = (0.5:mPeriod)'/mPeriod;  
a = 1./(2-cos(2*pi*z))
global microTimePeriod; microTimePeriod=0;
%{
\end{matlab}

To use a hierarchy of patches with \verb|nPatch| of 5, 9,17,
\ldots, we need up to \(N\)~patches plus one~\verb|dx| to
fit into the domain interval.  Cater for up to some
full-domain simulation---can compute \(\verb|log2Nmax|=129\)
(\(\epsilon=0.008\)) in a few seconds:
\begin{matlab}
%}
log2Nmax = 7 % 5 for plots, 7 for choice
nPatchMax=2^log2Nmax+1
%{
\end{matlab}
Set the periodicity~\(\epsilon\), and other microscale
parameters. 
\begin{matlab}
%}
nPeriodsPatch = 1 % any integer
nSubP = nPeriodsPatch*mPeriod+2 % for edgy int
epsilon = 1/(nPatchMax*nPeriodsPatch+1/mPeriod)
dx = epsilon/mPeriod
%{
\end{matlab}

\paragraph{For various numbers of patches}
Choose five to be the coarsest number of patches.  Want
place to store common results for the solutions.  Assign
\verb|Ps| to be the indices of the common patches: for
equispace set to the five common patches, but for chebyshev
the only common ones are the three centre and
boundary-adjacent patches.
\begin{matlab}
%}
us=[]; xs=[]; nPs=[];
for log2N=log2Nmax:-1:2
    if log2N==log2Nmax
         Ps=linspace(1,nPatchMax ...
           ,5-2*all(Dom.type=='chebyshev'))
    else Ps=(Ps+1)/2
    end
%{
\end{matlab}
Set the number of patches in \((0,1)\):
\begin{matlab}
%}
    nPatch = 2^log2N+1
%{
\end{matlab}
In the case of `usergiven', we choose standard Chebyshev
distribution of the centre of the patches, which involves
overlapping of patches near the boundaries! (instead of the
coded chebyshev which has a boundary layer of
non-overlapping patches and a Chebyshev within the
interior). 
\begin{matlab}
%}
    if all(Dom.type=='usergiven')
      halfWidth=dx*(nSubP-1)/2;
      X1 = 0+halfWidth; X2 = 1-halfWidth;
      Dom.X = (X1+X2)/2-(X2-X1)/2*cos(linspace(0,pi,nPatch));
    end
%{
\end{matlab}
Configure the patches:
\begin{matlab}
%}
    configPatches1(@heteroDiffF,[0 1],Dom,nPatch ...
        ,ordInt,dx,nSubP,'EdgyInt',true,'hetCoeffs',a);
%{
\end{matlab}

Set the forcing coefficients, either the original parabolic,
or sinusoidal.  At time \(t=1\) the resultant forcing we
actually apply here is simply the sum of the two components.
\begin{matlab}
%}
    if 0 %given forcing gives exact answers for ordInt=4 !
      patches.f1 = 2*( patches.x-patches.x.^2 );
      patches.f2 = 2*0.5+0*patches.x;
    else% simple exp-sine forcing 
      patches.f1 = sin(pi*patches.x).*exp(patches.x);
      patches.f2 = pi/2*sin(pi*patches.x).*exp(patches.x);
    end%if
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
      u0 = zeros(nSubP,1,1,nPatch);
    else u0 = u0(:,:,:,1:2:end);
    end
    u0([1 end],:) = nan;  
    patches.i = find(~isnan(u0));
    nVariables = numel(patches.i)
%{
\end{matlab}
Solve via \verb|fsolve| for simplicity and robustness (and
using \verb|optimoptions| to omit trace information), via
the generic patch system wrapper \verb|theRes|
(\cref{sec:theRes}).
\begin{matlab}
%}
    tic;
    [uSoln,resSoln] = fsolve(@theRes,u0(patches.i) ...
        ,optimoptions('fsolve','Display','off'));
    fsolveTime = toc
%{
\end{matlab}
Store the solution into the patches, and give
magnitudes---Inf norm is max(abs()).
\begin{matlab}
%}
    normSoln = norm(uSoln,Inf)
    normResidual = norm(resSoln,Inf)
    u0(patches.i) = uSoln;
    u0 = patchEdgeInt1(u0);
    u0( 1 ,:,:, 1 ) = 0;
    u0(end,:,:,end) = 0;
%{
\end{matlab}
Concatenate the solution on common patches into stores.
\begin{matlab}
%}
    us=cat(3,us,squeeze(u0(:,:,:,Ps)));
    xs=cat(3,xs,squeeze(patches.x(:,:,:,Ps)));
    nPs = [nPs;nPatch];
%{
\end{matlab}

End loop.  Check grids were aligned, then compute errors
compared to the full-domain solution.
\begin{matlab}
%}
end%for log2N
assert(max(abs(reshape(diff(xs,1,3),[],1)))<1e-12,'x-coord failure')
errs = us-us(:,:,1);
meanAbsErrs = mean(abs(reshape(errs,[],size(us,3))))
ratioErrs = meanAbsErrs(2:end)./meanAbsErrs(1:end-1)
%{
\end{matlab}



\paragraph{Plot solution in common patches}
First adjoin NaNs to separate patches, and reshape.
\begin{matlab}
%}
x=xs(:,:,1); u=us;
x(end+1,:)=nan; u(end+1,:)=nan; 
u=reshape(u,numel(x),[]);
%{
\end{matlab}
Reshape solution field.
\begin{matlab}
%}
figure(1),clf
plot(x(:),u,'.-'), legend(num2str(nPs))
xlabel('space $x$'),ylabel('equilibrium $u(x)$')
ifOurCf2tex([mfilename 'us'])%optionally save
%{
\end{matlab}

\paragraph{Plot errors}
Use quasi-log axis to separate the errors.
\begin{matlab}
%}
err = u(:,1)-u;  
figure(2), clf
plot(x(:),err,'.-'); legend(num2str(nPs))
quasiLogAxes(10,sqrt(prod(meanAbsErrs(2:3))))
xlabel('space $x$'), ylabel('errors in $u(x)$')
ifOurCf2tex(mfilename)%optionally save
%{
\end{matlab}
%}