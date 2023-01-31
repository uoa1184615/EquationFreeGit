% Test the accuracy of the patch scheme with Dirichlet boundaries. 
% Compare an example with various numbers of patches, the equilibrium in forced heterogeneous diffusion in
% 1D on patches.  Adapted the example from the second example of Eckhardt
% (2210.04536, sec 6.2.1).  Implement Dirichlet BCs using
% new facilities in the patch toolbox.   AJR, 27 Jan 2023
%!TEX root = doc.tex
%{
\section{\texttt{EckhartEquilibErrs}: find an equilibrium of a 1D
heterogeneous diffusion via small patches}
\label{sec:EckhartEquilibErrs}

\cref{fig:EckhardtEquilib} finds the equilibrium, of the forced heterogeneous
system with a forcing corresponding to that applied at time
\(t=1\).  Computational efficiency comes from only computing
the microscale heterogeneity on small spatially sparse
patches.  
Here we explore the errors as the number of patches increases.

Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches 
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
a = 1./(2-1.5*cos(2*pi*z))
global microTimePeriod; microTimePeriod=0;
%{
\end{matlab}

To use a hierarchy of patches with \verb|nPatch| of 5, 9,17,
\ldots, we need up to \(N\)~patches plus one~\verb|dx| to
fit into the domain interval.  Cater for up to some
full-domain simulation---can compute \(\verb|log2Nmax|=129\)
(\(\epsilon=0.008\)) within a minute:
\begin{matlab}
%}
log2Nmax = 5
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
Select maybe fourth-order interpolation.
\begin{matlab}
%}
ordInt = 4
%{
\end{matlab}

\paragraph{For various numbers of patches}
Assume five patches is the coarsest patches.  Want place to
store common results for the solutions.  Assign \verb|Ps| to
be the indices of the common patches
\begin{matlab}
%}
us=[]; xs=[]; nPs=[];
for log2N=log2Nmax:-1:2
    if log2N==log2Nmax
         Ps=1:2^(log2N-2):nPatchMax
    else Ps=(Ps+1)/2
    end
%{
\end{matlab}
Use patches in \((0,1)\),  and either `equispace' or `chebyshev':
\begin{matlab}
%}
nPatch = 2^log2N+1
configPatches1(@heteroDiffF,[0 1],'equispace',nPatch ...
    ,ordInt,dx,nSubP,'EdgyInt',true,'hetCoeffs',a);
%{
\end{matlab}

Set the forcing coefficients, either the original parabolic,
or sinusoidal.  At time \(t=1\) the resultant forcing we
actually apply here is simply the sum of the two components.
\begin{matlab}
%}
if 0 %given forcing gives exact answers for ordInt=4 !!!
  patches.f1 = 2*( patches.x-patches.x.^2 );
  patches.f2 = 2*0.5+0*patches.x;
else% simple sine forcing 
  patches.f1 = sin(pi*patches.x);
  patches.f2 = pi/2*sin(pi*patches.x);
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
Solve via iterative solver \verb|bicgstab| (or \verb|gmres|)
first, and if that fails then use \verb|fsolve| for
simplicity and robustness (and using \verb|optimoptions| to
omit trace information).
\begin{matlab}
%}
tic;
rhsb = theRes1(zeros(size(patches.i)));
%{
\end{matlab}
This preconditioner improves the condition number of the matrix by a factor of roughly three, but that is not enough to get convergence.
\begin{matlab}
%}
if nSubP-2<nPatch % try constructing a preconditioner
  v=eye(nSubP-2,nPatch); % the interior points
  Av=theRes1(v)-rhsb;
  microA=reshape(Av(1:(nSubP-2)^2),nSubP-2,nSubP-2)
  condMicroA=cond(microA)
  mFun = @(u) reshape( microA\ reshape(u ...
              ,nSubP-2,[]) ,size(u,1),[]);
% construct Jacobian explicitly to test preconditioner
nJac=(nSubP-2)*nPatch
Jac=nan(nJac);
mustbezero=norm(theRes1(zeros(nJac,1))-rhsb)
assert(mustbezero<1e-9,'failed zero residual')
for j=1:nJac;
  Jac(:,j)=theRes1((1:nJac)'==j)-rhsb;
end
figure(3),spy(abs(Jac)>1e-9)
condJac=cond(Jac)
relativeCondMJac=cond(mFun(Jac))/condJac
pause(1)
% finished with Jacobian
else mFun = @(u) u;
end% if nSubP 
maxIt = ceil(nVariables/10);
[uSoln,flag] = bicgstab(@(u) theRes1(u)-rhsb,rhsb,1e-9,maxIt,mFun);
%restart=ceil(nVariables/10);
%[uSoln,flag] = gmres(@(u) theRes1(u)-rhsb,rhsb,[],1e-9,[],mFun);
bicgstabFlag=flag
bicgTime = toc
if flag>0, disp('**** bicg failed, trying fsolve')
    tic;
    uSoln = fsolve(@theRes1,u0(patches.i) ...
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
normResidual = norm(theRes1(uSoln),Inf)
u0(patches.i) = uSoln;
u0 = patchEdgeInt1(u0);
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
First clear figure, and adjoin NaNs to separate patches.
\begin{matlab}
%}
figure(1),clf
x=xs(:,:,1); u=us;
x(end+1,:)=nan; u(end+1,:)=nan; 
%{
\end{matlab}
Reshape solution field.
\begin{matlab}
%}
u=reshape(u,numel(x),[]);
plot(x(:),u,'.-'), legend(num2str(nPs))
xlabel('space $x$'),ylabel('equilibrium $u(x)$')
%{
\end{matlab}

\paragraph{Plot errors}
Use quasi-log axis to separate the errors.
\begin{matlab}
%}
err = u(:,1)-u;  
maxAbsErr = max(abs(err(:)))
figure(2), clf
h=plot(x(:),err,'.-'); legend(num2str(nPs))
quasiLogAxes(h,10,1e-5)
xlabel('x'), ylabel('errors in $u(x)$')
%{
\end{matlab}

%Optionally write to graphic file.
%\begin{matlab}
%%}
%matlab2tikz('Figs/EckhardtEquilib.tex','showInfo',false ...
%,'noSize',true,'parseStrings',false,'showWarnings',false ...
%,'extraCode','\tikzsetnextfilename{Figs/EckhardtEquilib}' ...
%,'extraAxisOptions','\extraAxisOptions' )
%%{
%\end{matlab}





\subsection{\texttt{theRes1()}: function to zero}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system, and returns the vector of patch-interior time
derivatives.
\begin{matlab}
%}
function f=theRes1(u)
  global patches 
  v=nan(size(patches.x));
  v(patches.i)=u;
  f=patchSys1(1,v(:),patches);
  f=f(patches.i);
end%function theRes1
%{
\end{matlab}


Fin.
%}