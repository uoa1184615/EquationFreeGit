% Explore heterogeneous diffusion in 1D on patches. Use the
% Jacobian to find a finite domain method that preserves
% symmetry.  AJR, Jun 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{a2Dbeam}: computational
homogenisation of a 1D heterogeneous diffusion with
Dirichlet BCs}
\label{sec:a2Dbeam}

%Finds a first method that preserves symmetry of
%heterogeneous diffusion  in 1D with Dirichlet boundary
%conditions.

First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the lattice, and random
log-normal values, albeit normalised to have harmonic mean
one.  
\begin{matlab}
%}
clear all
mPeriod = 3
cHetr=exp(0.5*randn(mPeriod,1));
cHetr = cHetr*mean(1./cHetr) % normalise
%{
\end{matlab}

Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice diffusion
system~\cref{eq:hetroDiff} solved on 
\begin{matlab}
%}
global patches xi nJac u0 i Diri bcShft bcwt
patches.intTest = true; % do not nan the two extreme ends
bcShft = 1 % relative microgrid point computed by BCs
Diri = true  % false seems worse but only ...??
edgyInt = true % symmetric in interior iff true
nSubP = (2-edgyInt)*mPeriod+1+edgyInt
ordInterp = 2
nPatch = 5+ordInterp
Dom.type = 'equispace';
Dom.bcOffset = 1
Len=nPatch-1;
dxs = 0.01*(0:7);
wts = zeros(numel(dxs),1);
for idx = 2:numel(dxs)
dx = dxs(idx)
%{
\end{matlab}
Some beam parameters, then configure
\begin{matlab}
%}
ny = 4
viscosity = 1e-3
fn = @elastic2Dstaggered
yLim= [-0.5 0.5]*dx/(ny-1)
nVars = 4*ny-2
E  = 1*exp(0.*(2*rand(mPeriod,2*ny-1)-1) );
nu =  0.3 +0.*(2*rand(mPeriod,2*ny-1)-1);
EQuartiles  = quantile(E(:) ,0:0.25:1)
nuQuartiles = quantile(nu(:),0:0.25:1)
lambda=nu.*E./(1+nu)./(1-2*nu);
mu    =    E./(1+nu)./2;
cElas = [mu lambda];
%{
\end{matlab}
Configure the patches
\begin{matlab}
%}
configPatches1(fn,[0 Len],Dom,nPatch ...
    ,ordInterp,dx,nSubP,'EdgyInt',edgyInt,'hetCoeffs',cHetr);
patches.viscosity = viscosity; 
patches.yLim = yLim;
patches.nVars = nVars;
patches.f = 0*patches.x;
%{
\end{matlab}

\paragraph{Set \(y\)-coordinates} Set the cross-beam
microscale coordinates and indices in row vectors: one for
\(u,\sigma_{xx},\sigma_{yy}\)-level; and one for
\(v,\sigma_{xy}\)-level.  Then show the micro-grid spacing.
\begin{matlab}
%}
yv=linspace(yLim(1),yLim(2),ny);
yu=(yv(1:ny-1)+yv(2:ny))/2;
ju=1:ny-1; jv=1:ny;
dy=diff(yu(1:2)), dx=diff(patches.x(1:2))
%{
\end{matlab}

%Check interpolation.
%\begin{matlab}
%}
%u0=pi*min(patches.x,Len-patches.x);
%%u0=patches.x.*(1-patches.x);
%bcwt=1
%ui=aPatchEdgeInt1(u0);
%u0=squeeze(u0)
%ui=squeeze(ui)
%return%%%%%%%%%%%%%
%{
%\end{matlab}



\paragraph{Compute Jacobian and its symmetry}
Form the Jacobian matrix, linear operator,
by numerical construction about a zero field.  Use~\verb|i|
to store the indices of the micro-grid points that are
interior to the patches and hence are the system variables.
\begin{matlab}
%}
  u0 = zeros(nSubP,1,1,nPatch);
  u0([1 end],:,:,:) = nan; u0 = u0(:);
  i = find(~isnan(u0));
  nJac = length(i)
  
bad=[];
bcwts=[1 0];
for iTry=1:3
  bcwt=bcwts(iTry)
  % the following formula is the exact weights that we find
%  if iTry==3, m=nSubP-2; bcwt=m*dx/(2-(m-1)*dx/3), end
  Jac = nan(nJac);
  for j = 1:nJac
    u0(i) = ((1:nJac)==j);
    if Diri
      u0( 1+bcShft ,:,:, 1 )=0; % left-edge of leftmost is zero
      u0(end-bcShft,:,:,end)=0; % right-edge of rightmost is zero
    else
      u0( 1+bcShft ,:,:, 1 )=u0( 2+bcShft ,:,:, 1 ); % left-edge 
      u0(end-bcShft,:,:,end)=u0(end-1-bcShft,:,:,end); % right-edge
    end
    dudt = aPatchSys1(0,u0);
    Jac(:,j) = dudt(i);
  end
  nonSymmetric = sum(sum(abs(Jac-Jac')))
  bad(iTry) = Jac(mPeriod,mPeriod+1)-Jac(mPeriod+1,mPeriod)
  if iTry==2, bcwts(3)=bad(2)/diff(bad(1:2)),  end
 end
%  Jac(abs(Jac)<1e-12) = 0;
%  Jac = Jac
spy(Jac,'o'), hold on
spy(abs(Jac-Jac')>1e-9*norm(Jac,'fro'),'rx'), hold off
%{
\end{matlab}

End of dx for-loop.   Following loop we use divided
differences to construct approximate polynomials in the dx
dependence of the weights.  Then print the rational form of
the highest order polynomial.
\begin{matlab}
%}
wts(idx)=bcwt; % store divided differences 
end%for idx
wtp=wts(1,1); % cumulatively form symbolic polynomial form 
syms d
for j=1:numel(dxs)-1
  wts(:,j+1) = [diff(wts(:,j))...
  ./(dxs(min(j+1:end+j-1,end))-dxs(1:end-1))'; nan];
  wtp = wtp+wts(1,j+1)*prod(d-dxs(1:j));
  wtcs = fliplr(sym2poly(wtp)); % store (output) coeffs
end% for j
%wts(abs(wts)>10)=inf;
[Num,Den]=rat(wtcs,1e-5); % approximate coeffs by rationals
%{
\end{matlab}
Consider the above computed numerators and denominators to
discover they come from simple rational functions, which we
can code and verify. Verify that the following coefficents
are all zero except for a one in the second position (first
power of dx).
\begin{matlab}
%}
switch mPeriod
case 2, cfs=fliplr(sym2poly(expand(wtp*(1-d/6)/1)))
case 3, cfs=fliplr(sym2poly(expand(wtp*(1-d/3)/(3/2))))
case 4, cfs=fliplr(sym2poly(expand(wtp*(1-d/2)/2)))
case 5, cfs=fliplr(sym2poly(expand(wtp*(1-2*d/5)/(5/2))))
end%switch
%{
\end{matlab}



%\input{../Patch/aHeteroDiff.m}
%}