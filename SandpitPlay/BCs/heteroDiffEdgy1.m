% Explore heterogeneous diffusion in 1D on patches.
% explore the Jacobian and eigenvalues.  AJR, May 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{heteroDiffEdgy1}: computational
homogenisation of a 1D heterogeneous diffusion by simulation
on small patches}
\label{sec:heteroDiffEdgy1}

First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the lattice, and random
log-normal values, albeit normalised to have harmonic mean
one.  This normalisation then means that macroscale
diffusion on a domain of length~\(2\pi\) should have near
integer decay rates, the squares of \(0,1,2,\ldots\). Then
the heterogeneity is repeated to fill each patch, and
phase-shifted for an ensemble.
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
configPatches1(@heteroDiffF,[0 Len],Dom,nPatch ...
    ,ordInterp,dx,nSubP,'EdgyInt',edgyInt,'hetCoeffs',cHetr);
%{
\end{matlab}

Check interpolation.
\begin{matlab}
%}
%u0=pi*min(patches.x,Len-patches.x);
%%u0=patches.x.*(1-patches.x);
%bcwt=1
%ui=patchEdgeInt1(u0);
%u0=squeeze(u0)
%ui=squeeze(ui)
%return%%%%%%%%%%%%%
%{
\end{matlab}


Forcing to zero
\begin{matlab}
%}
global microTimePeriod
microTimePeriod = 0;
patches.f1 = 0*patches.x;
patches.f2 = 0*patches.x;
%{
\end{matlab}




\paragraph{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics for a range of orders of
interpolation, all for the same patch design and
heterogeneity.  Form the Jacobian matrix, linear operator, by numerical
construction about a zero field.  Use~\verb|i| to store the
indices of the micro-grid points that are interior to the
patches and hence are the system variables.
\begin{matlab}
%}
  u0 = zeros(nSubP,1,1,nPatch);
  u0([1 end],:,:,:) = nan; u0 = u0(:);
  i = find(~isnan(u0));
  nJac = length(i)
  
%xi=patches.x(:,1,1,1); xi=xi-mean(xi)
%X0=squeeze(mean(patches.x(:,1,1,2:end-1)))
%nonSymX0=fun(X0)
%X0=X0+0.3*randn(size(X0))
%nonSymX0=fun(X0)
%nX=numel(X0)
%A= toeplitz([1 zeros(1,nX-2)],[1 -1 zeros(1,nX-2)])
%[X,nonSymVal] = fmincon(@fun,X0,A,xi(1)*ones(nX-1,1) ...
%                       ,[],[],0*X0,Len+0*X0)  
%patches.x(:,:,:,2:end-1) = xi+shiftdim(X,-3);

bad=[];
bcwts=[1 0];
for iTry=1:3
  bcwt=bcwts(iTry)
  if iTry==3, m=nSubP-2; bcwt=m*dx/(2-(m-1)*dx/3), end
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
    dudt = patchSys1(0,u0);
    Jac(:,j) = dudt(i);
  end
%  nonnormal = sqrt(norm(Jac*Jac'-Jac'*Jac));
%  Jac(abs(Jac)<1e-12) = 0;
%  Jac = Jac
  spy(Jac,'o'), hold on
  spy(abs(Jac-Jac')>1e-9*norm(Jac,'fro'),'rx'), hold off
  nonSymmetric = sum(sum(abs(Jac-Jac')))
  bad(iTry) = Jac(mPeriod,mPeriod+1)-Jac(mPeriod+1,mPeriod)
  if iTry==2, bcwts(3)=bad(2)/diff(bad(1:2)),  end
 end
%disp("bcwt = "+num2str(bcwt,6))
%{
\end{matlab}

End of dx for-loop
\begin{matlab}
%}
wts(idx)=bcwt; % store divided differences 
wtp=wts(1,1); % cumulatively form symbolic polynomial form 
syms d
end%for idx
for j=1:numel(dxs)-1
  wts(:,j+1) = [diff(wts(:,j))...
  ./(dxs(min(j+1:end+j-1,end))-dxs(1:end-1))'; nan];
  wtp = wtp+wts(1,j+1)*prod(d-dxs(1:j));
  wtcs = fliplr(sym2poly(wtp)); % store (output) coeffs
end
%wts(abs(wts)>10)=inf;
[N,D]=rat(wtcs,1e-5); % approximate coeffs by rationals
% The following are all one in the second position (first power)
switch mPeriod
case 2, cfs=fliplr(sym2poly(expand(wtp*(1-d/6)/1)))
case 3, cfs=fliplr(sym2poly(expand(wtp*(1-d/3)/(3/2))))
case 4, cfs=fliplr(sym2poly(expand(wtp*(1-d/2)/2)))
case 5, cfs=fliplr(sym2poly(expand(wtp*(1-2*d/5)/(5/2))))
end

%{
\end{matlab}



%Function version
%\begin{matlab}
%%}
%function nonSym=fun(X)
%  global patches xi nJac u0 i Diri bcShft
%  patches.x(:,:,:,2:end-1) = xi+shiftdim(X,-3);
%%  x=squeeze(patches.x)
%  Jac = nan(nJac);
%  for j = 1:nJac
%    u0(i) = ((1:nJac)==j);
%    if Diri
%      u0( 1+bcShft ,:,:, 1 )=0; % left-edge of leftmost is zero
%      u0(end-bcShft,:,:,end)=0; % right-edge of rightmost is zero
%    else
%      u0( 1+bcShft ,:,:, 1 )=u0( 2+bcShft ,:,:, 1 ); % left-edge 
%      u0(end-bcShft,:,:,end)=u0(end-1-bcShft,:,:,end); % right-edge
%    end
%    dudt = patchSys1(0,u0);
%    Jac(:,j) = dudt(i);
%  end
%%  nonSym = norm(Jac-Jac');
%  nonSym = (Jac-Jac'); nonSym=sum(nonSym(:).^2);
%end%function
%%{
%\end{matlab}



%\input{../Patch/heteroDiff.m}
%}