% Simulate heterogeneous diffusion in 2D on patches as an
% example application of patches in space. Here the
% microscale is of known period so we interpolate
% next-to-edge values to get opposite edge values. Then
% explore the Jacobian and eigenvalues.  JEB, May 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffEdgy2}: computational homogenisation of a 2D diffusion by simulation on small patches}
\label{sec:homoDiffEdgy2}
%\localtableofcontents

This section extends the 1D code discussed in
\cref{sec:homoDiffEdgy1} to 2D.
\begin{matlab}
%}
clear all
mPeriod = 1+randi(3,1,2)
% set random diffusion coefficients
%rng('default'); rng(1); %
cHetr=exp(1*randn([mPeriod 2]));
cHetr = cHetr*mean(1./cHetr(:)) % normalisation hack
nPeriodsPatch=[1 1]

global patches
nPatch = [9 9]
ratio = [0.3 0.4];
edgyInt = true % true to use edges for interpolation
nSubP = (2-edgyInt)*nPeriodsPatch.*mPeriod+1+edgyInt
nEnsem = 1; % no. ensemble of configurations---not yet implemented??
if nEnsem>1 
    edgyInt=true; % nEnsem>1 implies EdgyInt    
   % nSubP = [3 3] % > [2 2]; when nEnsem>1, nSubP need not depend on mPeriod
end
configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
    ,0,ratio,nSubP ,'EdgyInt',edgyInt ,'nEnsem',nEnsem );
%{
\end{matlab}

Replicate the heterogeneous coefficients across the width of
each patch. First, replicate too many times, then truncate.
\begin{matlab}
%}
cx=repmat(cHetr(:,:,1),nSubP); 
cy=repmat(cHetr(:,:,2),nSubP);
patches.cx=cx(1:nSubP(1)-1,1:nSubP(2)-2);
patches.cy=cy(1:nSubP(1)-2,1:nSubP(2)-1);
%{
\end{matlab}
For $\verb|nEnsem|>1$ an ensemble of configurations 
is constructed.
\begin{matlab}
%}
if nEnsem>1
   shiftx=mod(bsxfun(@plus,0:(mPeriod(1)-1),(0:(nSubP(1)-2))'),mPeriod(1))+1;
   shifty=mod(bsxfun(@plus,0:(mPeriod(2)-1),(0:(nSubP(2)-2))'),mPeriod(2))+1;
   patches.cx=nan(nSubP(1)-1,nSubP(2)-2,mPeriod(1)*mPeriod(2));
   patches.cy=nan(nSubP(1)-2,nSubP(2)-1,mPeriod(1)*mPeriod(2));
   for p=1:mPeriod(1)
      for  q=1:mPeriod(2)
          patches.cx(:,:,(p-1)*mPeriod(2)+q) =cHetr(shiftx(:,p),shifty(1:(end-1),q),1);
          patches.cy(:,:,(p-1)*mPeriod(2)+q) =cHetr(shiftx(1:(end-1),p),shifty(:,q),2);
      end
   end
   patches.cx=permute(repmat(patches.cx,[1,1,1,nPatch]),[1,2,4,5,3]);
   patches.cy=permute(repmat(patches.cy,[1,1,1,nPatch]),[1,2,4,5,3]);
   % need to specify how configurations are coupled 
   patches.le =mod((0:(mPeriod(1)*mPeriod(2)-1))+mPeriod(2)*rem(nSubP(1)-2,mPeriod(1)),mPeriod(1)*mPeriod(2))+1;
   patches.ri =mod((0:(mPeriod(1)*mPeriod(2)-1))-mPeriod(2)*rem(nSubP(1)-2,mPeriod(1)),mPeriod(1)*mPeriod(2))+1;
   patches.bo =mod((1:(mPeriod(1)*mPeriod(2)))+nSubP(2)-3,mPeriod(2))+mPeriod(2)*floor((0:(mPeriod(1)*mPeriod(2)-1))/mPeriod(2))+1;
   patches.to =mod((1:(mPeriod(1)*mPeriod(2)))-nSubP(2)+1,mPeriod(2))+mPeriod(2)*floor((0:(mPeriod(1)*mPeriod(2)-1))/mPeriod(2))+1;
end
%{
\end{matlab}

\paragraph{Simulate}
Set initial conditions of a simulation.
\begin{matlab}
%}
u0 = cos(patches.x).*sin(patches.y) ...
     +0.2*randn([nSubP,1,1,nPatch]); 
u0 = repmat(u0,1,1,1,nEnsem,1,1); 
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSmooth2, 0.3*linspace(0,1).^2, u0(:));
else % octave version
    [ts,us] = odeOcts(@patchSmooth2, 0.3*linspace(0,1).^2, u0(:));
end
%{
\end{matlab}
Plot field solutions.
\begin{matlab}
%}
if ts(end)>0.099, disp('plot sequence of surfaces')
figure(1), clf, colormap(hsv)
x = squeeze(patches.x); y = squeeze(patches.y);
x(end+1,:)=nan;  y(end+1,:)=nan; % pad with nans
for i = 1:length(ts)
  u = squeeze( mean( patchEdgeInt2(us(i,:)) ,4));
  u(end+1,:,:)=nan; u(:,end+1,:)=nan;
  u = reshape(permute(u,[1 3 2 4]), [numel(x) numel(y)]);
  if i>1, set(hsurf,'ZData', u');
  else hsurf = surf(x(:),y(:),u'); view(60,40) 
       axis([-pi pi -pi pi -1 1])
       xlabel('x'), ylabel('y'), zlabel('u(x,y)')
  end
  legend(['time = ' num2str(ts(i),2)],'Location','north')
  caxis([-1 1])
  pause(0.1)
end%for
end%if
%{
\end{matlab}



\paragraph{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics for a range of orders of
interpolation, all for the same patch design and
heterogeneity.  Here use a smaller ratio, and more patches,
as we do not plot.
\begin{matlab}
%}
ratio = [0.1 0.1]
nSmallEvals=nPatch(1)*(nPatch(2)+1);
egs=[];     gravEvals=[];
maxords=10;


for ii=0:2:maxords
    
    ord=ii    
    configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
    ,ord,ratio,nSubP,'EdgyInt',edgyInt,'nEnsem',nEnsem);

    disp('Check linear characteristics of the patch scheme')
    u0 = zeros([nSubP,1,nEnsem,nPatch]);
    u0([1 end],:,:,:,:,:) = nan;
    u0(:,[1 end],:,:,:,:) = nan;
    i = find(~isnan(u0));

    disp('Construct the Jacobian, use large perturbations as linear')
    small = 1;
    jac = nan(length(i));
    sizeJacobian = size(jac)
    for j = 1:length(i)
      u = u0(:);
      u(i(j)) = u(i(j))+small;
      tmp = patchSmooth2(0,u)/small;
      jac(:,j) = tmp(i);
    end
    notSymmetric=norm(jac-jac')    
    if edgyInt, assert(notSymmetric<1e-7,'failed symmetry')
    elseif notSymmetric>1e-7, disp('failed symmetry')
    end 
    disp('Find the smallest real-part eigenvalues')
    if edgyInt, evals = eig((jac+jac')/2);
    else evals = eig(jac);
    end
    biggestImag=max(abs(imag(evals)))
    nEvals = length(evals);
    [~,k] = sort(abs(real(evals)));
    evals=evals(k);
    evalsWithSmallestRealPart = evals(1:nSmallEvals);
    gravEvals=[gravEvals real(evalsWithSmallestRealPart)];

    egs=[egs evalsWithSmallestRealPart];
end 

egs=egs

if maxords>2
noe=10;
err=abs(egs-egs(:,1))./(1e-7+abs(egs(:,1)));
figure(2);
semilogy(2:2:maxords,err(2:2:noe,2:end)','o:')
xlabel('coupling order')
ylabel('eigenvalue relative error')
%xlim([1 10]);
leg=legend(num2str(real(egs(2:2:noe,1)),'%.4f'),'Location','northeastoutside');
title(leg,'eigenvalues')
legend boxoff 
end


%{
\end{matlab}

\input{../Patch/heteroDiffEdgy2.m}


Fin.
%}