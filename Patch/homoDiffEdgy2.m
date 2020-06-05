%Simulate heterogeneous diffusion in 2D on patches as an
%example application of patches in space. Here the
%microscale is of known period so we interpolate
%next-to-edge values to get opposite edge values. Then
%explore the Jacobian and eigenvalues.  JEB, May 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffEdgy2}: computational homogenisation of a 2D diffusion by simulation on small patches}
\label{sec:homoDiffEdgy2}
%\localtableofcontents

This section extends the 1D code discussed in \cref{sec:homoDiffEdgy2} to 2D.

%}
clear all
mPeriod = [5 4]
% set random diffusion coefficients
rng('default');
rng(1); %
cHetr=exp(1*randn([mPeriod,2]));
cHetr = cHetr*mean(1./cHetr(:)) % normalise
nPeriodsPatch=[1 1]

global patches
nPatch = [5 5]
ratio = [0.5 0.4];
nSubP = [nPeriodsPatch(1)*mPeriod(1)+2 nPeriodsPatch(2)*mPeriod(2)+2]
patches.EdgyInt = 1; % one to use edges for interpolation
patches.EdgyEns=0; % one for ensemble of configurations
if patches.EdgyEns 
    patches.EdgyInt=1; % EdgyEns=1 implies EdgyInt=1     
   % nSubP = [3 3] % > [2 2]; when EdgyEns=1, nSubP need not depend on mPeriod
end
configPatches2(@heteroDiff2,[0 2*pi 0 2*pi],nan,nPatch ...
    ,0,ratio,nSubP);
%{
\end{matlab}

Replicate the heterogeneous coefficients across the width of
each patch. For \verb|patches.EdgyEns| an ensemble of configurations 
is constructed.
\begin{matlab}
%}
if patches.EdgyEns   
   shiftx=mod(bsxfun(@plus,0:(mPeriod(1)-1),(0:(nSubP(1)-2))'),mPeriod(1))+1;
   shifty=mod(bsxfun(@plus,0:(mPeriod(2)-1),(0:(nSubP(2)-2))'),mPeriod(2))+1;
   patches.cx=nan(nSubP(1)-1,nSubP(2)-2,mPeriod(1)*mPeriod(2));
   patches.cy=nan(nSubP(1)-2,nSubP(2)-1,mPeriod(1)*mPeriod(2));
   for p=1:mPeriod(1)
      for  q=1:mPeriod(2)
          patches.cx(:,:,(p-1)*mPeriod(2)+q)=cHetr(shiftx(:,p),shifty(1:(end-1),q),1);
          patches.cy(:,:,(p-1)*mPeriod(2)+q)=cHetr(shiftx(1:(end-1),p),shifty(:,q),2);
      end
   end
   patches.cx=permute(repmat(patches.cx,[1,1,1,nPatch]),[1,2,4,5,3]);
   patches.cy=permute(repmat(patches.cy,[1,1,1,nPatch]),[1,2,4,5,3]);
   % need to specify how configurations are coupled 
   patches.le=mod((0:(mPeriod(1)*mPeriod(2)-1))+mPeriod(2)*rem(nSubP(1)-2,mPeriod(1)),mPeriod(1)*mPeriod(2))+1;
   patches.ri=mod((0:(mPeriod(1)*mPeriod(2)-1))-mPeriod(2)*rem(nSubP(1)-2,mPeriod(1)),mPeriod(1)*mPeriod(2))+1;
   patches.bo=mod((1:(mPeriod(1)*mPeriod(2)))+nSubP(2)-3,mPeriod(2))+mPeriod(2)*floor((0:(mPeriod(1)*mPeriod(2)-1))/mPeriod(2))+1;
   patches.to=mod((1:(mPeriod(1)*mPeriod(2)))-nSubP(2)+1,mPeriod(2))+mPeriod(2)*floor((0:(mPeriod(1)*mPeriod(2)-1))/mPeriod(2))+1;
else  
   patches.cx=[repmat(cHetr(:,:,1),[(nSubP-2)./mPeriod,1]);repmat(cHetr(1,:,1),[1,(nSubP(2)-2)/mPeriod(2)])];
   patches.cy=[repmat(cHetr(:,:,2),[(nSubP-2)./mPeriod,1]),repmat(cHetr(:,1,2),[(nSubP(1)-2)/mPeriod(1),1])];
end
%{
\end{matlab}

\paragraph{Simulate}
Set the initial conditions of a simulation.
\begin{matlab}
%}
if patches.EdgyEns 
    u0=repmat(permute(reshape((cos(patches.x(:))*sin(patches.y(:))'),[nSubP(1),nPatch(1),nSubP(2),nPatch(2)]),[1,3,2,4])+0.2*randn([nSubP,nPatch]) ...
    ,1,1,1,1,mPeriod(1)*mPeriod(2));
else
   u0=permute(reshape((cos(patches.x(:))*sin(patches.y(:))'),[nSubP(1),nPatch(1),nSubP(2),nPatch(2)]),[1,3,2,4])+0.2*randn([nSubP,nPatch]);  
end
%{
\end{matlab}
Integrate using standard stiff integrators.
\begin{matlab}
%}

if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode15s(@patchSmooth2, [0 0.005 0.05 0.5], u0(:));
else % octave version
    [ts,us] = odeOcts(@patchSmooth2, [0 0.6], u0(:));
end

%{
\end{matlab}
Plot field solutions.
\begin{matlab}
%}

disp('plot sequence of surfaces')
figure(1), clf, colormap(hsv)
x = patches.x; y = patches.y;
x(end+1,:)=nan; y(end+1,:)=nan; % pad with nans
uPad=nan([nSubP+1,nPatch]);
for i = 1:length(ts)
    if  patches.EdgyEns
       uPad(1:(end-1),1:(end-1),:,:) = mean(reshape(patchEdgeInt2(us(i,:)'),[nSubP,nPatch,mPeriod(1)*mPeriod(2)]),5);
    else
        uPad(1:(end-1),1:(end-1),:,:) = reshape(patchEdgeInt2(us(i,:)'),[nSubP,nPatch]);
    end    
  u = reshape(permute(uPad,[1 3 2 4]), [numel(x) numel(y)]);
  if i==1,
  hsurf = surf(x(:),y(:),u');
  axis([0 2*pi 0 2*pi -1.3 1.3]), view(60,40) %view(60,40)
  xlabel('x'), ylabel('y'), zlabel('u')
  else   set(hsurf,'ZData', u');
  end
  legend(['time = ' num2str(ts(i),2)],'Location','north')
  caxis([0 1])
  pause(0.1)
end


%{
\end{matlab}

\paragraph{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics for a range of orders of
interpolation, all for the same patch design and
heterogeneity.  Here use a smaller ratio, and more patches,
as we do not plot.
\begin{matlab}
%}

nPatch = [10 11]
ratio = [0.1 0.1];

nSmallEvals=nPatch(1)*(nPatch(2)+1);
egs=nan(11,nSmallEvals);
ords=20;

configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
    ,0,ratio,nSubP);

if patches.EdgyEns   
   nVars=mPeriod(1)*mPeriod(2);
   shiftx=mod(bsxfun(@plus,0:(mPeriod(1)-1),(0:(nSubP(1)-2))'),mPeriod(1))+1;
   shifty=mod(bsxfun(@plus,0:(mPeriod(2)-1),(0:(nSubP(2)-2))'),mPeriod(2))+1;
   patches.cx=nan(nSubP(1)-1,nSubP(2)-2,mPeriod(1)*mPeriod(2));
   patches.cy=nan(nSubP(1)-2,nSubP(2)-1,mPeriod(1)*mPeriod(2));
   for p=1:mPeriod(1)
      for  q=1:mPeriod(2)
          patches.cx(:,:,(p-1)*mPeriod(2)+q)=cHetr(shiftx(:,p),shifty(1:(end-1),q),1);
          patches.cy(:,:,(p-1)*mPeriod(2)+q)=cHetr(shiftx(1:(end-1),p),shifty(:,q),2);
      end
   end
   patches.cx=permute(repmat(patches.cx,[1,1,1,nPatch]),[1,2,4,5,3]);
   patches.cy=permute(repmat(patches.cy,[1,1,1,nPatch]),[1,2,4,5,3]);
   % need to specify how configurations are coupled 
   patches.le=mod((0:(mPeriod(1)*mPeriod(2)-1))+mPeriod(2)*rem(nSubP(1)-2,mPeriod(1)),mPeriod(1)*mPeriod(2))+1;
   patches.ri=mod((0:(mPeriod(1)*mPeriod(2)-1))-mPeriod(2)*rem(nSubP(1)-2,mPeriod(1)),mPeriod(1)*mPeriod(2))+1;
   patches.bo=mod((1:(mPeriod(1)*mPeriod(2)))+nSubP(2)-3,mPeriod(2))+mPeriod(2)*floor((0:(mPeriod(1)*mPeriod(2)-1))/mPeriod(2))+1;
   patches.to=mod((1:(mPeriod(1)*mPeriod(2)))-nSubP(2)+1,mPeriod(2))+mPeriod(2)*floor((0:(mPeriod(1)*mPeriod(2)-1))/mPeriod(2))+1;
else  
   nVars=1; 
   patches.cx=[repmat(cHetr(:,:,1),[(nSubP-2)./mPeriod,1]);repmat(cHetr(1,:,1),[1,(nSubP(2)-2)/mPeriod(2)])];
   patches.cy=[repmat(cHetr(:,:,2),[(nSubP-2)./mPeriod,1]),repmat(cHetr(:,1,2),[(nSubP(1)-2)/mPeriod(1),1])];
end;

for ii=0:2:ords
    
ord=ii    
configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
    ,ord,ratio,nSubP);

    disp('Check linear characteristics of the patch scheme')
    u0 = zeros([nSubP, nPatch nVars]);
    u0([1 end],:,:,:,:) = nan;
    u0(:,[1 end],:,:,:) = nan;
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
    assert(notSymmetric<1e-7,'failed symmetry')
        
    
    gravEvals=[];

    disp('Find the smallest real-part eigenvalues')
    evals = eig(jac);
    biggestImag=max(abs(imag(evals)))
    nEvals = length(evals);
    [~,k] = sort(abs(real(evals)));
    evals=evals(k);
    evalsWithSmallestRealPart = evals(1:nSmallEvals);
    gravEvals=[gravEvals real(evalsWithSmallestRealPart)];

    egs(1+ord/2,:)=evalsWithSmallestRealPart;
end 


noe=20;
err=abs((egs(2:end,2:noe)-egs(1,2:noe))./egs(1,2:noe));
set(0,'defaultAxesLineStyleOrder','o-|x-|v-|s-|p-|h-|d-|^-')

plot(1:(ords/2),err(:,1:2:end))
xlabel('coupling order')
ylabel('eigenvalue relative error')
%xlim([1 10]);
set(gca, 'YScale', 'log')
leg=legend(num2str(real(egs(1,2:2:noe))','%.4f'),'Location','northeastoutside');
title(leg,'eigenvalues')
legend boxoff 



%{
\end{matlab}

\input{../Patch/heteroDiffEdgy2.m}


Fin.
%}