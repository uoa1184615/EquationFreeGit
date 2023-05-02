% Simulate heterogeneous diffusion in 2D on patches as an
% example application of patches in space. Here the
% microscale is of known period so we interpolate
% next-to-edge values to get opposite edge values. Then
% explore the Jacobian and eigenvalues.  
% JEB & AJR, May 2020 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffEdgy2}: computational
homogenisation of a 2D diffusion via simulation on small
patches}
\label{sec:homoDiffEdgy2}



This section extends to 2D the 1D code discussed in
\cref{sec:homoDiffEdgy1}. First set random heterogeneous
diffusivities of random period in each of the two
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
mPeriod = randi([2 3],1,2)
cHetr = exp(1*randn([mPeriod 2]));
cHetr = cHetr*mean(1./cHetr(:)) 
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  Use spectral interpolation as
we test other orders subsequently.  In 2D we appear to get
only real eigenvalues by using edgy interpolation.  What
happens for non-edgy interpolation is unknown.
\begin{matlab}
%}
edgyInt = true; 
nEnsem = 1 %prod(mPeriod) % or just set one
if nEnsem==1% use more patches
    nPatch = [9 9]
    nSubP = (2-edgyInt)*mPeriod+1+edgyInt
else % when nEnsem>1 use fewer patches
    nPatch = [5 5]
    nSubP = mPeriod+randi([1 4],1,2) % +2 is decoupled
end
ratio = 0.2+0.2*rand(1,2)   
configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
    ,0,ratio,nSubP ,'EdgyInt',edgyInt ,'nEnsem',nEnsem ...
    ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Simulate}
Set initial conditions of a simulation, replicated for each
in the ensemble.
\begin{matlab}
%}
global patches
u0 = 0.8*cos(patches.x).*sin(patches.y) ...
     +0.1*randn([nSubP,1,1,nPatch]); 
u0 = repmat(u0,1,1,1,nEnsem,1,1); 
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSys2, 0.3*linspace(0,1).^2, u0(:));
else % octave version
    [ts,us] = odeOcts(@patchSys2, 0.3*linspace(0,1).^2, u0(:));
end
%{
\end{matlab}

\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
if ts(end)>0.099, disp('plot animation of solution field')
figure(1), clf, colormap(hsv)
%{
\end{matlab}
Get spatial coordinates and pad them with NaNs to separate
patches.
\begin{matlab}
%}
x = squeeze(patches.x); y = squeeze(patches.y);
x(end+1,:)=nan;  y(end+1,:)=nan; % pad with nans
%{
\end{matlab}
For every time step draw the surface and pause for a short
display.
\begin{matlab}
%}
for i = 1:length(ts)
%{
\end{matlab}
Get the row vector of data,  form into the 6D array via the
interpolation to the edges, then pad with Nans between
patches, and reshape to suit the surf function.
\begin{matlab}
%}
  u = squeeze( mean( patchEdgeInt2(us(i,:)) ,4));
  u(end+1,:,:,:)=nan; u(:,end+1,:,:)=nan;
  u = reshape(permute(u,[1 3 2 4]), [numel(x) numel(y)]);
%{
\end{matlab}
If the initial time then draw the surface with labels,
otherwise just update the surface data.
\begin{matlab}
%}
  if i==1
       hsurf = surf(x(:),y(:),u'); view(60,40) 
       axis([-pi pi -pi pi -1 1]),   caxis([-1 1])
       xlabel('$x$'), ylabel('$y$'), zlabel('$u(x,y)$')
  else set(hsurf,'ZData', u');
  end
  legend(['time = ' num2str(ts(i),2)],'Location','north')
  pause(0.05)
%{
\end{matlab}
finish the animation loop and if-plot.
\begin{matlab}
%}
end%for over time
end%if-plot
%{
\end{matlab}






\subsection{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics for a range of orders of
interpolation, all for the same patch design and
heterogeneity.  Except here use a small ratio as we do not
plot.
\begin{matlab}
%}
ratio = [0.1 0.1]
nLeadEvals=prod(nPatch)+max(nPatch);
leadingEvals=[];
%{
\end{matlab}

Evaluate eigenvalues for spectral as the base case for
polynomial interpolation of order \(2,4,\ldots\).
\begin{matlab}
%}
maxords=10;
for ord=0:2:maxords
    ord=ord    
%{
\end{matlab} 
Configure with same parameters, then because they are reset
by this configuration, restore coupling.
\begin{matlab}
%}
    configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
        ,ord,ratio,nSubP,'EdgyInt',edgyInt,'nEnsem',nEnsem ...
        ,'hetCoeffs',cHetr);
%{
\end{matlab}
Find which elements of the 6D array are interior micro-grid
points and hence correspond to dynamical variables.
\begin{matlab}
%}
    u0 = zeros([nSubP,1,nEnsem,nPatch]);
    u0([1 end],:,:) = nan;
    u0(:,[1 end],:) = nan;
    i = find(~isnan(u0));
%{
\end{matlab}
Construct the Jacobian of the scheme as the matrix of the
linear transformation, obtained by transforming the standard
unit vectors.
\begin{matlab}
%}
    jac = nan(length(i));
    sizeJacobian = size(jac)
    for j = 1:length(i)
      u = u0(:)+(i(j)==(1:numel(u0))');
      tmp = patchSys2(0,u);
      jac(:,j) = tmp(i);
    end
%{
\end{matlab}
Test for symmetry, with error if we know it should be
symmetric.
\begin{matlab}
%}
    notSymmetric=norm(jac-jac')    
    if edgyInt, assert(notSymmetric<1e-7,'failed symmetry')
    elseif notSymmetric>1e-7, disp('failed symmetry')
    end 
%{
\end{matlab}
Find all the eigenvalues (as \verb|eigs| is unreliable).
\begin{matlab}
%}
    if edgyInt, [evecs,evals] = eig((jac+jac')/2,'vector');
    else evals = eig(jac);
    end
    biggestImag=max(abs(imag(evals)));
    if biggestImag>0, biggestImag=biggestImag, end
%{
\end{matlab}
Sort eigenvalues on their real-part with most positive
first, and most negative last. Store the leading eigenvalues
in \verb|egs|, and write out when computed all orders.
The number of zero eigenvalues, \verb|nZeroEv|, gives
the number of decoupled systems in this patch configuration.
\begin{matlab}
%}
    [~,k] = sort(-real(evals));
    evals=evals(k); evecs=evecs(:,k);
    if ord==0, nZeroEv=sum(abs(evals(:))<1e-5), end
    if ord==0, evec0=evecs(:,1:nZeroEv*nLeadEvals); 
    else % find evec closest to that of each leading spectral
        [~,k]=max(abs(evecs'*evec0));
        evals=evals(k); % sort in corresponding order
    end
    leadingEvals=[leadingEvals evals(nZeroEv*(1:nLeadEvals))];
end 
disp('     spectral    quadratic      quartic  sixth-order ...')
leadingEvals=leadingEvals
%{
\end{matlab}

Plot the errors in the eigenvalues using the spectral ones
as accurate.  Only plot every second,~\verb|iEv|, as all are
repeated eigenvalues.
\begin{matlab}
%}
if maxords>2
    iEv=2:2:12;
    figure(2);
    err=abs(leadingEvals-leadingEvals(:,1)) ...
        ./(1e-7+abs(leadingEvals(:,1)));
    semilogy(2:2:maxords,err(iEv,2:end)','o:')
    xlabel('coupling order')
    ylabel('eigenvalue relative error')
    leg=legend( ...
        strcat('$',num2str(real(leadingEvals(iEv,1)),'%.4f'),'$') ...
        ,'Location','northeastoutside');
    if ~exist('OCTAVE_VERSION','builtin')
        title(leg,'eigenvalues'), end
    legend boxoff 
end%if-plot
%{
\end{matlab}


\input{../Patch/heteroDiff2.m}

Fin.
%}