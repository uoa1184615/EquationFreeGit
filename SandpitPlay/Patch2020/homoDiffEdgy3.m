% Simulate heterogeneous diffusion in 3D on patches as an
% example application of patches in space. Here the
% microscale is of known period so we interpolate
% next-to-edge values to get opposite edge values.  
% JEB & AJR, May 2020 -- Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffEdgy3}: computational
homogenisation of a 3D diffusion via simulation on small
patches}
\label{sec:homoDiffEdgy3}



This section extends to 3D the 2D code discussed in
\cref{sec:homoDiffEdgy1}. First set random heterogeneous
diffusivities of random period in each of the two
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
mPeriod = [3 3 3]
cHetr = exp(1*randn([mPeriod 3]));
cHetr = cHetr*mean(1./cHetr(:)) 
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  Use spectral interpolation as
we test other orders subsequently.  In 3D we appear to get
only real eigenvalues by using edgy interpolation.  What
happens for non-edgy interpolation is unknown.
\begin{matlab}
%}
nPatch=mPeriod+2;
configPatches3(@heteroDiff3, [-pi pi], nan, 5 ...
    ,0, 0.3, nPatch, 'EdgyInt',true  ...
    ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Simulate}
Set initial conditions of a simulation, replicated for each
in the ensemble.
\begin{matlab}
%}
global patches
%u0 = 0.8*cos(patches.x).*sin(patches.y) ...
%     +0.1*randn([nSubP,1,1,nPatch]); 
u0 = exp(-patches.x.^2-patches.y.^2-patches.z.^2);
u0 = u0.*(0.9+0.1*rand(size(u0)));
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSmooth3, 0.5*linspace(0,1).^2, u0(:));
else % octave version
    [ts,us] = odeOcts(@patchSmooth3, 0.3*linspace(0,1).^2, u0(:));
end
%{
\end{matlab}


\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
%if ts(end)>0.099, disp('plot animation of solution field')
figure(1), clf
colormap(0.8*hsv), rgb=solarized(7);
%{
\end{matlab}
Get spatial coordinates.
\begin{matlab}
%}
x = reshape( patches.x([2:end-1],:,:,:) ,[],1); 
y = reshape( patches.y(:,[2:end-1],:,:) ,[],1);
z = reshape( patches.z(:,:,[2:end-1],:) ,[],1);
%{
\end{matlab}
For every time step draw the surface and pause for a short
display.
\begin{matlab}
%}
for i = 1:length(ts)
%{
\end{matlab}
Get the row vector of data, form into the 6D array, then omit patch edges, and reshape to suit the isosurface function.
We do not use interpolation to get edge values as the interpolation omits the corner edges (makes them nan, which often messes the isosurfaces).
\begin{matlab}
%}
  u = reshape( us(i,:) ,[mPeriod+2 nPatch]);
  u = u([2:end-1],[2:end-1],[2:end-1],:,:,:);
  u = reshape( permute(u,[1 4 2 5 3 6]) ...
      , [numel(x) numel(y) numel(z)]);
%{
\end{matlab}
Draw the isosurfaces with labels.  
Could jazz up with more colour using the following.
\begin{verbatim}
   [x,y,z,v] = flow;
   [faces,verts,colors] = isosurface(x,y,z,v,-3,x);
   patch('Vertices', verts, 'Faces', faces, ...
      'FaceVertexCData', colors, ...
      'FaceColor','interp', ...
      'edgecolor', 'interp')
\end{verbatim}
\begin{matlab}
%}
  clf;
  for p=1:2
    subplot(1,2,p)
    for iso=1:5
       isov=(1-iso/6).^2;
       hsurf(iso) = patch(isosurface(x,y,z,u,isov));  
       isonormals(x,y,z,u,hsurf(iso))
       set(hsurf(iso) ,'FaceColor',rgb(iso,:) ...
           ,'EdgeColor','none' ...
           ,'FaceAlpha',isov); 
       hold on
    end
    axis equal, view(65-7*p,35)
    axis(pi*[-1 1 -1 1 -1 1])
    xlabel('x'), ylabel('y'), zlabel('z')
    legend(['time = ' num2str(ts(i),'%4.2f')],'Location','north')
    camlight, lighting gouraud
    hold off
  end% each p
  pause(0.05)
%{
\end{matlab}
finish the animation loop and if-plot.
\begin{matlab}
%}
end%for over time
%end%if-plot
return%%%%%%%%%%%%%%%%%
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
      tmp = patchSmooth2(0,u);
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
    title(leg,'eigenvalues')
    legend boxoff 
end%if-plot
%{
\end{matlab}


\input{../Patch/heteroDiffEdgy2.m}

Fin.
%}