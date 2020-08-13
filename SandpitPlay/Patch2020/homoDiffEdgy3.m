% Simulate heterogeneous diffusion in 3D space on 3D patches
% as an example application. Here the microscale is of known
% period so we interpolate next-to-edge values to get
% opposite edge values.  
% JEB & AJR, May 2020 -- Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffEdgy3}: computational
homogenisation of a 3D diffusion via simulation on small
patches}
\label{sec:homoDiffEdgy3}



This section extends to 3D the 2D code discussed in
\cref{sec:homoDiffEdgy2}. First set random heterogeneous
diffusivities of random period in each of the three
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
mPeriod = randi([2 3],1,3)
cHetr = exp(0.3*randn([mPeriod 3]));
cHetr = cHetr*mean(1./cHetr(:)); 
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  Use spectral interpolation as
we test other orders subsequently.  In 3D we appear to get
only real eigenvalues by using edgy interpolation.  What
happens for non-edgy interpolation is unknown.
\begin{matlab}
%}
nSubP=mPeriod+2;
nPatch=[5 5 5];
configPatches3(@heteroDiff3, [-pi pi], nan, nPatch ...
    ,0, 0.3, nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Simulate}
Set initial conditions of a simulation, replicated for each
in the ensemble.
\begin{matlab}
%}
global patches
u0 = exp(-patches.x.^2/4-patches.y.^2/2-patches.z.^2);
u0 = u0.*(1+0.3*rand(size(u0)));
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSmooth3, 0.25*linspace(0,1,50).^2, u0(:));
else % octave version
    [ts,us] = odeOcts(@patchSmooth3, 0.25*linspace(0,1).^2, u0(:));
end
%{
\end{matlab}


\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
figure(1), clf
rgb=get(gca,'defaultAxesColorOrder');
colormap(0.8*hsv)
%{
\end{matlab}
Get spatial coordinates of patch interiors.
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
Get the row vector of data, form into a 6D array, then omit
patch faces, and reshape to suit the isosurface function. We
do not use interpolation to get face values as the
interpolation omits the corner edges and so breaks up the
isosurfaces.
\begin{matlab}
%}
  u = reshape( us(i,:) ,[nSubP nPatch]);
  u = u([2:end-1],[2:end-1],[2:end-1],:,:,:);
  u = reshape( permute(u,[1 4 2 5 3 6]) ...
      , [numel(x) numel(y) numel(z)]);
%{
\end{matlab}
Draw cross-eyed stereo view of some isosurfaces.  
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
       isov=(1-(iso-0.5)/5);
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
  if i==1, disp('Press any key to start'), pause
  else pause(0.05)
  end
%{
\end{matlab}
Finish the animation loop.
\begin{matlab}
%}
end%for over time
%{
\end{matlab}






\subsection{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics for a range of orders of
interpolation, all for the same random patch design and
heterogeneity.  Except here use a small ratio as we do not
plot and then the scale separation is clearest.
\begin{matlab}
%}
ratio = 0.025*(1+rand(1,3))
nSubP=randi([3 5],1,3)
nPatch=[3 3 3]
nEnsem = prod(mPeriod) % or just set one
%{
\end{matlab}
Find which elements of the 8D array are interior micro-grid
points and hence correspond to dynamical variables.
\begin{matlab}
%}
    u0 = zeros([nSubP,1,nEnsem,nPatch]);
    u0([1 end],:,:,:) = nan;
    u0(:,[1 end],:,:) = nan;
    u0(:,:,[1 end],:) = nan;
    i = find(~isnan(u0));
    sizeJacobian = length(i)
    if sizeJacobian>4000,
        warning('Jacobian is too big to generate and analyse')
        return
    end
%{
\end{matlab}
Store this many eigenvalues in array across different orders
of interpolation.
\begin{matlab}
%}
nLeadEvals=prod(nPatch)+max(nPatch);
leadingEvals=[];
%{
\end{matlab}

Evaluate eigenvalues for spectral as the base case for
polynomial interpolation of order \(2,4,\ldots\).
\begin{matlab}
%}
maxords=6;
for ord=0:2:maxords
    ord=ord    
%{
\end{matlab} 
Configure with same heterogeneity.
\begin{matlab}
%}
    configPatches3(@heteroDiff3,[-pi pi],nan,nPatch ...
        ,ord,ratio,nSubP,'EdgyInt',true,'nEnsem',nEnsem ...
        ,'hetCoeffs',cHetr);
%{
\end{matlab}
Construct the Jacobian of the scheme as the matrix of the
linear transformation, obtained by transforming the standard
unit vectors.
\begin{matlab}
%}
    jac = nan(length(i));
    for j = 1:length(i)
      u = u0(:)+(i(j)==(1:numel(u0))');
      tmp = patchSmooth3(0,u);
      jac(:,j) = tmp(i);
    end
%{
\end{matlab}
Test for symmetry, with error if we know it should be
symmetric.
\begin{matlab}
%}
    notSymmetric=norm(jac-jac')   
    if notSymmetric>1e-7, spy(abs(jac-jac')>1e-7), end%??
    assert(notSymmetric<1e-7,'failed symmetry')
    pause(1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%??
%{
\end{matlab}
Find all the eigenvalues (as \verb|eigs| is unreliable).
\begin{matlab}
%}
    [evecs,evals] = eig((jac+jac')/2,'vector');
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


\input{../Patch/heteroDiffEdgy3.m}

Fin.
%}