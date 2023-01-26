% Finds equilibrium of forced heterogeneous diffusion in 3D
% space on 3D patches as an example application. Here the
% microscale is of known period so we interpolate
% next-to-edge values to get opposite edge values.  
% AJR, Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffBdryEquil3}: equilibrium via
computational homogenisation of a 3D diffusion on small
patches}
\label{sec:homoDiffBdryEquil3}

Find the equilibrium of a forced heterogeneous diffusion in
3D space on 3D patches as an example application. 
\begin{matlab}
%}
clear all
global patches
%{
\end{matlab}

First set random heterogeneous diffusivities of random
(small) period in each of the three directions. Crudely
normalise by the harmonic mean so the decay time scale is
roughly one. 
\begin{matlab}
%}
mPeriod = randi([2 3],1,3)
cDiff = exp(0.3*randn([mPeriod 3]));
cDiff = cDiff*mean(1./cDiff(:))
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
square domain, patches, and micro-grid spacing~\(0.05\). 
Use high order interpolation as few patches in each
direction.  Configure for Dirichlet boundaries except for
Neumann on the right \(x\)-face.
\begin{matlab}
%}
nSubP=mPeriod+2;
nPatch=5;
Dom.type='equispace';
Dom.bcOffset=zeros(2,3); Dom.bcOffset(2)=0.5;
configPatches3(@microDiffBdry3, [-1-0*rand 1+0*rand], Dom ...
    , nPatch, 0, 0.05, nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',cDiff );
%{
\end{matlab}

Set forcing, and store in global patches for access by the
microcode
\begin{matlab}
%}
patches.fu = 10*exp(-patches.x.^2-patches.y.^2-patches.z.^2);
patches.fu = patches.fu.*(1+rand(size(patches.fu)));
%{
\end{matlab}


\paragraph{Solve for steady state}
Set initial guess of zero, with \verb|NaN| to indicate
patch-edge values.  Index~\verb|i| are the indices of
patch-interior points, store in global patches for access by
\verb|theRes3|, and the number of unknowns is then its
number of elements.
\begin{matlab}
%}
u0 = zeros([nSubP,1,1,nPatch,nPatch,nPatch]);
u0([1 end],:,:,:) = nan;  
u0(:,[1 end],:,:) = nan;
u0(:,:,[1 end],:) = nan;
patches.i = find(~isnan(u0));
nVariables = numel(patches.i)
%{
\end{matlab}
Solve by iteration.  Use \verb|fsolve| for simplicity and
robustness (optionally \verb|optimoptions| to omit trace
information).
\begin{matlab}
%}
tic;
uSoln = fsolve(@theRes3,u0(patches.i));% ...
%        ,optimoptions('fsolve','Display','off')); 
solveTime = toc
%{
\end{matlab}
Store the solution into the patches, and give magnitudes.
\begin{matlab}
%}
u0(patches.i) = uSoln;
normSoln = norm(uSoln)
normResidual = norm(theRes3(uSoln))
%{
\end{matlab}



\paragraph{Plot isosurfaces of the solution}
\begin{matlab}
%}
figure(1), clf
rgb=get(gca,'defaultAxesColorOrder');
%{
\end{matlab}
Reshape spatial coordinates of patches.
\begin{matlab}
%}
x = patches.x(:); 
y = patches.y(:);
z = patches.z(:);
%{
\end{matlab}
Draw  isosurfaces.  Get the solution with interpolated
faces, form into a 6D array, and reshape and transpose~\(x\)
and~\(y\) to suit the isosurface function. 
\begin{matlab}
%}
  u = squeeze( patchEdgeInt3(u0) );
  u = reshape( permute(u,[2 5 1 4 3 6]) ...
      , [numel(y) numel(x) numel(z)]);
  maxu=max(u(:)),  minu=min(u(:))
%{
\end{matlab}
Optionally cut-out the front corner so we can see inside.
\begin{matlab}
%}
  u( (x'>0) & (y<0) & (shiftdim(z,-2)>0) ) = nan;
%{
\end{matlab}
Draw cross-eyed stereo view of some isosurfaces. 
\begin{matlab}
%}
  clf;
  for p=1:2
    subplot(1,2,p)
    for iso=5:-1:1
       isov=(iso-0.5)/5*(maxu-minu)+minu;
       hsurf(iso) = patch(isosurface(x,y,z,u,isov));  
       isonormals(x,y,z,u,hsurf(iso))
       set(hsurf(iso) ,'FaceColor',rgb(iso,:) ...
           ,'EdgeColor','none' ...
           ,'FaceAlpha',iso/5); 
       hold on
    end
    axis tight, axis equal, view(45-7*p,25)
    xlabel('x'), ylabel('y'), zlabel('z')
    camlight, lighting gouraud
    hold off
  end% each p
%{
\end{matlab}





\subsection{\texttt{microDiffBdry3()}: 3D forced heterogeneous diffusion with boundaries}
\label{sec:microDiffBdry3}

This function codes the lattice forced heterogeneous
diffusion inside the 3D patches.  For 8D input
array~\verb|u| (via edge-value interpolation of
\verb|patchEdgeInt3|, such as by \verb|patchSys3|,
\cref{sec:patchSys3}), computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.  The
three 3D array of diffusivities,~$c^x_{ijk}$, $c^y_{ijk}$
and~$c^z_{ijk}$, have previously been stored
in~\verb|patches.cs| (4D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = microDiffBdry3(t,u,patches)
  if nargin<3, global patches, end
%{
\end{matlab}
Microscale space-steps.  
%Q: is using \verb|i,j,k| slower than \verb|2:end-1|??
\begin{matlab}
%}
  dx = diff(patches.x(2:3));  % x micro-scale step
  dy = diff(patches.y(2:3));  % y micro-scale step
  dz = diff(patches.z(2:3));  % z micro-scale step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  k = 2:size(u,3)-1; % y interior points in a patch
%{
\end{matlab}
Code microscale boundary conditions of say Neumann on right,
and Dirichlet on left, top, bottom, front, and back (viewed
along the \(z\)-axis).
\begin{matlab}
%}
u( 1 ,:,:,:,:, 1 ,:,:)=0; %left face of leftmost patch
u(end,:,:,:,:,end,:,:)=u(end-1,:,:,:,:,end,:,:); %right face of rightmost
u(:, 1 ,:,:,:,:, 1 ,:)=0; %bottom face of bottommost 
u(:,end,:,:,:,:,end,:)=0; %top face of topmost
u(:,:, 1 ,:,:,:,:, 1 )=0; %front face of frontmost 
u(:,:,end,:,:,:,:,end)=0; %back face of backmost
%{
\end{matlab}
Reserve storage and then assign interior patch values to the
heterogeneous diffusion time derivatives. Using \verb|nan+u|
appears quicker than \verb|nan(size(u),patches.codist)|
\begin{matlab}
%}
  ut = nan+u; % reserve storage
  ut(i,j,k,:) ...
  = diff(patches.cs(:,j,k,1).*diff(u(:,j,k,:),1,1),1,1)/dx^2 ...
   +diff(patches.cs(i,:,k,2).*diff(u(i,:,k,:),1,2),1,2)/dy^2 ...
   +diff(patches.cs(i,j,:,3).*diff(u(i,j,:,:),1,3),1,3)/dz^2 ...
   +patches.fu(i,j,k); 
end% function
%{
\end{matlab}




\subsection{\texttt{theRes3()}: function to zero}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system, and returns the vector of patch-interior time
derivatives.
\begin{matlab}
%}
function f=theRes3(u)
  global patches 
  v=nan(size(patches.x+patches.y+patches.z));
  v(patches.i)=u;
  f=patchSys3(0,v(:),patches);
  f=f(patches.i);
end%function theRes
%{
\end{matlab}


Fin.
%}