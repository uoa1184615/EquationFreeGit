% Finds equilibrium of forced heterogeneous diffusion in 3D
% cube on 3D patches as an example application.  Boundary
% conditions are Neumann on the right face of the cube, and
% Dirichlet on the other faces.  The microscale is of known
% period so we interpolate next-to-edge values to get
% opposite edge values.  AJR, 1 Feb 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoDiffBdryEquil3}: equilibrium via
computational homogenisation of a 3D heterogeneous diffusion
on small patches}
\label{sec:homoDiffBdryEquil3}

Find the equilibrium of a forced heterogeneous diffusion in
3D space on 3D patches as an example application. Boundary
conditions are Neumann on the right face of the cube, and
Dirichlet on the other faces.  \cref{fig:homoDiffBdryEquil3}
shows five isosurfaces of the 3D solution field.
\begin{figure}\centering
\caption{\label{fig:homoDiffBdryEquil3}% Equilibrium of the
macroscale of the random heterogeneous diffusion in 3D with
boundary conditions of zero on all faces except for the
Neumann condition on \(x=1\)
(\cref{sec:homoDiffBdryEquil3}). The small patches are
equispaced in space. }
\includegraphics[scale=0.8]{Figs/homoDiffBdryEquil3.png}
\end{figure}


Clear variables, and establish globals.
\begin{matlab}
%}
clear all
global patches
%global OurCf2eps, OurCf2eps=true %option to save plots
%{
\end{matlab}


Set random heterogeneous diffusivities of random
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
cubic domain, patches, and micro-grid spacing~\(0.05\). 
Use high order interpolation as few patches in each
direction.  Configure for Dirichlet boundaries except for
Neumann on the right \(x\)-face.
\begin{matlab}
%}
nSubP = mPeriod+2;
nPatch = 5;
Dom.type = 'equispace';
Dom.bcOffset = zeros(2,3);  Dom.bcOffset(2) = 0.5;
configPatches3(@microDiffBdry3, [-1 1], Dom ...
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
\verb|theRes|, and the number of unknowns is then its
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
information), via the generic patch system wrapper
\verb|theRes| (\cref{sec:theRes}).
\begin{matlab}
%}
disp('Solving system, takes 10--40 secs'),tic
uSoln = fsolve(@theRes,u0(patches.i) ...
        ,optimoptions('fsolve','Display','off')); 
solveTime = toc
normResidual = norm(theRes(uSoln))
normSoln = norm(uSoln)
%{
\end{matlab}
Store the solution into the patches, and give magnitudes.
\begin{matlab}
%}
u0(patches.i) = uSoln;
u0 = patchEdgeInt3(u0);
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
x = patches.x(:);  y = patches.y(:);  z = patches.z(:);
%{
\end{matlab}
Draw  isosurfaces.  Get the solution with interpolated
faces, form into a 6D array, and reshape and transpose~\(x\)
and~\(y\) to suit the isosurface function. 
\begin{matlab}
%}
u = reshape( permute(squeeze(u0),[2 5 1 4 3 6]) ...
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
Draw some isosurfaces. 
\begin{matlab}
%}
clf;
for iso=5:-1:1
   isov=(iso-0.5)/5*(maxu-minu)+minu;
   hsurf(iso) = patch(isosurface(x,y,z,u,isov));  
   isonormals(x,y,z,u,hsurf(iso))
   set(hsurf(iso) ,'FaceColor',rgb(iso,:) ...
       ,'EdgeColor','none' ,'FaceAlpha',iso/5); 
   hold on
end
hold off
axis equal, axis([-1 1 -1 1 -1 1]), view(35,25)
xlabel('$x$'), ylabel('$y$'), zlabel('$z$')
camlight, lighting gouraud
ifOurCf2eps(mfilename) %optionally save plot
if exist('OurCf2eps') && OurCf2eps,  print('-dpng',['Figs/' mfilename]), end
%{
\end{matlab}





\subsection{\texttt{microDiffBdry3()}: 3D forced
heterogeneous diffusion with boundaries}
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
%}