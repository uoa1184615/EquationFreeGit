% Simulate heterogeneous diffusion in 3D on one patch in at
% least one direction to test the code for one patch width.  
% AJR, 18 Oct 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{onePatchDiff3}: computational homogenise a
3D diffusion via simulation on one/many patches}
\label{sec:onePatchDiff3}



First set random heterogeneous
diffusivities of random period in each of the two
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
clear all
mPeriod = randi([5 10],1,3)
cHetr = exp(1*randn([mPeriod 3]));
cHetr = cHetr*mean(1./cHetr(:));
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  
\begin{matlab}
%}
global patches
for i=1:9 % random choice at least one direction with one patch
    nPatch = randi(4,1,3);%2*randi(2,1,3)-1;
    if min(nPatch)==1, break, end
end%for i
nPatch = nPatch
nSubP = mPeriod+1 
% randomly select Domain types
Doms = ['periodic '; 'equispace'; 'chebyshev'];
%Dom = Doms(randi(2,3,1),:) % chooses from the first two
Dom = Doms(randi(3,3,1),:) % chooses from all three
configPatches3(@heteroDiff30,[-pi pi],Dom,nPatch ...
    ,4,1./(nSubP-1),nSubP ,'EdgyInt',true ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Simulate}
Set initial conditions of a simulation.
\begin{matlab}
%}
x0=min(patches.x(:)); x1=max(patches.x(:))-x0;
y0=min(patches.y(:)); y1=max(patches.y(:))-y0;
z0=min(patches.z(:)); z1=max(patches.z(:))-y0;
u0 = 0.5*sin((patches.x-x0)*pi/x1)+0.5*cos((patches.y-y0)*pi/y1) ...
    +0.5*cos((patches.z-z0)*2*pi/z1)-0.7+0.5*rand([nSubP,1,1,nPatch]);
%u0int = patchEdgeInt3(u0,patches)
%du0dt = reshape(patchSys3(0,u0(:)), [nSubP,nPatch])
%return%%%%%%%%%%%%%%%%%%%
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
    [ts,us] = ode23(@patchSys3, min(nPatch)*0.2*linspace(0,1,51).^2, u0(:));
%{
\end{matlab}

\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
disp('plot animation of solution field')
figure(1), clf, colormap(parula)
%{
\end{matlab}
Get spatial coordinates and pad them with NaNs to ensure patches are separated (works with function \verb|slices()|).
\begin{matlab}
%}
x = squeeze(patches.x); 
y = squeeze(shiftdim(patches.y,1));
z = squeeze(patches.z);
x(end+1,:)=nan;  y(end+1,:)=nan;  z(end+1,:)=nan; % pad with nans
[X,Y,Z]=ndgrid(x(:),y(:),z(:));
%{
\end{matlab}
For every time step draw the surface and pause for a short
display.
\begin{matlab}
%}
i = round(nSubP.*[.35;.65])
for l = 1:length(ts)
%{
\end{matlab}
Get the row vector of data,  form into the 6D array via the
interpolation to the edges, optionally apply macroscale BCs, then pad 
with Nans between patches, and reshape to suit the surf function.
\begin{matlab}
%}
  u = patchEdgeInt3(us(l,:));
  % apply BCs, copied from heteroDiff30()
  if ~patches.periodic(1)
  u( 1 ,:,:,:,:, 1 ,:,:)=0; % left-edge of leftmost is zero
  u(end,:,:,:,:,end,:,:)=0; % right-edge of rightmost is zero
  end;
  if ~patches.periodic(2)
  u(:, 1 ,:,:,:,:, 1 ,:)=0; % bottom-edge of bottommost
  u(:,end,:,:,:,:,end,:)=0; % top-edge of topmost
%  u(:, 1 ,:,:,:,:, 1 ,:)=u(:,  2  ,:,:,:,:, 1 ,:); % bottom-edge of bottommost
%  u(:,end,:,:,:,:,end,:)=u(:,end-1,:,:,:,:,end,:); % top-edge of topmost
  end;
  if ~patches.periodic(3)
  u(:,:, 1 ,:,:,:,:, 1 )=0; % back-edge of rearmost is one
  u(:,:,end,:,:,:,:,end)=1; % front-edge of frontmost is zero
  end;
  % pad with nans to ensure patches separated
  u(end+1,:,:,:)=nan; u(:,end+1,:,:)=nan; u(:,:,end+1,:)=nan;
  u = reshape(permute(u,[1 6 2 7 3 8 4 5]), [numel(x) numel(y) numel(z)]);
%{
\end{matlab}
If the initial time then draw the surface with labels,
otherwise just update the surface data. The ``+1'' in 
indices is due to the padding by NaNs.
\begin{matlab}
%}
  slices(X,Y,Z,u ...
  , i(:,1)+(nSubP(1)+1)*(0:nPatch(1)-1) ...
  , i(:,2)+(nSubP(2)+1)*(0:nPatch(2)-1) ...
  , i(:,3)+(nSubP(3)+1)*(0:nPatch(3)-1) );  
  axis equal, colorbar    
  xlabel('$x$'), ylabel('$y$'), zlabel('$z$')
  legend(['time = ' num2str(ts(l),2)],'Location','north')
  if l==1, pause, else pause(0.05), end
%{
\end{matlab}
finish the animation loop and if-plot.
\begin{matlab}
%}
end%for over time
%{
\end{matlab}



\subsection{\texttt{heteroDiff30()}: heterogeneous diffusion}
\label{sec:heteroDiff30}

This function codes the lattice heterogeneous diffusion
inside 3D patches.  For 8D input array~\verb|u| (via
edge-value interpolation of \verb|patchEdgeInt3|, such as by
\verb|patchSys3|, \cref{sec:patchSys3}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
three 3D array of diffusivities,~$c^x_{ijk}$, $c^y_{ijk}$
and~$c^z_{ijk}$, have previously been stored
in~\verb|patches.cs| (4+D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = heteroDiff30(t,u,patches)
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
If not periodic, then the macroscale boundary conditions are Dirichlet zero at the
extreme edges of left-right and front-back extreme patches, and maybe Neumann zero 
at extreme edges of top-bottom extreme patches.
\begin{matlab}
%}
  if ~patches.periodic(1)
  u( 1 ,:,:,:,:, 1 ,:,:)=0; % left-edge of leftmost is zero
  u(end,:,:,:,:,end,:,:)=0; % right-edge of rightmost is zero
  end;
  if ~patches.periodic(2)
  u(:, 1 ,:,:,:,:, 1 ,:)=0; % bottom-edge of bottommost
  u(:,end,:,:,:,:,end,:)=0; % top-edge of topmost
%  u(:, 1 ,:,:,:,:, 1 ,:)=u(:,  2  ,:,:,:,:, 1 ,:); % bottom-edge of bottommost
%  u(:,end,:,:,:,:,end,:)=u(:,end-1,:,:,:,:,end,:); % top-edge of topmost
  end;
  if ~patches.periodic(3)
  u(:,:, 1 ,:,:,:,:, 1 )=0; % back-edge of rearmost is one
  u(:,:,end,:,:,:,:,end)=1; % front-edge of frontmost is zero
  end;
%{
\end{matlab}
Reserve storage and then assign interior patch values to the
heterogeneous diffusion time derivatives. Using \verb|nan+u|
appears quicker than \verb|nan(size(u),patches.codist)|
\begin{matlab}
%}
  ut = nan+u; % reserve storage
  ut(i,j,k,:,:,:,:,:) ...
  = diff(patches.cs(:,j,k,1,:).*diff(u(:,j,k,:,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(i,:,k,2,:).*diff(u(i,:,k,:,:,:,:,:),1,2),1,2)/dy^2 ...
   +diff(patches.cs(i,j,:,3,:).*diff(u(i,j,:,:,:,:,:,:),1,3),1,3)/dz^2; 
end% function
%{
\end{matlab}
%}