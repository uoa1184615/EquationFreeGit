% Simulate heterogeneous diffusion in 2D on one patch in at
% least one direction to test the code for one patch.  Then
% explore the Jacobian and eigenvalues.  AJR, 19 Oct 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{onePatchDiff2}: computational homogenise a
2D diffusion via simulation on one/many patches}
\label{sec:onePatchDiff2}



First set random heterogeneous
diffusivities of random period in each of the two
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
clear all
mPeriod = randi([6 10],1,2)
cHetr = exp(1*randn([mPeriod 2]));
cHetr = cHetr*mean(1./cHetr(:));
facRngHetero = exp(diff(log(quantile(cHetr(:),0:1))))
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  
\begin{matlab}
%}
global patches
nPatch = randi(4,1,2)
nSubP = mPeriod+1 
% randomly select Domain types
Doms = ['periodic '; 'equispace'; 'chebyshev'];
Dom = Doms(randi(3,2,1),:) % chooses from all three
configPatches2(@heteroDiff20,[-pi pi 0 2*pi],Dom,nPatch ...
    ,4,1./(nSubP-1),nSubP ,'EdgyInt',true ,'hetCoeffs',cHetr );
%{
\end{matlab}


\paragraph{Simulate}
Set initial conditions of a simulation.
\begin{matlab}
%}
x0=min(patches.x(:)); x1=max(patches.x(:))-x0;
y0=min(patches.y(:)); y1=max(patches.y(:))-y0;
u0 = 0.7*sin((patches.x-x0)*pi/x1)+0.7*cos((patches.y-y0)*pi/y1) ...
     -0.4+0*rand([nSubP,1,1,nPatch]);
%u0int = patchEdgeInt2(u0,patches)
%du0dt = reshape(patchSys2(0,u0(:)), [nSubP,nPatch])
%return%%%%%%%%%%%%%%%%%%%
%{
\end{matlab}
Integrate using standard integrators, unevenly spaced in
time to better display transients.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSys2, min(nPatch)*0.2*linspace(0,1).^2, u0(:));
else % octave version
    [ts,us] = odeOcts(@patchSys2, 0.2*linspace(0,1).^2, u0(:));
end
%{
\end{matlab}

\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
disp('plot animation of solution field')
figure(1), clf, colormap(parula)
%{
\end{matlab}
Get spatial coordinates and pad them with NaNs to separate
patches.
\begin{matlab}
%}
x = squeeze(patches.x); y = squeeze(shiftdim(patches.y,1));
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
interpolation to the edges, apply the macroscale BCs, then pad 
with Nans between patches, and reshape to suit the surf function.
\begin{matlab}
%}
  u = patchEdgeInt2(us(i,:));
  if ~patches.periodic(1)
    u( 1 ,:,:,:, 1 ,:)=0; % left-edge of leftmost is zero
    u(end,:,:,:,end,:)=0; % right-edge of rightmost is zero
  end;
  if ~patches.periodic(2)
    u(:, 1 ,:,:,:, 1 )=u(:,  2  ,:,:,:, 1 ); % bottom-edge of bottommost
    u(:,end,:,:,:,end)=u(:,end-1,:,:,:,end); % top-edge of topmost
  end;
  u(end+1,:,:)=nan; u(:,end+1,:)=nan;
  u = reshape(permute(u,[1 5 2 6 3 4]), [numel(x) numel(y)]);
%{
\end{matlab}
If the initial time then draw the surface with labels,
otherwise just update the surface data.
\begin{matlab}
%}
  if i==1
       hsurf = surf(x(:),y(:),u'); view(60,40) 
       axis equal, zlim([-1 1])
       xlabel('$x$'), ylabel('$y$'), zlabel('$u(x,y)$')
       caxis(quantile(u(:),0:1)), colorbar
       pause
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
%{
\end{matlab}




\subsection{\texttt{heteroDiff20()}: heterogeneous diffusion}
\label{sec:heteroDiff20}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 6D input arrays~\verb|u|, \verb|x|,
and~\verb|y| (via edge-value interpolation of
\verb|patchSys2|, \cref{sec:patchSys2}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
two 2D array of diffusivities,~$c^x_{ij}$ and~$c^y_{ij}$,
have previously been stored in~\verb|patches.cs| (3D). 
\begin{matlab}
%}
function ut = heteroDiff20(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
The macroscale boundary conditions are Dirichlet zero at the
extreme edges of left-right extreme patches, and Neumann zero 
at extreme edges of top-bottom extreme patches.
\begin{matlab}
%}
  if ~patches.periodic(1)
  u( 1 ,:,:,:, 1 ,:)=0; % left-edge of leftmost is zero
  u(end,:,:,:,end,:)=0; % right-edge of rightmost is zero
  end;
  if ~patches.periodic(2)
  u(:, 1 ,:,:,:, 1 )=u(:,  2  ,:,:,:, 1 ); % bottom-edge of bottommost
  u(:,end,:,:,:,end)=u(:,end-1,:,:,:,end); % top-edge of topmost
  end;
%{
\end{matlab}
Code the microscale diffusion.
\begin{matlab}
%}
  ut(ix,iy,:,:,:,:) ...
  = diff(patches.cs(:,iy,1,:).*diff(u(:,iy,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(ix,:,2,:).*diff(u(ix,:,:,:,:,:),1,2),1,2)/dy^2; 
end% function
%{
\end{matlab}
%}

