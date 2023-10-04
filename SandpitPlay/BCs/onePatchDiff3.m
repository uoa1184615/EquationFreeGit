% Simulate heterogeneous diffusion in 3D on one patch in at
% least one direction to test the code for one patch width.  
% AJR, 4 Oct 2023
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
nPatch = 2*randi(2,1,3)-1
nSubP = mPeriod+1 
if 1, Dom='equispace', else Dom='chebyshev', end
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
    +0.5*cos((patches.z-z0)*2*pi/z1)-0.5+0*rand([nSubP,1,1,nPatch]);
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
if ts(end)>0.099, disp('plot animation of solution field')
figure(1), clf, colormap(parula)
%{
\end{matlab}
Get spatial coordinates and pad them with NaNs to separate
patches (works with function \verb|slices()|).
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
  if 1 % apply BCs because slice chokes on NaNs??
    u( 1 ,:,:,:,:, 1 ,:,:)=0; % left-edge of leftmost is zero
    u(end,:,:,:,:,end,:,:)=0; % right-edge of rightmost is zero
    u(:, 1 ,:,:,:,:, 1 ,:)=u(:,  2  ,:,:,:,:, 1 ,:); % bottom-edge of bottommost
    u(:,end,:,:,:,:,end,:)=u(:,end-1,:,:,:,:,end,:); % top-edge of topmost
    u(:,:, 1 ,:,:,:,:, 1 )=0; % back-edge of rearmost is one
    u(:,:,end,:,:,:,:,end)=1; % front-edge of frontmost is zero
  end%if BC option
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
  axis equal, colorbar     %, clim([-1 1])  %view(60,40)
  xlabel('$x$'), ylabel('$y$'), zlabel('$z$')
  legend(['time = ' num2str(ts(l),2)],'Location','north')
  if l==1, pause, else pause(0.05), end
%{
\end{matlab}
finish the animation loop and if-plot.
\begin{matlab}
%}
end%for over time
end%if-plot
%{
\end{matlab}


\input{../Patch/heteroDiff30.m}
%}