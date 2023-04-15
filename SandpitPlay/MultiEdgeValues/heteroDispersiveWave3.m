% Simulate in 3D on patches the heterogeneous dispersive
% waves in a fourth-order wave PDE.  AJR, 16 Apr 2023 
%!TEX root = doc.tex
%{
\section{\texttt{heteroDispersiveWave3}: heterogeneous
Dispersive Waves from 4th order PDE}
\label{sec:heteroDispersiveWave3}

This uses small spatial patches to simulate heterogeneous
dispersive waves in 3D.  The wave equation 
for~\(u(x,y,z,t)\) is the fourth-order in space \pde
\begin{equation*}
u_{tt}=-\delsq(C\delsq u)
\end{equation*}
for microscale variations in scalar~\(C(x,y,z)\).


Initialise some Matlab aspects.
\begin{matlab}
%}
clear all
cMap=jet(64); cMap=0.8*cMap(7:end-7,:); % set colormap
basename = [num2str(floor(1e5*rem(now,1))) mfilename]
%global OurCf2eps, OurCf2eps=true %optional to save plots
%{
\end{matlab}


Set random heterogeneous coefficients of period two in each
of the three directions. Crudely normalise by the harmonic
mean so the macro-wave time scale is roughly one. 
\begin{matlab}
%}
mPeriod = [2 2 2];
cHetr = exp(0.9*randn(mPeriod));
cHetr = cHetr*mean(1./cHetr(:)) 
%{
\end{matlab}

Establish global patch data struct to interface with a
function coding a fourth-order heterogeneous wave \pde: to
be solved on $[-\pi,\pi]^3$-periodic domain, with
$5^3$~patches, spectral interpolation~($0$) couples the
patches, each patch with micro-grid spacing~$0.22$
(relatively large for visualisation), and with $6^3$~points
forming each patch.  (Six because two edge layers on each of
two faces, and two interior points for the \pde.)
\begin{matlab}
%}
global patches
patches = configPatches3(@heteroDispWave3,[-pi pi] ...
    ,'periodic', 5, 0, 0.22, mPeriod+4 ,'EdgyInt',true ...
    ,'hetCoeffs',cHetr ,'nEdge',2);
%{
\end{matlab}
Set a wave initial state using auto-replication of the
spatial grid, and as \cref{fig:heteroDispersiveWave3ic}
shows. This wave propagates diagonally across space.
Concatenate the two \(u,v\)-fields to be the two components
of the fourth
dimension.
\begin{matlab}
%}
u0 = 0.5+0.5*sin(patches.x+patches.y+patches.z);
v0 =    -0.5*cos(patches.x+patches.y+patches.z)*3;
uv0 = cat(4,u0,v0);
%{
\end{matlab}
\begin{figure}\centering
\caption{\label{fig:heteroDispersiveWave3ic} initial
field~$u(x,y,z,t)$ at time $t=0$ of the patch scheme applied
to a heterogeneous dispersive wave~\pde:
\cref{fig:heteroDispersiveWave3fin} plots the computed field
at time $t=6$.}
\includegraphics[scale=0.9]{24168heteroDispersiveWave3ic}
\end{figure}
Integrate in time to $t=6$ using standard functions. In
Matlab \verb|ode15s| would be natural as the patch scheme is
naturally stiff, but \verb|ode23| is much quicker
\cite[Fig.~4]{Maclean2020a}.
\begin{matlab}
%}
disp('Simulate heterogeneous wave u_tt=delsq[C*delsq(u)]')
tic
[ts,us] = ode23(@patchSys3,linspace(0,6),uv0(:));
simulateTime=toc
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:heteroDispersiveWave3fin}.  Use
\verb|patchEdgeInt3| to obtain patch-face values in order to
most easily reconstruct the array data structure.

Replicate $x$, $y$, and~$z$ arrays to get individual
spatial coordinates of every data point.  Then, optionally,
set faces to~\verb|nan| so the plot just shows
patch-interior data. 
\begin{matlab}
%}
%%
figure(1), clf, colormap(cMap)
xs = patches.x+0*patches.y+0*patches.z;
ys = patches.y+0*patches.x+0*patches.z;
zs = patches.z+0*patches.y+0*patches.x;
if 1, xs([1:2 end-1:end],:,:,:)=nan;  
      xs(:,[1:2 end-1:end],:,:)=nan;  
      xs(:,:,[1:2 end-1:end],:)=nan;  
end;%option
j=find(~isnan(xs));
%{
\end{matlab}
In the scatter plot, \verb|col()| maps the $u$-data values
to the colour of the dots.
\begin{matlab}
%}
col = @(u) sign(u).*abs(u);
%{
\end{matlab}
Loop to plot at each and every time step.
\begin{matlab}
%}
for i = 1:length(ts)
  uv = patchEdgeInt3(us(i,:));
  u = uv(:,:,:,1,:);
  for p=1:2
    subplot(1,2,p)
    if (i==1)
      scat(p) = scatter3(xs(j),ys(j),zs(j),'.'); 
      axis equal, caxis(col([0 1])), view(45-4*p,42)
      xlabel('x'), ylabel('y'), zlabel('z')
    end 
    title(['view cross-eyed:  time = ' num2str(ts(i),'%4.2f')])
    set( scat(p),'CData',col(u(j)) );
  end
%{
\end{matlab}
Optionally save the initial condition to graphic file for
\cref{fig:heteroDispersiveWave3ic}, and optionally save the
last plot.
\begin{matlab}
%}
  if i==1,
    ifOurCf2eps([basename 'ic'])
    disp('Type space character to animate simulation')
    pause
  else pause(0.1)
  end
end% i-loop over all times
ifOurCf2eps([basename 'fin'])
%{
\end{matlab}
\begin{figure}\centering
\caption{\label{fig:heteroDispersiveWave3fin}
field~$u(x,y,z,t)$ at time $t=6$ of the patch scheme applied
to the heterogeneous dispersive wave~\pde\ with initial
condition in \cref{fig:heteroDispersiveWave3ic}.}
\includegraphics[scale=0.9]{24168heteroDispersiveWave3fin}
\end{figure}





\subsection{\texttt{heteroDispWave3()}: PDE function of 
4th-order heterogeneous dispersive waves}
\label{sec:heteroDispWave3}

This function codes the lattice heterogeneous waves inside
the patches.  The wave \pde\  for \(u(x,y,z,t)\) and
`velocity'~\(v(x,y,z,t)\) is
\begin{equation*}
u_t=v,\quad v_t=-\delsq(C\delsq u)
\end{equation*}
for microscale variations in scalar~\(C(x,y,z)\). For 8D
input arrays~\verb|u|, \verb|x|, \verb|y|, and~\verb|z| (via
edge-value interpolation of \verb|patchSys3|,
\cref{sec:patchSys3}), computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.  The
3D array of heterogeneous coefficients,~$C_{ijk}$,
$c^y_{ijk}$ and~$c^z_{ijk}$, have been stored
in~\verb|patches.cs| (3D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = heteroDispWave3(t,u,patches)
  if nargin<3, global patches, end
%{
\end{matlab}
Micro-grid space steps.  
\begin{matlab}
%}
dx = diff(patches.x(2:3));  
dy = diff(patches.y(2:3));  
dz = diff(patches.z(2:3));  
%{
\end{matlab}
First, compute \(C\delsq u\) into say~\verb|u|, using
indices for all but extreme micro-grid points.  We use a
single colon to represent the last four array dimensions
because the result arrays are already dimensioned.
\begin{matlab}
%}  
I = 2:size(u,1)-1; J = 2:size(u,2)-1;  K = 2:size(u,3)-1;
u(I,J,K,1,:) = patches.cs(I,J,K,1,:).*( diff(u(:,J,K,1,:),2,1)/dx^2 ...
    +diff(u(I,:,K,1,:),2,2)/dy^2 +diff(u(I,J,:,1,:),2,3)/dz^2 ); 
%{
\end{matlab}
Reserve storage, set lowercase indices to non-edge interior,
and then assign interior patch values to the heterogeneous
diffusion time derivatives. 
\begin{matlab}
%}
ut = nan+u;  % preallocate output array
i = I(2:end-1); j = J(2:end-1);  k = K(2:end-1);
ut(i,j,k,1,:) = u(i,j,k,2,:);  % du/dt=v
% dv/dt=delta^2 of above C*delta^2
ut(i,j,k,2,:) = -( diff(u(I,j,k,1,:),2,1)/dx^2 ...
  +diff(u(i,J,k,1,:),2,2)/dy^2 +diff(u(i,j,K,1,:),2,3)/dz^2 ); 
end% function
%{
\end{matlab}
%}
