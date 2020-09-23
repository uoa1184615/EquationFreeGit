% Simulate heterogeneous diffusion in 1D space on 3D patches
% as an example application. Here the microscale is of known
% period so we interpolate next-to-edge values to get
% opposite edge values.  Then compute macroscale eigenvalues
% of the patch scheme applied to this heterogeneous
% diffusion to validate and to compare various orders of
% inter-patch interpolation.
% JEB & AJR, May 2020 -- Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{spmdHomoDiffEdgy31}: computational
homogenisation of a 1D diffusion via parallel simulation on small
3D patches}
\label{sec:spmdHomoDiffEdgy31}

Simulate heterogeneous diffusion along 1D space on 3D patches
as an example application. Then compute macroscale
eigenvalues of the patch scheme applied to this
heterogeneous diffusion to validate and to compare various
orders of inter-patch interpolation.

This code extends to 3D the 2D code discussed in
\cref{sec:homoDiffEdgy2}. First set random heterogeneous
diffusivities of random (small) period in each of the three
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
clear all% seem to need to clear old composites etc
mPeriod = [2 3 4]%randi([3 4],1,3)
cHetr = exp(0.3*randn([mPeriod 3]));
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
global patches %%??
nSubP=mPeriod+2;
nPatch=[7 1 1];
xLim=[-pi pi 0 0.5 0 0.5];
patches = configPatches3(@heteroDiff3, xLim, nan ...
    , nPatch, 0, [0.3 1 1], nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',cHetr ,'parallel',true );
warning('**** finished configPatches3')
%{
\end{matlab}


\subsection{Simulate heterogeneous diffusion}
Set initial conditions of a simulation as shown in
\cref{fig:spmdHomoDiffEdgy31t0}.
\begin{matlab}
%}
spmd
%warning('**** restarting spmd')
mpiprofile on
u0 = exp(-patches.x.^2/4-patches.y.^2/2-patches.z.^2);
%warning('**** finished first u0 =')
u0 = u0.*(1+0.3*rand(size(u0),patches.codist));
%warning('**** finished second u0 =')
du0dt = spmdPatchSmooth3(0,u0,patches);
warning('**** finished du0dt =')
end%spmd
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:spmdHomoDiffEdgy31t0}initial
field~\(u(x,y,z,0)\) of the patch scheme applied to a
heterogeneous diffusion~\pde.  Plotted are the isosurfaces
at field values \(u=0.1,0.3,\ldots,0.9\), with the front
quadrant omitted so you can see inside.
\cref{fig:spmdHomoDiffEdgy31tFin} plots the isosurfaces of the
computed field at time \(t=0.3\).
}
\includegraphics[width=\linewidth]{spmdHomoDiffEdgy31t0}
\end{figure}
Integrate using a distributed integrator. For initial
parameters, need micro-time steps of roughly~\(0.0001\) for
stability, recalling that, by default, \verb|spmdRK2mesoPatch|
does twenty micro-steps for each specified step in~\verb|ts|.
Use non-uniform time-steps for fun, and to show more of the
initial rapid transients. 
\begin{matlab}
%}
warning('**** Integrating system in time, wait')
ts=0.2*linspace(0,1).^2;
ts=ts(1:6)
[us,uerrs]=spmdRK2mesoPatch(ts,u0);
%warning('**** starting maxErrors')
maxError=max(uerrs{1}) % somehow uerrs is composite, but not us
%warning('**** finished and returning')
spmd,mpiprofile viewer,end
return%%%%%%%%%%%%%%%
%{
\end{matlab}


\paragraph{Plot the solution} 
First show the errors
\begin{matlab}
%}
figure(2), clf
semilogy(ts,uerrs)
xlabel('time'), ylabel('estimated step error')
%{
\end{matlab}
Then the solution field as an animation over time.
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
x = reshape(gather( patches.x([2:end-1],:,:,:) ),[],1); 
y = reshape(gather( patches.y(:,[2:end-1],:,:) ),[],1);
z = reshape(gather( patches.z(:,:,[2:end-1],:) ),[],1);
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
  u = reshape( gather(us{i}) ,[nSubP nPatch]);
  u = u([2:end-1],[2:end-1],[2:end-1],:,:,:);
  u = reshape( permute(u,[1 4 2 5 3 6]) ...
      , [numel(x) numel(y) numel(z)]);
%{
\end{matlab}
%Optionally cut-out the front corner so we can see inside.
%\begin{matlab}
%%}
%  u( (x>0) & (y'<0) & (shiftdim(z,-2)>0) ) = nan;
%%{
%\end{matlab}
The \verb|isosurface| function requires us to transpose~\(x\)
and~\(y\).
\begin{matlab}
%}
v = permute(u,[2 1 3]);
%{
\end{matlab}
Draw cross-eyed stereo view of some isosurfaces. 
\begin{matlab}
%}
  clf;
  for p=1:2
    subplot(1,2,p)
    for iso=5:-1:1
       isov=(iso-0.5)/5;
       hsurf(iso) = patch(isosurface(x,y,z,v,isov));  
       isonormals(x,y,z,v,hsurf(iso))
       set(hsurf(iso) ,'FaceColor',rgb(iso,:) ...
           ,'EdgeColor','none' ...
           ,'FaceAlpha',iso/5); 
       hold on
    end
    axis equal, axis(xLim), view(45-7*p,25)
    xlabel('x'), ylabel('y'), zlabel('z')
    legend(['time = ' num2str(ts(i),'%5.3f')],'Location','north')
    camlight, lighting gouraud
    hold off
  end% each p
  if i==1 % pause for the viewer
       makeJpeg=false;
       if makeJpeg, print(['Figs/' mfilename 't0'],'-djpeg'), end
       disp('Press any key to start animation')
       pause
  else pause(0.05)
  end
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:spmdHomoDiffEdgy31tFin}final
field~\(u(x,y,z,0.3)\) of the patch scheme applied to a
heterogeneous diffusion~\pde.  Plotted are the isosurfaces
at field values \(u=0.1,0.3,\ldots,0.9\), with the front
quadrant omitted so you can see inside.}
\includegraphics[width=\linewidth]{spmdHomoDiffEdgy31tFin}
\end{figure}
Finish the animation loop, and optionally output the
isosurfaces of the final field,
\cref{fig:spmdHomoDiffEdgy31tFin}.
\begin{matlab}
%}
end%for over time
if makeJpeg, print(['Figs/' mfilename 'tFin'],'-djpeg'), end
%{
\end{matlab}



\input{../Patch/heteroDiff3.m}

Fin.
%}