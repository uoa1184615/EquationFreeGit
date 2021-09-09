% Simulate heterogeneous waves in 2D on patches as an
% example application of patches in space. Here the
% microscale is of known period so we interpolate
% next-to-edge values to get opposite edge values.  Then
% optionally explore the Jacobian and eigenvalues.  We
% compare with "Numerical upscaling for wave equations with
% time-dependent multiscale coefficients", Bernhard Maier
% and Barbara Verfurth, arxiv.org:2107.14069 This code was
% adapted from homoDiffEdgy2.m by AJR, Aug 2021
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{homoWaveEdgy2}: computational
homogenisation of a forced, non-autonomous, 2D wave via
simulation on small patches}
\label{sec:homoWaveEdgy2}


This section extends to 2D waves, in a microscale
heterogeneous media, the 2D diffusion code discussed in
\cref{sec:homoDiffEdgy2}.  It favourably compares to the
examples of \cite{Maier2021}.

\begin{figure}
\centering \caption{\label{fig:homoWaveEdgy2fig}results for
the computational homogenisation of a forced,
non-autonomous, 2D wave (\cref{sec:homoWaveEdgy2}).
(left)~relative \textsc{rms} error of the patch scheme, each
patch of width~\(1/128\), as a function of patch
spacing~\(H\).  The unfilled symbols are those of the energy
norm from \cite{Maier2021} (their Figure~5.1).  (right)~the
relative compute time decreases very quickly in~\(H\) as
there are fewer patches spaced further apart.}
\input{../Patch/Figs/homoWaveEdgy2fig}
\end{figure}
\cref{fig:homoWaveEdgy2fig} summarises the results here. The
left (larger) graph shows the error in the patch scheme
decreasing with decreasing patch spacing~\(H\) (increasing
number of patches). Forcing~\(f_1\) and~\(f_2\) are as
specified by \S5.1 of \cite{Maier2021}, whereas~\(f_3\) here
is~\(f\) in their \S5.2. For the case of forcing~\(f_1\)
which is discontinuous in space (at \(x=0.4\)), the errors
are similar to that of \cite{Maier2021}---compare the filled
with unfilled circles. For the case of forcing~\(f_2\) which
is continuous in the spatial domain, except for a second
derivative discontinuity in its odd-periodic extension, the
errors of the patch scheme are an order of magnitude better
that that of \cite{Maier2021}---compare the filled with
unfilled squares. For the case of forcing~\(f_3\) which is
smooth in the domain and in its odd-periodic extension, the
patch scheme errors, roughly~\(10^{-8}\), are at the
tolerance of the time integration. Two caveats in a
comparison with \cite{Maier2021} are the slightly different
norms used, and that they also address errors in the time
integration, whereas here we use a standard adaptive
integrator in order to focus purely on the spatial errors of
the patch scheme.

Now let's code the simulation of the forced, non-autonomous,
2D wave. \cite{Maier2021} have Dirichlet BCs of zero around
the unit square, so replicate here by the odd periodic
extension to the spatial domain~\([-1,1]^2\).  In their
\S5.1, their microscale mesh step is~\(1/512=2^{-9}\). 
Coding that here results in a compute time of roughly
90~minutes, so here I provide a much coarser case that
computes in only a few minutes: change as you please.
\begin{matlab}
%}
clear all
dx = 1/128  % 1/512=2^{-9} is the original, but takes 90 mins
%{
\end{matlab}
The heterogeneity is of period four on the microscale
lattice, so code a minimal patch size that covers one
period.
\begin{matlab}
%}
epsilon = 4*dx
nPeriodsPatch = 1
mPeriod = round(epsilon/dx)
nSubP = mPeriod*nPeriodsPatch+2
%{
\end{matlab}
Choose which of three forcing functions to use
\begin{matlab}
%}
fn=2
%{
\end{matlab}


\cite{Maier2021} use varying number of macroscale grid steps
from~\(4\) to~\(64\) on~\([0,1]\) so here on~\([-1,1]\) we
use double the number patches in each direction.  Loop over
the number of patches used, starting with the full domain
simulation, and then progressively coarsening the macroscale
grid of patches.
\begin{matlab}
%}
nPatch = 2/epsilon/nPeriodsPatch
for iPat=0:9
if iPat>0, nPatch=nPatch/2, end
if nPatch<8, break, end
%{
\end{matlab}
Set the periodic heterogeneous coefficient, isotropic:  
\begin{equation*}
a_\epsilon(t,x)=\big[3+\sin(2\pi x/\epsilon)+\sin(2\pi t)\big]
\cdot\big[3+\sin(2\pi y/\epsilon)+\sin(2\pi t)\big],
\end{equation*}
which being in product form with two time-dependencies we
store as the two spatially varying factors---although to
preserve odd symmetry we phase shift the heterogeneity from
sines to cosines.  It is a user's choice whether to code such
spatial dependencies here with~\verb|cHetr| or within the
time derivative function itself.  In this case, I choose to
code microscale heterogeneous coefficients here
via~\verb|cHetr|, and the macroscale variation of~\(f_i\) in
the time derivative function.

Here the period of the heterogeneity is only four microscale
lattice points in each direction (which is pretty inaccurate
on the microscale, but immaterial as we and \cite{Maier2021}
only compare to the coded system on the microscale lattice,
not to the \pde).  With the following careful choices we
ensure all the hierarchy of patch schemes both maintain odd
symmetry, and also compute on grid points that are common
with the full domain. 
\begin{matlab}
%}
ratio = (nSubP-2)*dx/(2/nPatch)   
Xleft=(1-ratio)/nPatch;
xmid=Xleft+dx*(0:mPeriod-1)';    % half-points
xi = Xleft+dx*(-0.5:mPeriod-1)'; % grid-points
% two components for ax, the x-dirn interactions
cHetr(:,:,1) = (3+cos(2*pi*xmid/epsilon))+0*xi';
cHetr(:,:,2) = 0*xmid+(3+cos(2*pi*xi'/epsilon));
% two components for ay, the y-dirn interactions
cHetr(:,:,3) = (3+cos(2*pi*xi/epsilon))+0*xmid';
cHetr(:,:,4) = 0*xi+(3+cos(2*pi*xmid'/epsilon));
%{
\end{matlab}

Configure patches using spectral interpolation.  Quadratic
interpolation did not seem significantly different for the
case of discontinuous forcing~\(f_1\).
\begin{matlab}
%}
configPatches2(@heteroWave2,[-1 1 -1 1],nan,nPatch ...
    ,0,ratio,nSubP ,'EdgyInt',true ,'hetCoeffs',cHetr );
%{
\end{matlab}
A check on the spatial geometry.
\begin{matlab}
%}
global patches
dxPat=diff(patches.x(1:2));
assert(abs(dx-dxPat)<1e-9,"dx mismatch")
%{
\end{matlab}


\paragraph{Simulate}
Set the particular forcing function to use, and the zero
initial conditions of a simulation.
\begin{matlab}
%}
patches.eff=fn;
clear uv0
uv0(:,:,1,1,:,:) = 0*patches.x+0*patches.y;
uv0(:,:,2,1,:,:) = 0*patches.x+0*patches.y;
%{
\end{matlab}
Integrate using standard integrators. \cite{Maier2021} use a
scheme with fixed time-step of \(\tau=2^{-7}=1/128\)\,. 
Here \verb|ode23| uses variable steps of about~\(0.0003\)\,,
and takes~7\,s for \verb|nPatch=2*4| (whereas \verb|ode15s|
takes~149\,s---even for the dissipating case), and
takes~287\,s for \verb|nPatch=2*32| and roughly~4000\,s for
full domain \verb|nPatch=2*128|.
\begin{matlab}
%}
disp('Now simulate over time')
tic
[ts,us] = ode23(@patchSys2, linspace(0,1,11), uv0(:));
if iPat==0,  odeTime0=toc
else relOdeTime(iPat)=toc/odeTime0
end
%{
\end{matlab}


\paragraph{Compute error compared to full domain simulation}
Get spatial coordinates of patch-interior points, and
reshape to column vectors.
\begin{matlab}
%}
i = 2:nSubP-1;
x = squeeze(patches.x(i,:,:,:,:,:)); 
y = squeeze(patches.y(:,i,:,:,:,:));
x=x(:); y=y(:);
%{
\end{matlab}
At the final time of \(t=1\)\,, get the row vector of data,
form into the 6D array via the interpolation to the edges,
and reshape patch-interior points to 2D spatial array.
\begin{matlab}
%}
uv = squeeze( patchEdgeInt2(us(end,:)));
u = squeeze( uv(i,i,1,:,:) );
u = reshape(permute(u,[1 3 2 4]),[numel(x) numel(y)]);
%{
\end{matlab}
If this is the full domain simulation, then store as the
reference solution.
\begin{matlab}
%}
if iPat==0
  x0=x; y0=y; u0=u;
  rms0=sqrt(mean(u0(:).^2))
else
%{
\end{matlab}
Else compute the error compared to the full domain solution.
First find the indices of the full domain that match the
spatial locations of the patch scheme.
\begin{matlab}
%}
  [i,k] = find(abs(x0-x')<1e-9);
  assert(length(i)==length(x),'find error in index i')
  [j,k] = find(abs(y0-y')<1e-9);
  assert(length(j)==length(y),'find error in index j')
%{
\end{matlab}
The \textsc{rms} error over the surface is
\begin{matlab}
%}
  errs=u-u0(i,j);
  relrmserr(iPat)=sqrt(mean(errs(:).^2))/rms0
  H(iPat)=2/nPatch
end%if iPat
%{
\end{matlab}


End the loop over the various number of patches, and return.
Further, here not executed, code in the file animates the
solution over time, and computes spectrum of the system.
\begin{matlab}
%}
end%for iPat
figure(1), clf
loglog(H,relrmserr,'o:'),  grid on
xlabel('H'), ylabel('relative error')
return
%{
\end{matlab}
\input{../Patch/heteroWave2.m}
\endinput









\paragraph{Plot the solution} as an animation over time.
\begin{matlab}
%}
disp('plot animation of solution field')
figure(1), clf, colormap(hsv)
%{
\end{matlab}
Get spatial coordinates, pad them with NaNs to separate
patches, and reshape to vectors.
\begin{matlab}
%}
x = squeeze(patches.x); y = squeeze(patches.y);
x(end+1,:)=nan;  y(end+1,:)=nan; % pad with nans
x=x(:); y=y(:);
%{
\end{matlab}
Estimate the range of~\(u\) for plotting from the final
\(u\)-field.
\begin{matlab}
%}
uv = squeeze( patchEdgeInt2(us(end,:)));
u = squeeze( uv(:,:,1,:,:) );
maxu = max(abs(u(:))), maxu=1.05*maxu;
%{
\end{matlab}
For up to about 100 time steps, draw the surface and pause
for a short display.
\begin{matlab}
%}
di = ceil(length(ts)/100)
for i = [1:di:numel(ts)-1 numel(ts)]
%{
\end{matlab}
Get the row vector of data, form into the 6D array via the
interpolation to the edges, check odd symmetry is preserved,
then pad with Nans between patches, and reshape to suit the
surf function.
\begin{matlab}
%}
  uv = squeeze( patchEdgeInt2(us(i,:)));
  u = squeeze( uv(:,:,1,:,:) );
  v = reshape(permute(u,[1 3 2 4]),[numel(patches.x) numel(patches.y)]);
  w = v+flipud(v); udsymmetry = norm(w(~isnan(w)));
  w = v+fliplr(v); lrsymmetry = norm(w(~isnan(w)));
  assert(udsymmetry+lrsymmetry<1e-8,'not odd symmetric in space')
  u(end+1,:,:,:)=nan; u(:,end+1,:,:)=nan;
  u = reshape(permute(u,[1 3 2 4]), [numel(x) numel(y)]);
%{
\end{matlab}
If the initial time then draw the surface with labels,
otherwise just update the surface data.
\begin{matlab}
%}
  if i==1
       hsurf = surf(x,y,u'); view(60,40) 
%       axis([-1 1 -1 1 -maxu maxu]), caxis(maxu*[-1 1])
       axis([0 1 0 1 0 maxu]), caxis(maxu*[0 1])
       xlabel('x'), ylabel('y'), zlabel('u(x,y)')
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






\subsection{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics, with no forcing.
\begin{matlab}
%}
patches.eff=0;
%{
\end{matlab}
Find which elements of the 6D array are interior micro-grid
points and hence correspond to dynamical variables.
\begin{matlab}
%}
uv0([1 end],:,:) = nan;
uv0(:,[1 end],:) = nan;
i = find(~isnan(uv0));
if length(i)>4000
  disp('System too large to compute eigenvalues')
else
  disp('Computing Jacobian, and its eigenvalues')
%{
\end{matlab}
Construct the Jacobian of the scheme as the matrix of the
linear transformation, obtained by transforming the standard
unit vectors (as the system is linear).
\begin{matlab}
%}
    jac = nan(length(i));
    sizeJacobian = size(jac)
    for j = 1:length(i)
      uv = uv0(:)+(i(j)==(1:numel(uv0))');
      tmp = patchSys2(0,uv);
      jac(:,j) = tmp(i);
    end
%{
\end{matlab}
It is not symmetric, but nonetheless the eigenvalues are
good (presumably because the system represents nice wave
dynamics).
\begin{matlab}
%}
    notSymmetric=norm(jac-jac')    
%{
\end{matlab}
Find all the eigenvalues (as \verb|eigs| is unreliable).
\begin{matlab}
%}
    [evecs,evals] = eig(jac,'vector');
    biggestImag=max(abs(imag(evals)))
    biggestReal=max(abs(real(evals)))
    smallestMag=min(abs(evals))
%{
\end{matlab}
Sort eigenvalues on their real-part with most positive
first, and most negative last (ties determined by the
imaginary part). Store the leading eigenvalues in
\verb|egs|, and write out when computed all orders. The
number of zero eigenvalues, \verb|nZeroEv|, gives the number
of decoupled systems in this patch configuration.
\begin{matlab}
%}
    [~,k] = sort(-real(evals)+1e-5*abs(imag(evals)));
    evals=evals(k); evecs=evecs(:,k);
    nZeroEv=sum(abs(evals(:))<1e-4)
    leadingEvals=evals(1+4*(0:nPatch^2/1.5));
    figure(2),clf
    plot(asinh(real(evals)),asinh(imag(evals)),'.')
    xlabel('asinh(real)'), ylabel('asinh(imag)')
    pause
end%if compute eigenvals
%{
\end{matlab}
Fin.
%}