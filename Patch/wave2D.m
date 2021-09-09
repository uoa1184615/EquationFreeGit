% Simulate the linear wave PDE in 2D on patches.
% First it checks the spectrum of the system.
% AJR, Nov 2018 -- 17 Apr 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{wave2D}: example of a wave on patches in 2D}
\label{sec:wave2D}
\localtableofcontents

For \(u(x,y,t)\), test and simulate the simple wave PDE in
2D space:
\begin{equation*}
\DD tu=\delsq u\,.
\end{equation*}
This script shows one way to get started: a user's script
may have the following three steps (left-right arrows denote
function recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches2 
\item ode15s integrator \into patchSys2 \into wavePDE
\item process results
\end{enumerate}


\begin{devMan}
Establish the global data struct \verb|patches| to interface
with a function coding the wave \pde: to be solved on
\(2\pi\)-periodic domain, with \(9\times9\) patches,
spectral interpolation~(\(0\)) couples the patches, each
patch of half-size ratio~\(0.25\) (big enough for
visualisation), and with a \(5\times5\) micro-grid within
each patch.
\begin{matlab}
%}
global patches
nSubP = 5;
nPatch = 9;
configPatches2(@wavePDE,[-pi pi], nan, nPatch, 0, 0.25, nSubP);
%{
\end{matlab}



\subsection{Check on the linear stability of the wave PDE}
Construct the systems Jacobian via numerical differentiation.
Set a zero equilibrium as basis. Then find the indices of
patch-interior points as the only ones to vary in order to
construct the Jacobian.
\begin{matlab}
%}
disp('Check linear stability of the wave scheme')
uv0 = zeros(nSubP,nSubP,2,1,nPatch,nPatch);
uv0([1 end],:,:,:,:,:) = nan;
uv0(:,[1 end],:,:,:,:) = nan;
i = find(~isnan(uv0));
%{
\end{matlab}
Now construct the Jacobian. Since this is a \emph{linear}
wave \pde, use large perturbations.
\begin{matlab}
%}
small = 1;
jac = nan(length(i));
sizeJacobian = size(jac)
for j = 1:length(i)
  uv = uv0(:);
  uv(i(j)) = uv(i(j))+small;
  tmp = patchSys2(0,uv)/small;
  jac(:,j) = tmp(i);
end
%{
\end{matlab}
Now explore the eigenvalues a little: find the ten with the
biggest real-part; if these are small enough, then the
method may be good.
\begin{matlab}
%}
evals = eig(jac);
nEvals = length(evals)
[~,k] = sort(-abs(real(evals)));
evalsWithBiggestRealPart = evals(k(1:10))
if abs(real(evals(k(1))))>1e-4
    warning('eigenvalue failure: real-part > 1e-4')
    return, end
%{
\end{matlab}
Check that the eigenvalues are close to true waves of the
\pde\ (not yet the micro-discretised equations).
\begin{matlab}
%}
kwave = 0:(nPatch-1)/2;
freq = sort(reshape(sqrt(kwave'.^2+kwave.^2),1,[]));
freq = freq(diff([-1 freq])>1e-9);
freqerr = [freq; min(abs(imag(evals)-freq))]
%{
\end{matlab}





\subsection{Execute a simulation}
Set a Gaussian initial condition using auto-replication of
the spatial grid: here \verb|u0| and~\verb|v0| are in the
form required for computation: \(n_x\times n_y\times 1\times 1\times
N_x\times N_y\).
\begin{matlab}
%}
u0 = exp(-patches.x.^2-patches.y.^2); 
v0 = zeros(size(u0));
%{
\end{matlab}
Initiate a plot of the simulation using only the microscale
values interior to the patches: set \(x\)~and \(y\)-edges to
\verb|nan| to leave the gaps. Start by showing the initial
conditions of \cref{fig:configPatches2ic} while the
simulation computes. To mesh/surf plot we need to ??
`transpose' to size \(n_x\times N_x\times n_y\times N_y\),
then reshape to size  \(n_x\cdot N_x\times n_y\cdot N_y\).
\begin{matlab}
%}
x = squeeze(patches.x); y = squeeze(patches.y);
x([1 end],:) = nan; y([1 end],:) = nan;
u = reshape(permute(squeeze(u0),[1 3 2 4]), [numel(x) numel(y)]);
usurf = surf(x(:),y(:),u');
axis([-3 3 -3 3 -0.5 1]), view(60,40)
xlabel('space x'), ylabel('space y'), zlabel('u(x,y)')
legend('time = 0','Location','north')
colormap(hsv)
drawnow
ifOurCf2eps([mfilename 'ic'])
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:wave2Dic}initial
field~\(u(x,y,t)\) at time \(t=0\) of the patch scheme
applied to the simple wave~\pde: \cref{fig:wave2Dt6} plots
the computed field at time \(t=2\).}
\includegraphics[scale=0.9]{wave2Dic}
\end{figure}
Integrate in time using standard functions.
\begin{matlab}
%}
disp('Wait while we simulate u_t=v, v_t=u_xx+u_yy')
uv0 = cat(3,u0,v0);
if ~exist('OCTAVE_VERSION','builtin')
[ts,uvs] = ode23( @patchSys2,[0 6],uv0(:));
else % octave version is slower
[ts,uvs] = odeOcts(@patchSys2,linspace(0,6),uv0(:));
end
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:wave2Dt6}.  Because of the very small time-steps,
subsample to plot at most 100 times.
\begin{matlab}
%}
di = ceil(length(ts)/100);
for i = [1:di:length(ts)-1 length(ts)]
  uv = patchEdgeInt2(uvs(i,:)); 
  uv = reshape(permute(uv,[1 5 2 6 3 4]), [numel(x) numel(y) 2]);
  set(usurf,'ZData', uv(:,:,1)');
  legend(['time = ' num2str(ts(i),2)])
  pause(0.1)
end
ifOurCf2eps([mfilename 't' num2str(ts(end))])
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:wave2Dt6}field~\(u(x,y,t)\) at time
\(t=6\) of the patch scheme applied to the simple wave~\pde\
with initial condition in \cref{fig:wave2Dic}.}
\includegraphics[scale=0.9]{wave2Dt6}
\end{figure}

\input{../Patch/wavePDE.m}
\input{../Patch/odeOcts.m}

\end{devMan}
%}