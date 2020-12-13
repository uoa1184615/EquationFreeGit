% RK2mesoPatch() is a simple example of Runge--Kutta, 2nd
% order, integration of a given deterministic system on
% patches. AJR, Sept 2020 -- Dec 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{RK2mesoPatch()}}
\label{sec:RK2mesoPatch}

This is a Runge--Kutta, 2nd order, integration of a given
deterministic system of \ode{}s on patches. It invokes
meso-time updates of the patch-edge values in order to
reduce interpolation costs, and uses a linear variation in
edge-values over the meso-time-step \cite[case
\(Q=2\)]{Bunder2015a}.  This function is aimed primarily for
large problems executed on a computer cluster to markedly
reduce expensive communication between computers.

If using within projective integration, it appears quite
tricky to get all the time-steps chosen appropriately.  One
has to choose times for: the micro-scale time-step, the
meso-time interval between communications, the longer
meso-time burst length, and the macro-scale integration
time-step.

\begin{matlab}
%}
function [xs,errs] = RK2mesoPatch(ts,x0,nMicro,patches)
if nargin<4, global patches, end
%{
\end{matlab}



\paragraph{Input}
\begin{itemize}
\item \verb|patches.fun()| is a function such as
\verb|dxdt=fun(t,x,patches)| that computes the right-hand
side of the \ode\ \(d\xv/dt=\fv(t,\xv)\) where \xv~is a
vector\slash array, \(t\)~is a scalar, and the result~\fv\
is a correspondingly sized vector\slash array.
\item \verb|x0| is an vector\slash array of initial values at
the time \verb|ts(1)|.
\item \verb|ts| is a vector of meso-scale times to compute
the approximate solution, say in~\(\RR^\ell\) for
\(\ell\geq2\)\,.
\item \verb|nMicro|, optional, default~\(10\), is the number
of micro-time-steps taken for each meso-scale time-step.
\item \verb|patches| struct set by \verb|configPatches|\(n\) and
provided as either as parameter, or as a global variable.
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item  \verb|xs|,  5/7/9D (depending upon~\verb|nD|) array
of length~\(\ell \times\cdots\) of approximate solution
vector\slash array at the specified times. But, if using
parallel computing via \verb|spmd|, then \verb|xs| is a
\emph{composite} 5/7/9D~array, so outside of an
\verb|spmd|-block access a single copy of the array
via~\verb|xs{1}|.  Similarly for~\verb|errs|.
\item \verb|errs|, column vector in \(\RR^{\ell}\) of local
error estimate for the step from~\(t_{k-1}\) to~\(t_k\).
\end{itemize}




\begin{devMan}
\paragraph{Code of RK2 integration}
Set default number of micro-scale time-steps in each
requested meso-scale step of~\verb|ts|. Cannot use
\verb|nargin| inside explicit \verb|spmd|, but can use it if
the \verb|spmd| is already active from the code that invokes
this function.
\begin{matlab}
%}
if nargin<3|isempty(nMicro), nMicro=10; end
%{
\end{matlab}
If patches are set to be in parallel (there must be a
parallel pool), but only one worker available (i.e., not
already inside \verb|spmd|), then invoke function
recursively inside \verb|spmd|. Q:~is \verb|numlabs| defined
without the parallel computing toolbox??
\begin{matlab}
%}
%warning('checking spmd status in RK2mesoPatch')
%disp(['class-patches = ',class(patches)])
if isequal(class(patches),'Composite') && numlabs==1
    spmd, %warning('recursing into RK2mesoPatch')
      [xs,errs] = RK2mesoPatch(ts,x0,nMicro,patches); 
      %warning('finished recursion into RK2mesoPatch')
    end% spmd
    assert(isequal(class(xs)  ,'Composite'),' xs  not composite')
    assert(isequal(class(errs),'Composite'),'errs not composite')
    return
end
%warning('bypassed start spmd in RK2mesoPatch')
%{
\end{matlab}

Set the number of space dimensions from the number stored
patch-size ratios.
\begin{matlab}
%}
nD = length(patches.ratio);
%{
\end{matlab}
Set the micro-time-steps and create storage for outputs.
\begin{matlab}
%}
dt = diff(ts)/nMicro;
xs = nan([numel(ts) size(x0)]);
errs = nan(numel(ts),1);
%{
\end{matlab}
Initialise first result to the given initial condition, and
evaluate the initial time derivative into~\verb|f1|. Use
inter-patch interpolation to ensure edge values of the
initial condition are defined and are reasonable.
\footnote{These \texttt{gather()} functions cause all-to-all
interprocessor communication once every meso-step.  Maybe
better to use distributed array instead, (although need to
then need to put time index last instead of first??), but we
need to do some inter-cpu communication in order to estimate
errors.}
\begin{matlab}
%}
%warning('RK2mesoPatch: x0 = patchEdgeInt3(x0)')
switch nD
  case 1, x0 = patchEdgeInt1(x0,patches);
          xs(1,:,:,:,:) = gather(x0);
  case 2, x0 = patchEdgeInt2(x0,patches);
          xs(1,:,:,:,:,:,:) = gather(x0);
  case 3, x0 = patchEdgeInt3(x0,patches);
          xs(1,:,:,:,:,:,:,:,:) = gather(x0);
  end;%switch nD
%assert(iscodistributed(x0),'x0 not codist one')
errs(1) = 0;
%warning('RK2mesoPatch: patches.fun one')
f1 = patches.fun(ts(1),x0,patches);
%{
\end{matlab}
Compute the meso-time-steps from~\(t_k\) to~\(t_{k+1}\),
copying the derivative~\verb|f1| at the end of the last
micro-time-step to be the derivative at the start of this
one.
\begin{matlab}
%}
for k = 1:numel(dt)
%{
\end{matlab}
Perform meso-time burst with the new interpolation for edge
values, and an interpolation of the time derivatives to get
derivative estimates of the edge-values.
\begin{matlab}
%}
  switch nD
    case 1, dx0 = patchEdgeInt1(f1,patches);
    case 2, dx0 = patchEdgeInt2(f1,patches);
    case 3, dx0 = patchEdgeInt3(f1,patches);
    end;%switch nD
%{
\end{matlab}
Perform the micro-time steps.
\begin{matlab}
%}
  for m=1:nMicro
    f0 = f1; 
%    assert(iscodistributed(f0),'f0 not codist')
%{
\end{matlab}
For all micro-time derivative evaluations, include that the
edge values are varying according to the estimate made at
the start of the meso-time-step.
\begin{matlab}
%}
    switch nD
      case 1, f0([1 end],:,:,:)=dx0([1 end],:,:,:);
      case 2, f0([1 end],:,:,:,:,:)=dx0([1 end],:,:,:,:,:);
              f0(:,[1 end],:,:,:,:)=dx0(:,[1 end],:,:,:,:);
      case 3
        f0([1 end],:,:,:,:,:,:,:)=dx0([1 end],:,:,:,:,:,:,:);
        f0(:,[1 end],:,:,:,:,:,:)=dx0(:,[1 end],:,:,:,:,:,:);
        f0(:,:,[1 end],:,:,:,:,:)=dx0(:,:,[1 end],:,:,:,:,:);
    end;%switch nD
%    assert(iscodistributed(f0),'f0 not codist two')
%{
\end{matlab}
Simple second-order accurate Runge--Kutta micro-scale
time-step.
\begin{matlab}
%}
    xh = x0+f0*dt(k)/2;
%    assert(iscodistributed(xh),'xh not codist')
    fh = patches.fun(ts(k)+dt(k)*(m-0.5),xh,patches);
%    assert(iscodistributed(fh),'fh not codist one')
    switch nD
      case 1, fh([1 end],:,:,:)=dx0([1 end],:,:,:);
      case 2, fh([1 end],:,:,:,:,:)=dx0([1 end],:,:,:,:,:);
              fh(:,[1 end],:,:,:,:)=dx0(:,[1 end],:,:,:,:);
      case 3
        fh([1 end],:,:,:,:,:,:,:)=dx0([1 end],:,:,:,:,:,:,:);
        fh(:,[1 end],:,:,:,:,:,:)=dx0(:,[1 end],:,:,:,:,:,:);
        fh(:,:,[1 end],:,:,:,:,:)=dx0(:,:,[1 end],:,:,:,:,:);
    end;%switch nD
%    assert(iscodistributed(fh),'fh not codist two')
    x0 = x0+fh*dt(k);
%    assert(iscodistributed(x0),'x0 not codist two')
%{
\end{matlab}
End the burst of micro-time-steps.
\begin{matlab}
%}
  end
%{
\end{matlab}
At the end of each meso-step burst, refresh the interpolate
of the edge values, evaluate time-derivative, and
temporarily fill-in edges of derivatives (to ensure error
estimate is reasonable).
\begin{matlab}
%}
  switch nD
    case 1, x0 = patchEdgeInt1(x0,patches);
            xs(k+1,:,:,:,:) = gather(x0);
    case 2, x0 = patchEdgeInt2(x0,patches);
            xs(k+1,:,:,:,:,:,:) = gather(x0);
    case 3, x0 = patchEdgeInt3(x0,patches);
            xs(k+1,:,:,:,:,:,:,:,:) = gather(x0);
    end;%switch nD
%  assert(iscodistributed(x0),'x0 not codist three')
  f1 = patches.fun(ts(k+1),x0,patches);
  switch nD
    case 1, f1([1 end],:,:,:)=dx0([1 end],:,:,:);
    case 2, f1([1 end],:,:,:,:,:)=dx0([1 end],:,:,:,:,:);
            f1(:,[1 end],:,:,:,:)=dx0(:,[1 end],:,:,:,:);
    case 3
      f1([1 end],:,:,:,:,:,:,:)=dx0([1 end],:,:,:,:,:,:,:);
      f1(:,[1 end],:,:,:,:,:,:)=dx0(:,[1 end],:,:,:,:,:,:);
      f1(:,:,[1 end],:,:,:,:,:)=dx0(:,:,[1 end],:,:,:,:,:);
    end;%switch nD
%{
\end{matlab}
Use the time derivative at~\(t_{k+1}\) to estimate an error
by storing the difference with what Simpson's rule would
estimate over the last micro-time step performed.
\begin{matlab}
%}
  f0=f0-2*fh+f1;
%  assert(iscodistributed(f0),'f2ndDeriv not codist')
  errs(k+1) = sqrt(gather(mean(f0(:).^2,'omitnan')))*dt(k)/6;
end%for-loop
end%function
%{
\end{matlab}
End of the function with results returned in~\verb|xs|
and~\verb|errs|.
\end{devMan}
%}
