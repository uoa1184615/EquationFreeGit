% spmdHomoDiff31 simulates heterogeneous diffusion in 1D
% space on 3D patches as a Proof of Principle example of
% parallel computing with spmd.  The interest here is on
% using spmd and comparing it with code not using spmd. The
% discussion here only addresses issues with spmd and
% parallel computing. For discussion on the 3D patch scheme
% with heterogeneous diffusion, see code and documentation
% for homoDiffEdgy3.  AJR, Aug--Sep 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{spmdHomoDiff31}: computational
homogenisation of a 1D dispersion via parallel simulation on
small 3D patches of diffusion}
\label{sec:spmdHomoDiff31}
\localtableofcontents

Simulate heterogeneous dispersion along 1D space on 3D
patches as a Proof of Principle example of parallel
computing with \verb|spmd|.  The discussion here only
addresses issues with \verb|spmd| parallel computing.  For
discussion on the 3D patch scheme with heterogeneous
diffusion, see code and documentation for
\verb|homoDiffEdgy3| in \cref{sec:homoDiffEdgy3}.


Choose one of four cases:
\begin{itemize}
\item \verb|theCase=1| is corresponding code without
parallelisation (in this toy problem it is much the quickest
because there is no expensive communication);
\item \verb|theCase=2| for minimising coding by a user of
\verb|spmd|-blocks;
\item \verb|theCase=3| is for users happier to explicitly
invoke \verb|spmd|-blocks.
\item \verb|theCase=4| invokes projective integration for
long-time simulation via short bursts of the
micro-computation, bursts done within \verb|spmd|-blocks for
parallel computing.
\end{itemize}
First, clear all to remove any existing globals, old
composites, etc---although a parallel pool persists.
Then choose the case.
\begin{matlab}
%}
clear all
theCase = 2
%{
\end{matlab}

Set micro-scale heterogeneity with various periods.
\begin{matlab}
%}
mPeriod = [4 3 2] %1+randperm(3) 
cHetr = exp(0.3*randn([mPeriod 3]));
cHetr = cHetr*mean(1./cHetr(:))
%{
\end{matlab}

Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios---here each patch is a unit
cube in space.   Set \verb|patches| information to be global
so the info can be used for Case~1 without being explicitly
passed as arguments.  Choose the parallel option if not
Case~1, which invokes \verb|spmd|-block internally, so that
field variables become \emph{distributed} across cpus.
\begin{matlab}
%}
if any(theCase==[1 2]), global patches, end
nSubP=mPeriod+2
nPatch=[9 1 1]
ratio=0.3
xLim=[0 nPatch(1)/ratio 0 1 0 1] 
disp('**** Setting configPatches3')
patches = configPatches3(@heteroDiff3, xLim, nan ...
    , nPatch, 0, [ratio 1 1], nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',cHetr ,'parallel',(theCase>1) );
%{
\end{matlab}


\subsection{Simulate heterogeneous diffusion}
Set initial conditions of a simulation as shown in
\cref{fig:spmdHomoDiff31t0}.
\begin{matlab}
%}
disp('**** Set initial condition and testing du0dt =')
if theCase==1
%{
\end{matlab}
Without parallel processing, invoke the usual operations.
\begin{matlab}
%}
    u0 = exp( -(patches.x-xLim(2)/2).^2/xLim(2) ...
               -patches.y.^2/2-patches.z.^2 );
    u0 = u0.*(1+0.2*rand(size(u0)));
    du0dt = patchSmooth3(0,u0);
%{
\end{matlab}
With parallel, must use an \verb|spmd|-block for
computations: there is no difference in cases~2--4 here.
Also, we must sometimes explicitly code how to distribute
some new arrays over the cpus. Now \verb|patchSmooth3| does
not invoke \verb|spmd| so higher level code must, as here.
Even if \verb|patches| is global, inside \verb|spmd|-block
we must pass it explicitly as a parameter to
\verb|patchSmooth3|.
\begin{matlab}
%}
else, spmd
    u0 = exp( -(patches.x-xLim(2)/2).^2/xLim(2) ...
               -patches.y.^2/2-patches.z.^2/4 );
    u0 = u0.*(1+0.2*rand(size(u0),patches.codist));
    du0dt = patchSmooth3(0,u0,patches);
    end%spmd
end%if theCase
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:spmdHomoDiff31t0}initial
field~\(u(x,y,z,0)\) of the patch scheme applied to a
heterogeneous diffusion~\pde. The vertical spread indicates 
the extent of the structure in~\(u\) in the cross-section 
variables~\(y,z\).  \cref{fig:spmdHomoDiff31tFin}
plots the nearly smooth field values at time \(t=0.4\). }
\includegraphics[scale=0.8]{spmdHomoDiff31t0}
\end{figure}

Integrate in time.  Use non-uniform time-steps for
fun, and to show more of the initial rapid transients. 

Alternatively, use \verb|RK2mesoPatch| which reduces
communication between patches, recalling that, by default,
\verb|RK2mesoPatch| does ten micro-steps for each specified
step in~\verb|ts|. For unit cube patches, need micro-steps
less than about~\(0.004\) for stability.
\begin{matlab}
%}
warning('Integrating system in time, wait patiently')
ts=0.4*linspace(0,1,21).^2;
%{
\end{matlab}
Go to the selected case.
\begin{matlab}
%}
switch theCase
%{
\end{matlab}
\begin{enumerate}
\item For non-parallel, we could use \verb|RK2mesoPatch| as
indicated below, but instead choose to use standard
\verb|ode23| as here \verb|patchSmooth3| accesses patch
information via global \verb|patches|. For post-processing,
reshape each and every row of the computed solution to the
correct array size---that of the initial condition. 
\begin{matlab}
%}
case 1
%    [us,uerrs] = RK2mesoPatch(ts,u0);
    [ts,us] = ode23(@patchSmooth3,ts,u0(:));
    us=reshape(us,[length(ts) size(u0)]);
%{
\end{matlab}

\item In the second case, \verb|RK2mesoPatch| detects a
parallel patch code has been requested, but has only one cpu
worker, so it auto-initiates an \verb|spmd|-block for the
integration. Both this and the next case return
\emph{composite} results, so just keep one version of the
results.
\begin{matlab}
%}
case 2
    us = RK2mesoPatch(ts,u0);
    us = us{1};  
%{
\end{matlab}

\item In this third case, a user could merge this explicit
\verb|spmd|-block with the previous one that sets the
initial conditions.
\begin{matlab}
%}
case 3,spmd
    us = RK2mesoPatch(ts,u0,[],patches);
    end%spmd
    us = us{1};
%{
\end{matlab}

\item In this fourth case, use Projective Integration (PI)
over long times (\verb|PIRK4| also works). Currently the PI
is done serially, with parallel \verb|spmd|-blocks only
invoked inside function \verb|aBurst()| (\cref{secmBfPI}) to
compute each burst of the micro-scale simulation. A
macro-scale time-step of about~\(3\) seems good to resolve
the decay of the macro-scale `homogenised' diffusion.
\footnote{Curiously, \texttt{PIG()} appears to suffer
unrecoverable instabilities with its variable step size!}
The function \verb|microBurst()| here interfaces to
\verb|aBurst()| (\cref{secmBfPI}) in order to provide shaped
initial states, and to provide the patch information.
\begin{matlab}
%}
case 4
    microBurst = @(tb0,xb0,bT) ...
        aBurst(tb0 ,reshape(xb0,size(u0)) ,patches);
    ts = 0:3:51
    us = PIRK2(microBurst,ts,gather(u0(:)));
    us = reshape(us,[length(ts) size(u0)]);
%{
\end{matlab}
\end{enumerate}
End the four cases.
\begin{matlab}
%}
end%switch theCase
%{
\end{matlab}




\subsection{Plot the solution} 
Optionally save some plots to file.
\begin{matlab}
%}
if 0, global OurCf2eps, OurCf2eps=true, end
%{
\end{matlab}
Animate the solution field over time. Since the spatial
domain is long in~\(x\) and thin in~\(y,z\), just plot field
values as a function of~\(x\).
\begin{matlab}
%}
figure(1), clf
if theCase==1
    x = reshape( patches.x(2:end-1,:,:,:) ,[],1); 
else, spmd
    x = reshape(gather( patches.x(2:end-1,:,:,:) ),[],1); 
    end%spmd
    x = x{1};
end
%{
\end{matlab}
For every time step draw the field values as dots and pause
for a short display.
\begin{matlab}
%}
nTimes = length(ts)
for l = 1:length(ts)
%{
\end{matlab}
At each time, squeeze interior point data into a 4D array,
permute to get all the \(x\)-variation in the first two
dimensions, and reshape into \(x\)-variation for each and
every~\((y,z)\). 
\begin{matlab}
%}
  u = reshape( permute( squeeze( ...
      us(l,2:end-1,2:end-1,2:end-1,:) ) ,[1 4 2 3]) ,numel(x),[]);
%{
\end{matlab}
Draw point data to show spread at each cross-section, as
well as macro-scale variation in the long space direction. 
\begin{matlab}
%}
  if l==1
    hp = plot(x,u,'.');
    axis([xLim(1:2) 0 max(u(:))])
    xlabel('space x'), ylabel('u(x,y,z,t)')
    ifOurCf2eps([mfilename 't0'])
    legend(['time = ' num2str(ts(l),'%4.2f')])
    disp('**** pausing, press blank to animate')
    pause
 else
    for p=1:size(u,2), hp(p).YData=u(:,p); end
    legend(['time = ' num2str(ts(l),'%4.2f')])
    pause(0.1)
 end
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:spmdHomoDiff31tFin}final
field~\(u(x,y,z,0.4)\) of the patch scheme applied to a
heterogeneous diffusion~\pde.  }
\includegraphics[scale=0.8]{spmdHomoDiff31tFin}
\end{figure}
Finish the animation loop, and optionally output the final
plot, \cref{fig:spmdHomoDiff31tFin}.
\begin{matlab}
%}
end%for over time
ifOurCf2eps([mfilename 'tFin'])
%{
\end{matlab}




\subsection{\texttt{microBurst} function for Projective Integration}
\label{secmBfPI}
Projective Integration stability seems to need bursts longer
than~\(0.2\). Here take ten meso-steps, each with default
ten micro-steps so the micro-scale step is~\(0.002\). With
macro-step~\(3\), these parameters usually give stable
projective integration (but not always).
\begin{matlab}
%}
function [tbs,xbs] = aBurst(tb0,xb0,patches) 
    normx=max(abs(xb0(:)));
    disp(['aBurst t = ' num2str(tb0) '  |x| = ' num2str(normx)])
    assert(normx<10,'solution exploding')
    tbs = tb0+(0:0.02:0.2);
    spmd
      xb0 = codistributed(xb0,patches.codist);
      xbs = RK2mesoPatch(tbs,xb0,[],patches);
    end%spmd
    xbs=reshape(xbs{1},length(tbs),[]);
end%function
%{
\end{matlab}



%\input{../Patch/heteroDiff3.m}
Fin.
%}
