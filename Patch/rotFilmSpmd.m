% rotFilmSpmd simulates 2D fluid film flow on a rotating
% substrate with 2D patches as a Proof of Principle example
% of parallel computing with spmd.   AJR, Dec 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{rotFilmSpmd}: simulation of a 2D shallow
water flow on a rotating heterogeneous substrate}
\label{sec:rotFilmSpmd}
\localtableofcontents


As an example application, consider the flow of a shallow
layer of fluid on a solid flat rotating substrate, such as
in spin coating \citep[\S II.K, e.g.]{Wilson00, Oron97} or
large-scale shallow water waves \cite[e.g.]{Dellar2005,
Hereman2009}. Let $\xv=(x,y)$ parametrise location on the
rotating substrate, and let the fluid layer have
thickness~$h(\xv , t)$ and move with depth-averaged
horizontal velocity $\vv (\xv , t)=(u,v)$. We take as given
(with its simplified physics) that the (non-dimensional)
governing set of \pde{}s is the nonlinear system
\cite[eq.~(1), e.g.]{Bunder2018a}
\begin{subequations}\label{eqs:spinddt}%
\begin{align}
\D th&=-{\nabla}\cdot (h\vv ), \label{eq:spindhdt}\\
\D t\vv&=\begin{bmatrix} -b & f\\-f &-b\end{bmatrix}\vv 
-(\vv \cdot\nabla)\vv -g{\nabla}h+\divv(\nu\grad\vec{v})\,,
\label{eq:spindvdt}
\end{align}
\end{subequations}
where $b(\xv)$~represents heterogeneous `bed' drag, $f$~is
the Coriolis coefficient,  $g$~is the acceleration due to
gravity, $\nu(\xv)$~is a heterogeneous `kinematic
viscosity', and we neglect surface tension. 

The aim is to simulate the macroscale dynamics which (for
constant~$b$) is approximately that of the nonlinear
diffusion $\D th\approx \frac{gb}{b^2+f^2}\divv(h\grad h)$
\cite[eq.~(2)]{Bunder2018a}. But there is no known algebraic
closure for the macroscale in the case of
heterogeneous~$b(\xv)$ and~$\nu(\xv)$, nonetheless the patch
scheme automatically predicts a sensible macroscale for such
heterogeneous dynamics (\cref{fig:rotFilmSpmdtFin}).

For the microscale computation, \cref{sec:rotFilmMicro}
discretises the \pde{}s~\eqref{eqs:spinddt} in space with
$x,y$-spacing~$\delta x,\delta y$.



Choose one of four cases:
\begin{itemize}
\item \verb|theCase=1| is corresponding code without
parallelisation (in this toy problem it is much the quickest
because there is no expensive communication);
\item \verb|theCase=2| illustrates that \verb|RK2mesoPatch|
invokes \verb|spmd| computation if parallel has been
configured.
\item \verb|theCase=3| shows how users explicitly invoke
\verb|spmd|-blocks around the time integration.
\item \verb|theCase=4| invokes projective integration for
long-time simulation via short bursts of the
micro-computation, bursts done within \verb|spmd|-blocks for
parallel computing.
\end{itemize}
First, clear all to remove any existing globals, old
composites, etc---although a parallel pool persists. Then
choose the case.
\begin{matlab}
%}
clear all
theCase = 1
%{
\end{matlab}
Set micro-scale bed drag (array~1) and diffusivity
(arrays~2--3) to be a heterogeneous log-normal factor with
specified period: modify the strength of the heterogeneity
by the coefficient of~\verb|randn| from zero to perhaps one:
coefficient~$0.3$ appears a good moderate value. 
\begin{matlab}
%}
mPeriod = 5
bnu = shiftdim([1 0.5 0.5],-1) ...
      .*exp(0.3*randn([mPeriod mPeriod 3]));
%{
\end{matlab}


Configure the patch scheme with these choices of domain,
patches, size ratios---here each patch is square in space.
In Cases~1--2, set \verb|patches| information to be global
so the info can be used without being explicitly passed as
arguments.
\begin{matlab}
%}
if theCase<=2, global patches, end
%{
\end{matlab}
In Case~4, double the size of the domain and use more
separated patches accordingly, to maintain the spatial
microscale grid spacing to be~$0.055$.  Here use fourth
order edge-based coupling between patches.  Choose the
parallel option if not Case~1, which invokes
\verb|spmd|-block internally, so that field variables become
\emph{distributed} across cpus.
\begin{matlab}
%}
nSubP = 2+mPeriod
nPatch = 9
ratio = 0.2+0.2*(theCase<4)
Len = 2*pi*(1+(theCase==4))
disp('**** Setting configPatches2')
patches = configPatches2(@rotFilmMicro, [0 Len], nan ...
    , nPatch, 4, ratio, nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',bnu ,'parallel',(theCase>1) );
%{
\end{matlab}
When using parallel, any additional parameters to
\verb|patches|, such as physical parameters for the
microcode, must be set within a \verb|spmd| block (because
\verb|patches| is a co-distributed structure). Here set
frequency of substrate rotation, and strength of gravity.
\begin{matlab}
%}
f = 5,  g = 1
if theCase==1, patches.f = f; patches.g = g;
else     spmd, patches.f = f; patches.g = g; end
end
%{
\end{matlab}



\subsection{Simulate heterogeneous advection-diffusion}
Set initial conditions of a simulation as shown in
\cref{fig:rotFilmSpmdt0}.  Here the initial condition is a
(periodic) quasi-Gaussian in~$h$ and zero velocity~\vv, with
additive random perturbations.
\begin{matlab}
%}
disp('**** Set initial condition and test dhuv0dt =')
if theCase==1
%{
\end{matlab}
When not parallel processing, invoke the usual operations.
Here add a random noise to the velocity field, but
keep~$h(x,y,0)$ smooth as shown by \cref{fig:rotFilmSpmdt0}.
The \verb|shiftdim(...,-1)| moves the given row-vector of
coefficients into the third dimension to become coefficients
of the fields~$(h,u,v)$, respectively.
\begin{matlab}
%}
    huv0 = shiftdim([0.5 0 0],-1) ...
       .*exp(-cos(patches.x)/2-cos(patches.y));
    huv0 = huv0+0.1*shiftdim([0 1 1],-1).*rand(size(huv0));
    dhuv0dt = patchSys2(0,huv0);
%{
\end{matlab}
With parallel, we must use an \verb|spmd|-block for
computations: there is no difference in Cases~2--4 here.
Also, we must sometimes explicitly tell functions how to
distribute some initial condition arrays over the cpus. Now
\verb|patchSys2| does not invoke \verb|spmd| so higher
level code must, as here. Even if \verb|patches| is global,
inside an \verb|spmd|-block we \emph{must} pass
\verb|patches| explicitly as a parameter to
\verb|patchSys2|.
\begin{matlab}
%}
else, spmd
    huv0 = shiftdim([0.5 0 0],-1) ...
       .*exp(-cos(patches.x)/2-cos(patches.y));
    huv0 = huv0+0.1*rand(size(huv0),patches.codist);
    dhuv0dt = patchSys2(0,huv0,patches)
    end%spmd
end%if theCase
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:rotFilmSpmdt0}initial
field~$h(x,y,0)$ of the patch scheme applied to the
heterogeneous, shallow water, rotating substrate,
\pde~\eqref{eqs:spinddt}. The micro-scale sub-patch colour
displays the initial $y$-direction velocity
field~$v(x,y,0)$. \cref{fig:rotFilmSpmdtFin} plots the
roughly smooth field values at time $t=6$.  In this
example the patches are relatively large, ratio~$0.4$, for
visibility.}
\includegraphics[scale=0.8]{rotFilmSpmdt0}
\end{figure}

Integrate in time, either via the automatic \verb|ode23| or
via \verb|RK2mesoPatch| which reduces communication between
patches. By default, \verb|RK2mesoPatch| does ten
micro-steps for each specified meso-step in~\verb|ts|.  For
stability: with noise up to~$0.3$, need micro-steps less
than~$0.0003$; with noise~$1$, need micro-steps less
than~$0.0001$.
\begin{matlab}
%}
warning('Integrating system in time, wait a minute')
ts=0:0.003:0.3; 
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
\verb|ode23| as here \verb|patchSys2| accesses patch
information via global \verb|patches|. For post-processing,
reshape each and every row of the computed solution to the
correct array size---namely that of the initial condition. 
\begin{matlab}
%}
case 1
%    tic,[huvs,uerrs] = RK2mesoPatch(ts,huv0);toc
    [ts,huvs] = ode23(@patchSys2,[0 4],huv0(:));
    huvs=reshape(huvs,[length(ts) size(huv0)]);
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
    huvs = RK2mesoPatch(ts,huv0);
    huvs = huvs{1};  
%{
\end{matlab}

\item In this third case, a user could merge this explicit
\verb|spmd|-block with the previous one that sets the
initial conditions.
\begin{matlab}
%}
case 3,spmd
    huvs = RK2mesoPatch(ts,huv0,[],patches);
    end%spmd
    huvs = huvs{1};
%{
\end{matlab}

\item In this fourth case, use Projective Integration (PI).
Currently the PI is done serially, with parallel
\verb|spmd|-blocks only invoked inside function
\verb|aBurst()| (\cref{secRF2BfPI}) to compute each burst of
the micro-scale simulation. The macro-scale time-step needs
to be less than about~$0.1$ (which here is not much
projection). The function \verb|microBurst()| interfaces to
\verb|aBurst()| (\cref{secRF2BfPI}) in order to provide
shaped initial states, and to provide the patch information.
\begin{matlab}
%}
case 4
    microBurst = @(tb0,xb0,bT) ...
        aBurst(tb0 ,reshape(xb0,size(huv0)) ,patches);
    ts = 0:0.1:1
    huvs = PIRK2(microBurst,ts,gather(huv0(:)));
    huvs = reshape(huvs,[length(ts) size(huv0)]);
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
Optionally set to save some plots to file.
\begin{matlab}
%}
if 0, global OurCf2eps, OurCf2eps=true, end
%{
\end{matlab}
\paragraph{Animate the computed solution field over time}
\begin{matlab}
%}
figure(1), clf, colormap(0.8*jet)
%{
\end{matlab}
First get the $x$-coordinates and omit the patch-edge values
from the plot (because they are not here interpolated).
\begin{matlab}
%}
if theCase==1, x = patches.x; 
               y = patches.y; 
else, spmd
    x = gather( patches.x ); 
    y = gather( patches.y ); 
    end%spmd
    x = x{1}; y = y{1};
end
x([1 end],:,:,:,:,:) = nan; 
y(:,[1 end],:,:,:,:) = nan; 
%{
\end{matlab}
Draw the field values as a patchy surface evolving over
100--200 time steps.
\begin{matlab}
%}
nTimes = length(ts)
for l = 1:ceil(nTimes/200):nTimes
%{
\end{matlab}
At each time, squeeze sub-patch data fields into three 4D
arrays, permute to get all the $x/y$-variations in the
first/last two dimensions, and and then reshape to~2D.   
\begin{matlab}
%}
  h = reshape( permute( squeeze( ...
      huvs(l,:,:,1,1,:,:) ) ,[1 3 2 4]) ,numel(x),numel(y));
  u = reshape( permute( squeeze( ...
      huvs(l,:,:,2,1,:,:) ) ,[1 3 2 4]) ,numel(x),numel(y));
  v = reshape( permute( squeeze( ...
      huvs(l,:,:,3,1,:,:) ) ,[1 3 2 4]) ,numel(x),numel(y));
%{
\end{matlab}
Draw surface of each patch, to show both micro-scale and
macro-scale variation in space. Colour the surface according
to the velocity~$v$ in the $y$-direction.
\begin{matlab}
%}
  if l==1
    hp = surf(x(:),y(:),h',v');
    axis([0 Len 0 Len 0 max(h(:))])
    c = colorbar; c.Label.String = 'v(x,y,t)';
    legend(['time = ' num2str(ts(l),'%4.2f')] ...
        ,'Location','north')
    axis equal
    xlabel('space x'), ylabel('space y'), zlabel('h(x,y,t)')
    ifOurCf2eps([mfilename 't0'])
    disp('**** pausing, press blank to begin animation')
    pause
 else
    hp.ZData = h'; hp.CData = v';
    legend(['time = ' num2str(ts(l),'%4.2f')])
    pause(0.1)
 end
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:rotFilmSpmdtFin}final
field~$h(x,y,6)$, coloured by~$v(x,y,6)$, of the patch
scheme applied to the heterogeneous, shallow water, rotating
substrate, \pde~\eqref{eqs:spinddt} with heterogeneous
factors log-normal, here distributed $\exp[\mathcal
N(0,1)]$.  }
\includegraphics[scale=0.8]{rotFilmSpmdtFin}
\end{figure}

Finish the animation loop, and optionally save the final
plot to file, \cref{fig:rotFilmSpmdtFin}.
\begin{matlab}
%}
end%for over time
ifOurCf2eps([mfilename 'tFin'])
%{
\end{matlab}




\subsection{\texttt{microBurst} function for Projective Integration}
\label{secRF2BfPI}
Projective Integration stability appears to require bursts
longer than~$0.01$.  Each burst is done in parallel
processing.  Here use \verb|RK2mesoPatch| to take take
meso-steps, each with default ten micro-steps so the
micro-scale step is~$0.0003$. With macro-step~$0.1$, these
parameters usually give stable projective integration.
\begin{matlab}
%}
function [tbs,xbs] = aBurst(tb0,xb0,patches) 
    normx=max(abs(xb0(:)));
    disp(['* aBurst t=' num2str(tb0) '  |x|=' num2str(normx)])
    assert(normx<20,'solution exploding')
    tbs = tb0+(0:0.003:0.015);
    spmd
      xb0 = codistributed(xb0,patches.codist);
      xbs = RK2mesoPatch(tbs,xb0,[],patches);
    end%spmd
    xbs=reshape(xbs{1},length(tbs),[]);
end%function
%{
\end{matlab}
Fin.



\input{../Patch/rotFilmMicro.m}
%}
