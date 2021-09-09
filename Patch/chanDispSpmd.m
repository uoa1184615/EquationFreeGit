% chanDispSpmd simulates 2D shear dispersion in a long thin
% channel with 1D patches as a Proof of Principle example of
% parallel computing with spmd.    AJR, Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{chanDispSpmd}: simulation of a 1D shear
dispersion via simulation on small patches across a channel}
\label{sec:chanDispSpmd}
\localtableofcontents

Simulate 1D shear dispersion along long thin channel,
dispersion that is emergent from micro-scale dynamics in 2D
space.  Use 1D patches as a Proof of Principle example of
parallel computing with \verb|spmd|. In this shear
dispersion, although the micro-scale diffusivities are
one-ish, the shear causes an effective longitudinal
`diffusivity' of the order of~$\Pe^2$\,---which is typically
much larger than the micro-scale diffusivity
\cite[e.g.]{Taylor53}.

The spatial domain is the channel (large) $L$-periodic
in~$x$ and $|y|<1$\,. Seek to predict a concentration
field~$c(x,y,t)$ satisfying the linear
advection-diffusion~\pde
\begin{equation}
\D tc = -\Pe u(y)\D xc+\D x{}\Big[\kappa_x(y)\D xc\Big]
+\D y{}\Big[\kappa_y(y)\D yc\Big].
\label{eq:pdeChanDisp}
\end{equation}
where \Pe\ denotes a Peclet number, parabolic advection
velocity $u(y)=\tfrac32(1-y^2)$ with noise, and parabolic
diffusivity  $\kappa_x(y)=\kappa_y(y)=(1-y^2)$ with noise.
The noise is to be multiplicative and log-normal to ensure
advection and diffusion are all positive, and to be periodic
in~$x$.

For a microscale computation we discretise in space with
$x$-spacing~$\delta x$, and $n_y$~points over $|y|<1$ with
spacing $\delta y:=2/n_y$ at $y_j:=-1+(j-\tfrac12)\delta
y$\,, $j=1:n_y$\,.  Our microscale discretisation of
\pde~\eqref{eq:pdeChanDisp} is then
\begin{align}&
\D t{c_{ij}}=-\Pe u(y_j)\frac{c_{i+1,j}-c_{i-1,j}}{2\delta x}
+\frac{d_{i,j+1/2}-d_{i,j-1/2}}{\delta y}
+\frac{D_{i+1/2,j}-D_{i-1/2,j}}{\delta x}\,,
\nonumber\\&\quad
d_{ij}:=\kappa_y(y_j)\frac{c_{i,j+1/2}-c_{i,j-1/2}}{\delta y}\,,\quad
D_{ij}:=\kappa_x(y_j)\frac{c_{i+1/2,j}-c_{i-1/2,j}}{\delta x}\,.
\label{eq:ddeChanDisp}
\end{align}
These are coded in \cref{sec:chanDispMicro} for the computation.



Choose one of four cases:
\begin{itemize}
\item \verb|theCase=1| is corresponding code without
parallelisation (in this toy problem it is much the quickest
because there is no expensive interprocessor communication);
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
The micro-scale \pde\ is evaluated at positions~$y_j$ across
the channel, $|y|<1$\,.  The even indexed points are the
collocation points for the \pde, whereas the odd indexed
points are the half-grid points for specification of
$y$-diffusivities.
\begin{matlab}
%}
ny = 7
y = linspace(-1,1,2*ny+1);
yj = y(2:2:end);
%{
\end{matlab}
Set micro-scale advection (array~1) and diffusivity
(array~2) with (roughly) parabolic shape
\cite[e.g.]{Watt94b, MacKenzie03}. Here modify the parabola
by a heterogeneous log-normal factor with specified period
along the channel: modify the strength of the heterogeneity
by the coefficient of~\verb|randn| from zero to perhaps one:
coefficient~$0.3$ appears a good moderate value. Remember
that \verb|configPatches1| reshapes \verb|cHetr| to~2D.
\begin{matlab}
%}
mPeriod = 4
cHetr = shiftdim([3/2 1],-1).*(1-y.^2) ...
        .*exp(0.3*randn([mPeriod 2*ny+1 2]));
%{
\end{matlab}


Configure the patch scheme with some arbitrary choices of
domain, patches, size ratios.  Choose some random order of
interpolation to see the alternatives. Set \verb|patches|
information to be global so the info can be used for
Cases~1--2 without being explicitly passed as arguments. 
Choose the parallel option if not Case~1, which invokes
\verb|spmd|-block internally, so that field variables become
\emph{distributed} across cpus.
\begin{matlab}
%}
if theCase<=2, global patches, end
nPatch=15
nSubP=2+mPeriod
ratio=0.2+0.2*(theCase<4)
Len=nPatch/ratio 
ordCC=2*randi([0 3])
disp('**** Setting configPatches1')
patches = configPatches1(@chanDispMicro, [0 Len], nan ...
    , nPatch, ordCC, ratio, nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',cHetr ,'parallel',(theCase>1) );
%{
\end{matlab}
When using parallel then additional parameters to
\verb|patches| should be set within a \verb|spmd| block
(because \verb|patches| is a co-distributed structure).
\begin{matlab}
%}
Peclet = 10
if theCase==1, patches.Pe = Peclet;
else     spmd, patches.Pe = Peclet; end
end
%{
\end{matlab}



\subsection{Simulate heterogeneous advection-diffusion}
Set initial conditions of a simulation as shown in
\cref{fig:chanDispSpmdt0}.
\begin{matlab}
%}
disp('**** Set initial condition and test dc0dt =')
if theCase==1
%{
\end{matlab}
Without parallel processing, invoke the usual operations.
\begin{matlab}
%}
    c0 = 10*exp(-(ratio*patches.x-2.5).^2/2) +0*yj;
    c0 = c0.*(1+0.2*rand(size(c0)));
    dc0dt = patchSys1(0,c0);
%{
\end{matlab}
With parallel, we must use an \verb|spmd|-block for
computations: there is no difference in cases~2--4 here.
Also, we must sometimes use \verb|patches.codist| to
explicitly code how to distribute new arrays over the cpus.
Now \verb|patchSys1| does not invoke \verb|spmd| so
higher level code must, as here. Even if \verb|patches| is
global, inside \verb|spmd|-block we \emph{must} pass it
explicitly as a parameter to \verb|patchSys1|.
\begin{matlab}
%}
else, spmd
    c0 = 10*exp(-(ratio*patches.x-2.5).^2/2) +0*yj;
    c0 = c0.*(1+0.2*rand(size(c0),patches.codist));
    dc0dt = patchSys1(0,c0,patches)
    end%spmd
end%if theCase
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:chanDispSpmdt0}initial
field~$u(x,y,0)$ of the patch scheme applied to a
heterogeneous advection-diffusion~\pde. 
\cref{fig:chanDispSpmdtFin} plots the roughly smooth field
values at time $t=4$.  In this example the patches are
relatively large, ratio~$0.4$, for visibility.}
\includegraphics[scale=0.8]{chanDispSpmdt0}
\end{figure}

Integrate in time, either via the automatic \verb|ode23| or
via \verb|RK2mesoPatch| which reduces communication between
patches. By default, \verb|RK2mesoPatch| does ten
micro-steps for each specified meso-step in~\verb|ts|.  For
stability: with noise up to~$0.3$, need micro-steps less
than~$0.005$; with noise~$1$, need micro-steps less
than~$0.0015$.
\begin{matlab}
%}
warning('Integrating system in time, wait patiently')
ts=4*linspace(0,1); 
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
\verb|ode23| as here \verb|patchSys1| accesses patch
information via global \verb|patches|. For post-processing,
reshape each and every row of the computed solution to the
correct array size---namely that of the initial condition. 
\begin{matlab}
%}
case 1
%    [cs,uerrs] = RK2mesoPatch(ts,c0);
    [ts,cs] = ode23(@patchSys1,ts,c0(:));
    cs=reshape(cs,[length(ts) size(c0)]);
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
    cs = RK2mesoPatch(ts,c0);
    cs = cs{1};  
%{
\end{matlab}

\item In this third case, a user could merge this explicit
\verb|spmd|-block with the previous one that sets the
initial conditions.
\begin{matlab}
%}
case 3,spmd
    cs = RK2mesoPatch(ts,c0,[],patches);
    end%spmd
    cs = cs{1};
%{
\end{matlab}

\item In this fourth case, use Projective Integration (PI)
over long times (\verb|PIRK4| also works). Currently the PI
is done serially, with parallel \verb|spmd|-blocks only
invoked inside function \verb|aBurst()| (\cref{secmBfPI}) to
compute each burst of the micro-scale simulation. For a
Peclet number of ten, the macro-scale time-step needs to be
less than about~$0.5$ (which here is very little
projection)---presumably the mean advection in a macro-step
needs to be less than about the patch spacing. The function
\verb|microBurst()| here interfaces to \verb|aBurst()|
(\cref{secCHS1mBfPI}) in order to provide shaped initial
states, and to provide the patch information.
\begin{matlab}
%}
case 4
    microBurst = @(tb0,xb0,bT) ...
        aBurst(tb0 ,reshape(xb0,size(c0)) ,patches);
    ts = 0:0.7:5
    cs = PIRK2(microBurst,ts,gather(c0(:)));
    cs = reshape(cs,[length(ts) size(c0)]);
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
figure(1), clf, colormap(0.8*hsv)
%{
\end{matlab}
First get the $x$-coordinates and omit the patch-edge
values from the plot (because they are not here
interpolated).
\begin{matlab}
%}
if theCase==1, x = patches.x;  
else, spmd
    x = gather( patches.x ); 
    end%spmd
    x = x{1};
end
x([1 end],:,:,:) = nan; 
%{
\end{matlab}
For every time step draw the concentration values as a set
of surfaces on 2D patches, with a short pause to display
animation.
\begin{matlab}
%}
nTimes = length(ts)
for l = 1:nTimes
%{
\end{matlab}
At each time, squeeze sub-patch data into a 3D array,
permute to get all the $x$-variation in the first two
dimensions, and reshape into $x$-variation for each and
every~$(y)$. 
\begin{matlab}
%}
  c = reshape( permute( squeeze( ...
      cs(l,:,:,:,:) ) ,[1 3 2]) ,numel(x),ny);
%{
\end{matlab}
Draw surface of each patch, to show both micro-scale and
macro-scale variation in space. 
\begin{matlab}
%}
  if l==1
    hp = surf(x(:),yj,c');
    axis([0 Len -1 1 0 max(c(:))])
    axis equal
    xlabel('space x'), ylabel('y'); zlabel('c(x,y,t)')
    ifOurCf2eps([mfilename 't0'])
    legend(['time = ' num2str(ts(l),'%4.2f')] ...
        ,'Location','north')
    disp('**** pausing, press blank to animate')
    pause
 else
    hp.ZData = c'; 
    legend(['time = ' num2str(ts(l),'%4.2f')])
    pause(0.1)
 end
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:chanDispSpmdtFin}final
field~$c(x,y,4)$ of the patch scheme applied to a
heterogeneous advection-diffusion
\pde~\eqref{eq:pdeChanDisp} with heterogeneous factor
log-normal, here distributed $\exp[\mathcal N(0,1)]$.  }
\includegraphics[scale=0.8]{chanDispSpmdtFin}
\end{figure}

Finish the animation loop, and optionally save the final
plot to file, \cref{fig:chanDispSpmdtFin}.
\begin{matlab}
%}
end%for over time
ifOurCf2eps([mfilename 'tFin'])
%{
\end{matlab}

\paragraph{Macro-scale view}
Plot a macro-scale mesh of the predictions: at each of a
selection of times, for every patch, plot the patch-mean
value at the mean-$x$. 
\begin{matlab}
%}
figure(2), clf, colormap(0.8*hsv)
X = squeeze(mean(x(2:end-1,:,:,:)));
C = squeeze(mean(mean(cs(:,2:end-1,:,:,:),2),3));
j = 1:ceil(nTimes/30):nTimes;
mesh(X,ts(j),C(j,:));
xlabel('space x'),ylabel('time t'),zlabel('C(X,t)')
zlim([-0.1 11])
ifOurCf2eps([mfilename 'Macro'])
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:chanDispSpmdMacro}macro-scale
view of heterogeneous advection-diffusion~\pde\ along a
(periodic) channel obtained via the patch scheme.  }
\includegraphics[scale=0.8]{chanDispSpmdMacro}
\end{figure}




\subsection{\texttt{microBurst} function for Projective Integration}
\label{secCHS1mBfPI}
Projective Integration stability appears to require bursts
longer than~$0.2$.  Each burst is done in parallel
processing.  Here use \verb|RK2mesoPatch| to take take
meso-steps, each with default ten micro-steps so the
micro-scale step is~$0.0033$. With macro-step~$0.5$,
these parameters usually give stable projective integration.
\begin{matlab}
%}
function [tbs,xbs] = aBurst(tb0,xb0,patches) 
    normx=max(abs(xb0(:)));
    disp(['* aBurst t=' num2str(tb0) '  |x|=' num2str(normx)])
    assert(normx<20,'solution exploding')
    tbs = tb0+(0:0.033:0.2);
    spmd
      xb0 = codistributed(xb0,patches.codist);
      xbs = RK2mesoPatch(tbs,xb0,[],patches);
    end%spmd
    xbs=reshape(xbs{1},length(tbs),[]);
end%function
%{
\end{matlab}
Fin.



\input{../Patch/chanDispMicro.m}
%}
