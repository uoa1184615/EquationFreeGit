% Creates a data struct of the design of 2D patches for
% later use by the patch functions such as smoothPatch2() 
% AJR, Nov 2018 -- Feb 2019
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches2()}: configures spatial patches in 2D}
\label{sec:configPatches2}
\localtableofcontents

\subsection{Introduction}

Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative function~\verb|patchSmooth2()|.
\cref{sec:configPatches2eg} lists an example of its use.


\begin{matlab}
%}
function configPatches2(fun,Xlim,BCs,nPatch,ordCC,ratio,nSubP,nEdge)
global patches
%{
\end{matlab}

\paragraph{Input}
If invoked with no input arguments, then executes an example
of simulating a nonlinear diffusion \pde\ relevant to the
lubrication flow of a thin layer of fluid---see
\cref{sec:configPatches2eg} for the example code.
\begin{itemize}

\item \verb|fun| is the name of the user function,
\verb|fun(t,u,x,y)|, that computes time derivatives (or
time-steps) of quantities on the patches.

\item \verb|Xlim| array/vector giving the macro-space domain
of the computation: patches are distributed equi-spaced over
the interior of the rectangle \([\verb|Xlim(1)|, 
\verb|Xlim(2)|] \times [\verb|Xlim(3)|, \verb|Xlim(4)|]\): 
if \verb|Xlim| is of length two, then use the same interval 
in both directions.

\item \verb|BCs| somehow will define the macroscale boundary
conditions.  Currently, \verb|BCs| is ignored and the system
is assumed macro-periodic in the domain.

\item \verb|nPatch| determines the number of equi-spaced
spaced patches: if scalar, then use the same number of
patches in both directions, otherwise \verb|nPatch(1:2)|
give the number in each direction.

\item \verb|ordCC| is the `order' of interpolation across
empty space of the macroscale mid-patch values to the edge
of the patches for inter-patch coupling: currently must be
in~\(\{0\}\).

\item \verb|ratio| (real) is the ratio of the half-width of
a patch to the spacing of the patch mid-points: so
\(\verb|ratio|=\tfrac12\) means the patches abut; and
\(\verb|ratio|=1\) would be overlapping patches as in
holistic discretisation: if scalar, then use the same ratio
in both directions, otherwise \verb|ratio(1:2)| give the
ratio in each direction.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch: if scalar, then use the same
number in both directions, otherwise \verb|nSubP(1:2)| gives
the number in each direction. Must be odd so that there is a
central lattice point.

\item \verb|nEdge| is, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch.  May be omitted. The default is one (suitable for
microscale lattices with only nearest neighbours.
interactions).

\end{itemize}



\paragraph{Output} The \emph{global} struct \verb|patches|
is created and set with the following components.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(u,t,x,y)| that computes the time derivatives (or
steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.alt| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling---not yet implemented.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the
\(\verb|ordCC|\)-vector of weights for the inter-patch
interpolation onto the right and left edges (respectively)
with patch:macroscale ratio as specified.

\item \verb|.x| is \(\verb|nSubP(1)|\times \verb|nPatch(1)|\)
array of the regular spatial locations~\(x_{ij}\) of the
microscale grid points in every patch.  

\item \verb|.y| is \(\verb|nSubP(2)|\times \verb|nPatch(2)|\)
array of the regular spatial locations~\(y_{ij}\) of the
microscale grid points in every patch.  

\item \verb|.nEdge| is, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch.

\end{itemize}






\subsection{If no arguments, then execute an example}
\label{sec:configPatches2eg}
\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
The code here shows one way to get started: a user's script
may have the following three steps (arrows indicate function
recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches2 
\item ode15s integrator \into patchSmooth2 \into user's nonDiffPDE
\item process results
\end{enumerate}

Establish global patch data struct to interface with a
function coding a nonlinear `diffusion' \pde: to be solved
on \(6\times4\)-periodic domain, with \(9\times7\) patches,
spectral interpolation~(\(0\)) couples the patches, each
patch of half-size ratio~\(0.25\), and with \(5\times5\)
points within each patch. \cite{Roberts2011a} established
that this scheme is consistent with the \pde\ (as the patch 
spacing decreases).
\begin{matlab}
%}
nSubP = 5;
configPatches2(@nonDiffPDE,[-3 3 -2 2], nan, [9 7], 0, 0.25, nSubP);
%{
\end{matlab}
Set a Gaussian initial condition using auto-replication of
the spatial grid.
\begin{matlab}
%}
x = reshape(patches.x,nSubP,1,[],1); 
y = reshape(patches.y,1,nSubP,1,[]);
u0 = exp(-x.^2-y.^2);
u0 = u0.*(0.9+0.1*rand(size(u0)));
%{
\end{matlab}
Initiate a plot of the simulation using only the microscale
values interior to the patches: set \(x\)~and \(y\)-edges to
\verb|nan| to leave the gaps. 
\begin{matlab}
%}
figure(1), clf
x = patches.x; y = patches.y;
x([1 end],:) = nan; y([1 end],:) = nan;
%{
\end{matlab}
Start by showing the initial conditions of
\cref{fig:configPatches2ic} while the simulation computes.
\begin{matlab}
%}
u = reshape(permute(u0,[1 3 2 4]), [numel(x) numel(y)]);
hsurf = surf(x(:),y(:),u');
axis([-3 3 -3 3 -0.001 1]), view(60,40)
legend('time = 0','Location','north')
xlabel('space x'), ylabel('space y'), zlabel('u(x,y)')
drawnow
%{
\end{matlab}
Save the initial condition to file for \cref{fig:configPatches2ic}.
\begin{matlab}
%}
set(gcf,'PaperPosition',[0 0 14 10])
print('-depsc2','configPatches2ic')
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:configPatches2ic}initial
field~\(u(x,y,t)\) at time \(t=0\) of the patch scheme
applied to a nonlinear `diffusion'~\pde:
\cref{fig:configPatches2t3} plots the computed field at time
\(t=3\).}
\includegraphics[scale=0.9]{configPatches2ic}
\end{figure}
Integrate in time using standard functions.
\begin{matlab}
%}
disp('Wait while we simulate h_t=(h^3)_xx+(h^3)_yy')
[ts,ucts] = ode15s(@patchSmooth2,[0 3],u0(:));
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:configPatches2t3}.
\begin{matlab}
%}
for i = 1:length(ts)
  u = patchEdgeInt2(ucts(i,:));
  u = reshape(permute(u,[1 3 2 4]), [numel(x) numel(y)]);
  hsurf.ZData = u';
  legend(['time = ' num2str(ts(i),2)])
  pause(0.1)
end
print('-depsc2','configPatches2t3')
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:configPatches2t3}field~\(u(x,y,t)\) at
time \(t=3\) of the patch scheme applied to a nonlinear
`diffusion'~\pde\ with initial condition in
\cref{fig:configPatches2ic}.}
\includegraphics[scale=0.9]{configPatches2t3}
\end{figure}

Upon finishing execution of the example, exit this function.
\begin{matlab}
%}
return
end%if no arguments
%{
\end{matlab}

\paragraph{Example of nonlinear diffusion PDE inside patches}
As a microscale discretisation of \(u_t=\delsq(u^3)\), code
\(\dot u_{ijkl} =\frac1{\delta x^2} (u_{i+1,j,k,l}^3
-2u_{i,j,k,l}^3 +u_{i-1,j,k,l}^3) + \frac1{\delta y^2}
(u_{i,j+1,k,l}^3 -2u_{i,j,k,l}^3 +u_{i,j-1,k,l}^3)\).
\begin{matlab}
%}
function ut = nonDiffPDE(t,u,x,y)
  dx = diff(x(1:2));  dy = diff(y(1:2));  % microscale spacing
  i = 2:size(u,1)-1;  j = 2:size(u,2)-1;  % interior points in patches
  ut = nan(size(u));  % preallocate storage
  ut(i,j,:,:) = diff(u(:,j,:,:).^3,2,1)/dx^2 ...
               +diff(u(i,:,:,:).^3,2,2)/dy^2;
end
%{
\end{matlab}




\begin{devMan}



\subsection{The code to make patches}

Initially duplicate parameters as needed.
\begin{matlab}
%}
if numel(Xlim)==2, Xlim = repmat(Xlim,1,2); end
if numel(nPatch)==1, nPatch = repmat(nPatch,1,2); end
if numel(ratio)==1, ratio = repmat(ratio,1,2); end
if numel(nSubP)==1, nSubP = repmat(nSubP,1,2); end
%{
\end{matlab}

Set one edge-value to compute by interpolation if not 
specified by the user. Store in the struct.
\begin{matlab}
%}
if nargin<8, nEdge = 1; end
if nEdge>1, error('multi-edge-value interp not yet implemented'), end
if 2*nEdge+1>nSubP, error('too many edge values requested'), end
patches.nEdge = nEdge;
%{
\end{matlab}


First, store the pointer to the time derivative function in
the struct.
\begin{matlab}
%}
patches.fun = fun;
%{
\end{matlab}

Second, store the order of interpolation that is to provide
the values for the inter-patch coupling conditions. Spectral
coupling is \verb|ordCC| of~\(0\) or~\(-1\).
\begin{matlab}
%}
if ~ismember(ordCC,[0])
    error('ordCC out of allowed range [0]')
end
%{
\end{matlab}
For odd~\verb|ordCC| do interpolation based upon odd
neighbouring patches as is useful for staggered grids.
\begin{matlab}
%}
patches.alt = mod(ordCC,2);
ordCC = ordCC+patches.alt;
patches.ordCC = ordCC;
%{
\end{matlab}
%Check for staggered grid and periodic case.
%\begin{matlab}
%%}
%  if patches.alt & (mod(nPatch,2)==1)
%    error('Require an even number of patches for staggered grid')
%  end
%%{
%\end{matlab}
Might as well precompute the weightings for the
interpolation of field values for coupling. (Could sometime
extend to coupling via derivative values.)
\begin{matlab}
%}
ratio = ratio(:)'; % force to be row vector
if patches.alt  % eqn (7) in \cite{Cao2014a}
  patches.Cwtsr = [1
    ratio/2
    (-1+ratio.^2)/8
    (-1+ratio.^2).*ratio/48
    (9-10*ratio.^2+ratio.^4)/384
    (9-10*ratio.^2+ratio.^4).*ratio/3840
    (-225+259*ratio.^2-35*ratio.^4+ratio.^6)/46080
    (-225+259*ratio.^2-35*ratio.^4+ratio.^6).*ratio/645120 ];
else % 
  patches.Cwtsr = [ratio
    ratio.^2/2
    (-1+ratio.^2).*ratio/6
    (-1+ratio.^2).*ratio.^2/24
    (4-5*ratio.^2+ratio.^4).*ratio/120
    (4-5*ratio.^2+ratio.^4).*ratio.^2/720
    (-36+49*ratio.^2-14*ratio.^4+ratio.^6).*ratio/5040
    (-36+49*ratio.^2-14*ratio.^4+ratio.^6).*ratio.^2/40320 ];
end
patches.Cwtsr = patches.Cwtsr(1:ordCC,:);
% should avoid this next implicit auto-replication
patches.Cwtsl = (-1).^((1:ordCC)'-patches.alt).*patches.Cwtsr;
%{
\end{matlab}


Third, set the centre of the patches in a the macroscale
grid of patches assuming periodic macroscale domain.
\begin{matlab}
%}
X = linspace(Xlim(1),Xlim(2),nPatch(1)+1);
X = X(1:nPatch(1))+diff(X)/2;
DX = X(2)-X(1);
Y = linspace(Xlim(3),Xlim(4),nPatch(2)+1);
Y = Y(1:nPatch(2))+diff(Y)/2;
DY = Y(2)-Y(1);
%{
\end{matlab}
Construct the microscale in each patch, assuming Dirichlet
patch edges, and a half-patch length of~\(\verb|ratio(1)|
\cdot \verb|DX|\) and~\(\verb|ratio(2)| \cdot \verb|DY|\).
\begin{matlab}
%}
nSubP = nSubP(:)'; % force to be row vector
if mod(nSubP,2)==[0 0], error('configPatches2: nSubP must be odd'), end
i0 = (nSubP(1)+1)/2;
dx = ratio(1)*DX/(i0-1);
patches.x = bsxfun(@plus,dx*(-i0+1:i0-1)',X); % micro-grid
i0 = (nSubP(2)+1)/2;
dy = ratio(2)*DY/(i0-1);
patches.y = bsxfun(@plus,dy*(-i0+1:i0-1)',Y); % micro-grid
end% function
%{
\end{matlab}
Fin.
\end{devMan}
%}
