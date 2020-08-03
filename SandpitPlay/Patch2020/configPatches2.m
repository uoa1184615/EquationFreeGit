% Creates a data struct of the design of 2D patches for
% later use by the patch functions such as smoothPatch2() 
% AJR, Nov 2018 -- Jul 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches2()}: configures spatial patches in 2D}
\label{sec:configPatches2}
\localtableofcontents

\subsection{Introduction}

Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative\slash step function
\verb|patchSmooth2()|. \cref{sec:configPatches2eg} lists an
example of its use.
\begin{matlab}
%}
function configPatches2(fun,Xlim,BCs,nPatch,ordCC,ratio,nSubP ...
    ,varargin)
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
\verb|fun(t,u,x,y)|, that computes time-derivatives (or
time-steps) of quantities on the 2D micro-grid within all
the 2D patches.

\item \verb|Xlim| array/vector giving the macro-space domain
of the computation: patches are distributed equi-spaced over
the interior of the rectangle \([\verb|Xlim(1)|,
\verb|Xlim(2)|] \times [\verb|Xlim(3)|, \verb|Xlim(4)|]\):
if \verb|Xlim| is of length two, then the domain is the
square domain of the same interval in both directions.

\item \verb|BCs| eventually and somehow will define the
macroscale boundary conditions.  Currently, \verb|BCs| is
ignored and the system is assumed macro-periodic in the
domain.

\item \verb|nPatch| determines the number of equi-spaced
spaced patches: if scalar, then use the same number of
patches in both directions, otherwise \verb|nPatch(1:2)|
gives the number of patches in each direction.

\item \verb|ordCC| is the `order' of interpolation for
inter-patch coupling across empty space of the macroscale
patch values to the edge-values of the patches:
currently must be~\(0\); where \(0\)~gives spectral
interpolation.
??

\item \verb|ratio| (real) is the ratio of (depending upon
\verb|EdgyInt|) either the half-width or full-width of a
patch to the spacing of the patch mid-points.  So either
\(\verb|ratio|=\tfrac12\) means the patches abut and
\(\verb|ratio|=1\) is overlapping patches as in holistic
discretisation, or \(\verb|ratio|=1\) means the patches
abut.  Small~\verb|ratio| should greatly reduce
computational time.  If scalar, then use the same ratio in
both directions, otherwise \verb|ratio(1:2)| gives the ratio
in each direction.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch: if scalar, then use the same
number in both directions, otherwise \verb|nSubP(1:2)| gives
the number in each direction. If not using \verb|EdgyInt|,
then must be odd so that there is a central micro-grid point
in each patch.

\item \verb|'nEdge'| (not yet implemented), \emph{optional},
default=1, for each patch, the number of edge values set by
interpolation at the edge regions of each patch.  The
default is one (suitable for microscale lattices with only
nearest neighbour interactions).

\item \verb|'EdgyInt'|, true/false, \emph{optional},
default=false.  If true, then interpolate to left\slash
right\slash top\slash bottom edge-values from right\slash
left\slash bottom\slash top next-to-edge values.  So far
only implemented for spectral interpolation,
\(\verb|ordCC|=0\).

\item \verb|'nEnsem'|,  \emph{optional-experimental},
default one, but if more, then an ensemble over this
number of realisations.

\end{itemize}



\paragraph{Output} The \emph{global} struct \verb|patches|
is created and set with the following components.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x,y)| that computes the time derivatives (or
steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.alt| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling---not yet
implemented.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the
\(\verb|ordCC|\)-vector of weights for the inter-patch
interpolation onto the right and left edges (respectively)
with patch:macroscale ratio as specified---not yet
implemented.

\item \verb|.x| (6D) is \(\verb|nSubP(1)| \times1 \times1
\times1 \times \verb|nPatch(1)| \times1\) array of the
regular spatial locations~\(x_{ij}\) of the microscale grid
points in every patch.  

\item \verb|.y| (6D) is \(1 \times \verb|nSubP(2)| \times1
\times1 \times1 \times \verb|nPatch(2)|\) array of the
regular spatial locations~\(y_{ij}\) of the microscale grid
points in every patch.  

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
\item ode15s integrator \into patchSmooth2 \into user's PDE
\item process results
\end{enumerate}

Establish global patch data struct to interface with a
function coding a nonlinear `diffusion' \pde: to be solved
on \(6\times4\)-periodic domain, with \(9\times7\) patches,
spectral interpolation~(\(0\)) couples the patches, each
patch of half-size ratio~\(0.4\) (relatively large for
visualisation), and with \(5\times5\) points forming each
patch. \cite{Roberts2011a} established that this scheme is
consistent with the \pde\ (as the patch spacing decreases).
\begin{matlab}
%}
nSubP = 5;
configPatches2(@nonDiffPDE,[-3 3 -2 2], nan, [9 7] ...
    , 0, 0.4, nSubP ,'EdgyInt',true);
%{
\end{matlab}
Set a  perturbed-Gaussian initial condition using
auto-replication of the spatial grid.
\begin{matlab}
%}
u0 = exp(-patches.x.^2-patches.y.^2);
u0 = u0.*(0.9+0.1*rand(size(u0)));
%{
\end{matlab}
Initiate a plot of the simulation using only the microscale
values interior to the patches: set \(x\)~and \(y\)-edges to
\verb|nan| to leave the gaps between patches. 
\begin{matlab}
%}
figure(1), clf, colormap(hsv)
x = squeeze(patches.x);  y = squeeze(patches.y);
if 1, x([1 end],:) = nan; y([1 end],:) = nan; end
%{
\end{matlab}
Start by showing the initial conditions of
\cref{fig:configPatches2ic} while the simulation computes.
\begin{matlab}
%}
u = reshape(permute(squeeze(u0) ...
    ,[1 3 2 4]), [numel(x) numel(y)]);
hsurf = surf(x(:),y(:),u');
axis([-3 3 -3 3 -0.03 1]), view(60,40)
legend('time = 0.00','Location','north')
xlabel('space x'), ylabel('space y'), zlabel('u(x,y)')
colormap(hsv)
%{
\end{matlab}
Optionally save the initial condition to graphic file for
\cref{fig:configPatches2ic}.
\begin{matlab}
%}
ourcf2eps([mfilename 'ic'])
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
Integrate in time using standard functions. In Matlab
\verb|ode15s| would be natural as the patch scheme is
naturally stiff, but \verb|ode23| is quicker.  Ask for
output at non-uniform times as the diffusion slows.
\begin{matlab}
%}
disp('Wait while simulating nonlinear diffusion h_t=(h^3)_xx+(h^3)_yy')
drawnow
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSmooth2,linspace(0,2).^2,u0(:));
else % octave version is quite slow for me
    lsode_options('absolute tolerance',1e-4);
    lsode_options('relative tolerance',1e-4);
    [ts,us] = odeOcts(@patchSmooth2,[0 1],u0(:));
end
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:configPatches2t3}.  Use \verb|patchEdgeInt2| 
to interpolate patch-edge values (even if not drawn).
\begin{matlab}
%}
for i = 1:length(ts)
  u = patchEdgeInt2(us(i,:));
  u = reshape(permute(squeeze(u) ...
      ,[1 3 2 4]), [numel(x) numel(y)]);
  set(hsurf,'ZData', u');
  legend(['time = ' num2str(ts(i),'%4.2f')])
  pause(0.1)
end
ourcf2eps([mfilename 't3'])
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


\input{../Patch/nonDiffPDE.m}



\begin{devMan}

\subsection{Parse input arguments and defaults}
\begin{matlab}
%}
p = inputParser;
addRequired(p,'fun'); % is there test for a function name??
addRequired(p,'Xlim',@isnumeric);
addRequired(p,'BCs'); % nothing yet decided
addRequired(p,'nPatch',@isnumeric);
addRequired(p,'ordCC',@isnumeric);
addRequired(p,'ratio',@isnumeric);
addRequired(p,'nSubP',@isnumeric);
addParameter(p,'nEdge',1,@isnumeric);
addParameter(p,'EdgyInt',false,@islogical);
addParameter(p,'nEnsem',1,@isnumeric);
%addParameter(p,'nCore',1,@isnumeric); % not yet implemented
parse(p,fun,Xlim,BCs,nPatch,ordCC,ratio,nSubP,varargin{:});
%{
\end{matlab}
Set the optional parameters. 
\begin{matlab}
%}
patches.nEdge = p.Results.nEdge;
patches.EdgyInt = p.Results.EdgyInt;
patches.nEnsem = p.Results.nEnsem;
%patches.nCore = p.Results.nCore;
%{
\end{matlab}

Initially duplicate parameters for both space dimensions as
needed.
\begin{matlab}
%}
if numel(Xlim)==2, Xlim = repmat(Xlim,1,2); end
if numel(nPatch)==1, nPatch = repmat(nPatch,1,2); end
if numel(ratio)==1, ratio = repmat(ratio,1,2); end
if numel(nSubP)==1, nSubP = repmat(nSubP,1,2); end
%{
\end{matlab}

Check parameters.
\begin{matlab}
%}
assert(patches.nEdge==1 ...
      ,'multi-edge-value interp not yet implemented')
assert(all(2*patches.nEdge<nSubP) ...
      ,'too many edge values requested')
%if patches.nCore>1
%    warning('nCore>1 not yet tested in this version')
%    end
%{
\end{matlab}





\subsection{The code to make patches}

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
assert((ordCC>=-1) & (floor(ordCC)==ordCC), ...
    'ordCC out of allowed range integer>=-1')
%{
\end{matlab}
For odd~\verb|ordCC| do interpolation based upon odd
neighbouring patches as is useful for staggered grids.
\begin{matlab}
%}
patches.alt = mod(ordCC,2);
assert(patches.alt==0,'staggered not yet implemented??')
ordCC = ordCC+patches.alt;
patches.ordCC = ordCC;
%{
\end{matlab}
Check for staggered grid and periodic case.
\begin{matlab}
%}
if patches.alt, assert(mod(nPatch,2)==0, ...
    'Require an even number of patches for staggered grid')
end
%{
\end{matlab}
Might as well precompute the weightings for the
interpolation of field values for coupling---not yet used
here?? (Could sometime extend to coupling via derivative
values.)
\begin{matlab}
%}
ratio = reshape(ratio,1,2); % force to be row vector
if ordCC>0
  if patches.alt  % eqn (7) in \cite{Cao2014a}
    patches.Cwtsr(1:2:ordCC,1:2)=[1 ...
      cumprod((ratio.^2-(1:2:(ordCC-2))'.^2)/4)./ ...
      factorial(2*(1:(ordCC/2-1))')];    
    patches.Cwtsr(2:2:ordCC,1:2)=[ratio/2 ...
      cumprod((ratio.^2-(1:2:(ordCC-2))'.^2)/4)./ ...
      factorial(2*(1:(ordCC/2-1))'+1).*ratio/2];
  else %disp('normal, not staggered, macro-grid')
    ks=1:2:ordCC; ps=(1:ordCC/2)'-1; 
    patches.Cwtsr(ks  ,1:2)=cumprod(ratio.^2-ps.^2,1) ...
        ./factorial(2*ps+1)./ratio;
    patches.Cwtsr(ks+1,1:2)=cumprod(ratio.^2-ps.^2,1) ...
        ./factorial(2*ps+2);
  end
  patches.Cwtsl = (-1).^((1:ordCC)'-patches.alt).*patches.Cwtsr;
end% if ordCC>0
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
\cdot \verb|DX|\) and~\(\verb|ratio(2)| \cdot \verb|DY|\),
unless \verb|patches.EdgyInt| is set in which case the
patches are of length \verb|ratio*DX+dx| and
\verb|ratio*DY+dy|.
\begin{matlab}
%}
nSubP = reshape(nSubP,1,2); % force to be row vector
assert(patches.EdgyInt | all(mod(nSubP,2)==1), ...
    'configPatches2: nSubP must be odd')
i0 = (nSubP(1)+1)/2;
if ~patches.EdgyInt, dx = ratio(1)*DX/(i0-1);
else                 dx = ratio(1)*DX/(nSubP(1)-2);
end
patches.x = dx*(-i0+1:i0-1)'+X;  % micro-grid
patches.x = reshape(patches.x,nSubP(1),1,1,1,nPatch(1),1);
i0 = (nSubP(2)+1)/2;
if ~patches.EdgyInt, dy = ratio(2)*DY/(i0-1);
else                 dy = ratio(2)*DY/(nSubP(2)-2);
end
patches.y = dy*(-i0+1:i0-1)'+Y;  % micro-grid
patches.y = reshape(patches.y,1,nSubP(2),1,1,1,nPatch(2));
%{
\end{matlab}

\paragraph{Set ensemble inter-patch communication}
%For \verb|EdgyInt| or centre interpolation respectively, the
%right-edge/centre realisations \verb|1:nEnsem| are to
%interpolate to left-edge~\verb|le|, and the left-edge/centre
%realisations \verb|1:nEnsem| are to interpolate
%to~\verb|re|. \verb|re| and \verb|li| are transposes or each
%other as \verb|re(li)=le(ri)| are both \verb|1:nEnsem|. One
%may use the statement
%\begin{verbatim}
%c=toeplitz(c(1:nSubP-1),c([1 nEnsem:-1:2]));
%\end{verbatim}
%to \emph{correspondingly} generates all phase shifted copies
%of microscale heterogeneity (see \verb|homoDiffEdgy1| of
%\cref{sec:homoDiffEdgy1}).
%\begin{matlab}
%}
nE=patches.nEnsem;
patches.le=1:nE;  patches.ri=1:nE;  % nothing shifty for now??
patches.bo=1:nE;  patches.to=1:nE; 
%if patches.EdgyInt, nP=nSubP-2;
%   else nP=(nSubP-1)/2; end
%patches.le=mod((0:nE-1)-mod(nP,nE),nE)+1;
%patches.ri=mod((0:nE-1)+mod(nP,nE),nE)+1;


end% function
%{
\end{matlab}
Fin.
\end{devMan}
%}
