% Creates a data struct of the design of 2D patches for
% later use by the patch functions such as smoothPatch2() 
% AJR, Nov 2018 -- Jul 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches2()}: configures spatial
patches in 2D}
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
patch values to the edge-values of the patches: currently
must be~\(0,2,4,\ldots\); where \(0\)~gives spectral
interpolation.


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
left\slash bottom\slash top next-to-edge values.  

\item \verb|'nEnsem'|,  \emph{optional-experimental},
default one, but if more, then an ensemble over this
number of realisations.

\item \verb|'hetCoeffs'|, \emph{optional}, default empty.
Supply an array of microscale heterogeneous coefficients to
be used by the given microscale \verb|fun| in each patch.
Say the given array~\verb|cs| is of size \(m_x\times
m_y\times n_c\), where \(n_c\)~is the number of different
sets of coefficients. For example, in heterogeneous
diffusion, \(n_c=2\) for the diffusivities in the \emph{two}
different spatial directions. The coefficients are to
be the same for each and every patch; however, macroscale
variations are catered for by the \(n_c\)~coefficients being
\(n_c\)~parameters in some macroscale formula.
\begin{itemize}

\item If \(\verb|nEnsem|=1\), then the array of coefficients
is just tiled across the patch size to fill up each patch,
starting from the \((1,1)\)-element.

\item If \(\verb|nEnsem|>1\) (value immaterial), then reset
\(\verb|nEnsem|:=m_x\cdot m_y\) and construct an ensemble of
all \(m_x\cdot m_y\) phase-shifts of the coefficients. In
this scenario, the inter-patch coupling couples different
members in the ensemble. When \verb|EdgyInt| is true, and
when the coefficients are diffusivities\slash elasticities
in~\(x\) and~\(y\) directions, respectively, then this
coupling cunningly preserves symmetry.

\end{itemize}

\end{itemize}



\paragraph{Output} The \emph{global} struct \verb|patches|
is created and set with the following components.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x,y)| that computes the time derivatives (or
steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.stag| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling---not yet
implemented.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the
\(\verb|ordCC|\times 2\)-array of weights for the inter-patch
interpolation onto the right\slash top and left\slash bottom
edges (respectively) with patch:macroscale ratio as
specified.

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

\item \verb|.le|, \verb|.ri|, \verb|.bo|, \verb|.to|
determine inter-patch coupling of members in an ensemble.
Each a column vector of length~\verb|nEnsem|.

\item \verb|.cs| either
\begin{itemize}
\item \verb|[]| 0D, or 
\item if \(\verb|nEnsem|=1\), \((\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times n_c\) 3D array of microscale
heterogeneous coefficients, or
\item if \(\verb|nEnsem|>1\), \((\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times n_c\times m_xm_y\) 4D array of
\(m_xm_y\)~ensemble of phase-shifts of the microscale
heterogeneous coefficients.
\end{itemize}

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
\item ode23 integrator \into patchSmooth2 \into user's PDE
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
    , 0, 0.4, nSubP ,'EdgyInt',false);
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
ifOurCf2eps([mfilename 'ic'])
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
disp('Wait to simulate nonlinear diffusion h_t=(h^3)_xx+(h^3)_yy')
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
ifOurCf2eps([mfilename 't3'])
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
fnValidation = @(f) isa(f, 'function_handle');
addRequired(p,'fun',fnValidation); % is there test for a function name??
addRequired(p,'Xlim',@isnumeric);
addRequired(p,'BCs'); % nothing yet decided
addRequired(p,'nPatch',@isnumeric);
addRequired(p,'ordCC',@isnumeric);
addRequired(p,'ratio',@isnumeric);
addRequired(p,'nSubP',@isnumeric);
addParameter(p,'nEdge',1,@isnumeric);
addParameter(p,'EdgyInt',false,@islogical);
addParameter(p,'nEnsem',1,@isnumeric);
addParameter(p,'hetCoeffs',[],@isnumeric);
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
cs = p.Results.hetCoeffs;
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
patches.stag = mod(ordCC,2);
assert(patches.stag==0,'staggered not yet implemented??')
ordCC = ordCC+patches.stag;
patches.ordCC = ordCC;
%{
\end{matlab}
Check for staggered grid and periodic case.
\begin{matlab}
%}
if patches.stag, assert(mod(nPatch,2)==0, ...
    'Require an even number of patches for staggered grid')
end
%{
\end{matlab}
Might as well precompute the weightings for the
interpolation of field values for coupling. 
(Could sometime extend to coupling via derivative
values.)
\begin{matlab}
%}
ratio = reshape(ratio,1,2); % force to be row vector
if ordCC>0, patchCwts(ratio,ordCC,patches.stag), end
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





\subsection{Set ensemble inter-patch communication}
For \verb|EdgyInt| or centre interpolation respectively, the
right-edge\slash centre realisations \verb|1:nEnsem| are to
interpolate to left-edge~\verb|le|, and the left-edge\slash
centre realisations \verb|1:nEnsem| are to interpolate
to~\verb|re|. \verb|re| and \verb|li| are `transposes' of
each other as \verb|re(li)=le(ri)| are both \verb|1:nEnsem|.
Similarly for bottom-edge\slash centre interpolation to
top-edge via~\verb|to| and top-edge\slash centre
interpolation to bottom-edge via~\verb|bo|.

The default is nothing shifty.  This setting reduces the
number of if-statements in function \verb|patchEdgeInt2()|.
\begin{matlab}
%}
nE = patches.nEnsem;
patches.le = 1:nE;  patches.ri = 1:nE;  
patches.bo = 1:nE;  patches.to = 1:nE; 
%{
\end{matlab}
However, if heterogeneous coefficients are supplied via
\verb|hetCoeffs|, then do some non-trivial replications. 
First, get microscale periods, patch size, and replicate
many times in order to subsequently sub-sample: 
\verb|nSubP| times should be enough.
\begin{matlab}
%}
if ~isempty(cs)
  [mx,my,nc] = size(cs);
  nx = nSubP(1); ny = nSubP(2);
  cs = repmat(cs,nSubP);
%{
\end{matlab}
If only one member of the ensemble is required, then
sub-sample to patch size, and store coefficients in
\verb|patches| as is.
\begin{matlab}
%}
  if nE==1, patches.cs = cs(1:nx-1,1:ny-1,:); else
%{
\end{matlab}
But for $\verb|nEnsem|>1$ an ensemble of
\(m_xm_y\)~phase-shifts of the coefficients is constructed
from the over-supply.  Here code phase-shifts over the
periods---the phase shifts are like Hankel-matrices.
\begin{matlab}
%}
    patches.nEnsem = mx*my;
    patches.cs = nan(nx-1,ny-1,nc,mx,my);
    for j = 1:my
        js = (j:j+ny-2);
        for i = 1:mx
            is = (i:i+nx-2);
            patches.cs(:,:,:,i,j) = cs(is,js,:);
        end
    end
    patches.cs = reshape(patches.cs,nx-1,ny-1,nc,[]);
%{
\end{matlab}
Further, set a cunning left\slash right\slash bottom\slash
top realisation of inter-patch coupling.  The aim is to
preserve symmetry in the system when also invoking
\verb|EdgyInt|.  What this coupling does without
\verb|EdgyInt| is unknown.  Use auto-replication.
\begin{matlab}
%}
    le = mod((0:mx-1)+mod(nx-2,mx),mx)+1;
    patches.le = reshape(  le'+mx*(0:my-1)  ,[],1);
    ri = mod((0:mx-1)-mod(nx-2,mx),mx)+1;
    patches.ri = reshape(  ri'+mx*(0:my-1)  ,[],1);
    bo = mod((0:my-1)+mod(ny-2,my),my)+1;
    patches.bo = reshape( (1:mx)'+mx*(bo-1) ,[],1);
    to = mod((0:my-1)-mod(ny-2,my),my)+1;
    patches.to = reshape( (1:mx)'+mx*(to-1) ,[],1);
%{
\end{matlab}
Issue warning if the ensemble is likely to be affected by
lack of scale separation.  Need to justify this and the
arbitrary threshold more carefully??
\begin{matlab}
%}
if prod(ratio)*patches.nEnsem>0.9, warning( ...
'Probably poor scale separation in ensemble of coupled phase-shifts')
scaleSeparationParameter=ratio*patches.nEnsem
end
%{
\end{matlab}
End the two if-statements.
\begin{matlab}
%}
  end%if nEnsem=1 or >1  
end%if not-empty(cs)
%{
\end{matlab}


\paragraph{Fin}
\begin{matlab}
%}
end% function
%{
\end{matlab}
\end{devMan}
%}
