% Creates a data struct of the design of patches for later
% use by the patch functions such as smoothPatch1().  
% AJR, Nov 2017 -- Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches1()}: configures spatial
patches in 1D}
\label{sec:configPatches1}
\localtableofcontents

\subsection{Introduction}

Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative\slash step function
\verb|patchSmooth1()|. \cref{sec:configPatches1eg} lists an
example of its use.
\begin{matlab}
%}
function configPatches1(fun,Xlim,BCs,nPatch,ordCC,ratio,nSubP ...
                       ,varargin)
global patches
%{
\end{matlab}

\paragraph{Input}
If invoked with no input arguments, then executes an example
of simulating Burgers' \pde---see \cref{sec:configPatches1eg}
for the example code.
\begin{itemize}

\item \verb|fun| is the name of the user function,
\verb|fun(t,u,x)|, that computes time derivatives (or
time-steps) of quantities on the patches.

\item \verb|Xlim| give the macro-space spatial domain of the
computation: patches are equi-spaced over the interior of
the interval~\([\verb|Xlim(1)|,\verb|Xlim(2)|]\).

\item \verb|BCs| somehow will define the macroscale boundary
conditions. Currently, \verb|BCs| is ignored and the system
is assumed macro-periodic in the spatial domain.

\item \verb|nPatch| is the number of equi-spaced spaced
patches.

\item \verb|ordCC|, must be~\(\geq -1\), is the `order' of
interpolation across empty space of the macroscale patch
values to the edge of the patches for inter-patch coupling:
where \verb|ordCC| of~\(0\) or~\(-1\) gives spectral
interpolation; and \verb|ordCC| being odd is for staggered
spatial grids.

\item \verb|ratio| (real) is the ratio of (depending upon
\verb|EdgyInt|) either the half-width or full-width of a
patch to the spacing of the patch mid-points.  So either
\(\verb|ratio|=\tfrac12\) means the patches abut and
\(\verb|ratio|=1\) is overlapping patches as in holistic
discretisation, or \(\verb|ratio|=1\) means the patches
abut.  Small~\verb|ratio| should greatly reduce
computational time.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch. If not using \verb|EdgyInt|,
then must be odd so that there is a central lattice point.

\item \verb|'nEdge'| (not yet implemented), \emph{optional},
default=1, for each patch, the number of edge values set by
interpolation at the edge regions of each patch.  The
default is one (suitable for microscale lattices with only
nearest neighbour interactions).

\item \verb|'EdgyInt'|, true/false, \emph{optional},
default=false.  If true, then interpolate to left\slash right
edge-values from right\slash left next-to-edge values.  

\item \verb|'nEnsem'|,  \emph{optional-experimental},
default one, but if more, then an ensemble over this
number of realisations.

\item \verb|'hetCoeffs'|, \emph{optional}, default empty.
Supply an array of microscale heterogeneous coefficients to
be used by the given microscale \verb|fun| in each patch.
Say the given array~\verb|cs| is of size \(m_x\times n_c\),
where \(n_c\)~is the number of different sets of
coefficients. The coefficients are to be the same for each
and every patch; however, macroscale variations are catered
for by the \(n_c\)~coefficients being \(n_c\)~parameters in
some macroscale formula.
\begin{itemize}
\item If \(\verb|nEnsem|=1\), then the array of coefficients
is just tiled across the patch size to fill up each patch,
starting from the first element.
\item If \(\verb|nEnsem|>1\) (value immaterial), then reset
\(\verb|nEnsem|:=m_x\) and construct an ensemble of all
\(m_x\) phase-shifts of the coefficients. In this scenario,
the inter-patch coupling couples different members in the
ensemble. When \verb|EdgyInt| is true, and when the
coefficients are diffusivities\slash elasticities, then this
coupling cunningly preserves symmetry .
\end{itemize}

\item \verb|'nCore'|,  \emph{optional-experimental}, default
one, but if more, and only for non-EdgyInt, then
interpolates from an average over the core of a patch, a
core of size ??.  Then edge values are set according to
interpolation of the averages?? or so that average at edges
is the interpolant??

\end{itemize}


\paragraph{Output} The \emph{global} struct \verb|patches|
is created and set with the following components.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x)| that computes the time derivatives (or
steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.stag| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the
\(\verb|ordCC|\)-vector of weights for the inter-patch
interpolation onto the right and left edges (respectively)
with patch:macroscale ratio as specified.

\item \verb|.x| (4D) is \(\verb|nSubP| \times1 \times1
\times \verb|nPatch|\) array of the regular spatial
locations~\(x_{ij}\) of the \(i\)th~microscale grid point in
the \(j\)th~patch.  

\item \verb|.nEdge| is, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch.

\item \verb|.le|, \verb|.ri|
determine inter-patch coupling of members in an ensemble.
Each a column vector of length~\verb|nEnsem|.

\item \verb|.cs| either
\begin{itemize}
\item \verb|[]| 0D, or 
\item if \(\verb|nEnsem|=1\), \((\verb|nSubP(1)|-1)\times
n_c\) 2D array of microscale heterogeneous coefficients, or
\item if \(\verb|nEnsem|>1\), \((\verb|nSubP(1)|-1) \times
n_c\times m_x\) 3D array of \(m_x\)~ensemble of phase-shifts
of the microscale
heterogeneous coefficients.
\end{itemize}

\end{itemize}




\subsection{If no arguments, then execute an example}
\label{sec:configPatches1eg}
\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
The code here shows one way to get started: a user's script
may have the following three steps (left-right arrows denote
function recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches1 
\item ode15s integrator \into patchSmooth1 \into user's PDE
\item process results
\end{enumerate}
Establish global patch data struct to point to and interface
with a function coding Burgers' \pde: to be solved on
\(2\pi\)-periodic domain, with eight patches, spectral
interpolation couples the patches, each patch of half-size
ratio~\(0.2\), and with seven microscale points forming each
patch.
\begin{matlab}
%}
configPatches1(@BurgersPDE,[0 2*pi], nan, 8, 0, 0.2, 7);
%{
\end{matlab}
Set an initial condition, with some microscale randomness.
\begin{matlab}
%}
u0=0.3*(1+sin(patches.x))+0.1*randn(size(patches.x));
%{
\end{matlab}
Simulate in time using a standard stiff integrator and the
interface function \verb|patchsmooth1()|
(\cref{sec:patchSmooth1}).
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
[ts,us] = ode15s( @patchSmooth1,[0 0.5],u0(:));
else % octave version
[ts,us] = odeOcts(@patchSmooth1,[0 0.5],u0(:));
end
%{
\end{matlab}
Plot the simulation using only the microscale values
interior to the patches: either set \(x\)-edges to
\verb|nan| to leave the gaps; or use \verb|patchEdgyInt1| to
re-interpolate correct patch edge values and thereby join
the patches.  \autoref{fig:config1Burgers} illustrates an
example simulation in time generated by the patch scheme
applied to Burgers'~\pde.
\begin{matlab}
%}
figure(1),clf
if 1, patches.x([1 end],:,:,:)=nan;  us=us.';
else us=reshape(patchEdgyInt1(us.'),[],length(ts));  
end
surf(ts,patches.x(:),us)
view(60,40), colormap(hsv)
title('Burgers PDE: patches in space, continuous time')
xlabel('time t'), ylabel('space x'), zlabel('u(x,t)')
ifOurCf2eps(mfilename)
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:config1Burgers}field
\(u(x,t)\) of the patch scheme applied to Burgers'~\pde.}
\includegraphics[scale=0.85]{configPatches1}
\end{figure}
Upon finishing execution of the example, exit this function.
\begin{matlab}
%}
return
end%if no arguments
%{
\end{matlab}

\input{../Patch/BurgersPDE.m}
\input{../Patch/odeOcts.m}





\begin{devMan}

\subsection{Parse input arguments and defaults}
\begin{matlab}
%}
p = inputParser;
fnValidation = @(f) isa(f, 'function_handle'); %test for a function name
addRequired(p,'fun',fnValidation); 
addRequired(p,'Xlim',@isnumeric);
addRequired(p,'BCs'); % nothing yet decided
addRequired(p,'nPatch',@isnumeric);
addRequired(p,'ordCC',@isnumeric);
addRequired(p,'ratio',@isnumeric);
addRequired(p,'nSubP',@isnumeric);
addParameter(p,'nEdge',1,@isnumeric);
addParameter(p,'EdgyInt',false,@islogical);
addParameter(p,'nEnsem',1,@isnumeric);
addParameter(p,'nCore',1,@isnumeric);
addParameter(p,'hetCoeffs',[],@isnumeric);
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
patches.nCore = p.Results.nCore;
%{
\end{matlab}
Check parameters.
\begin{matlab}
%}
assert(patches.nEdge==1 ...
      ,'multi-edge-value interp not yet implemented')
assert(2*patches.nEdge+1<=nSubP ...
      ,'too many edge values requested')
if patches.nCore>1
    warning('nCore>1 not yet tested in this version')
    end
%{
\end{matlab}




\subsection{The code to make patches and interpolation}
First, store the pointer to the time derivative function in
the struct.
\begin{matlab}
%}
patches.fun=fun;
%{
\end{matlab}

Second, store the order of interpolation that is to provide
the values for the inter-patch coupling conditions. Spectral
coupling is \verb|ordCC| of~\(0\) and~\(-1\).
\begin{matlab}
%}
assert((ordCC>=-1) & (floor(ordCC)==ordCC), ...
    'ordCC out of allowed range integer>=-1')
%{
\end{matlab}
For odd~\verb|ordCC|, interpolate based upon odd
neighbouring patches as is useful for staggered grids.
\begin{matlab}
%}
patches.stag=mod(ordCC,2);
ordCC=ordCC+patches.stag;
patches.ordCC=ordCC;
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
interpolation of field values for coupling. (Could sometime
extend to coupling via derivative values.)
\begin{matlab}
%}
if ordCC>0, patchCwts(ratio,ordCC,patches.stag), end
%{
\end{matlab}


Third, set the centre of the patches in a the macroscale
grid of patches assuming periodic macroscale domain.
\begin{matlab}
%}
X=linspace(Xlim(1),Xlim(2),nPatch+1);
X=X(1:nPatch)+diff(X)/2;
DX=X(2)-X(1);
%{
\end{matlab}
Construct the microscale grid in each patch, assuming
Dirichlet patch edges, and a half-patch length of
\(\verb|ratio| \cdot \verb|DX|\), unless
\verb|patches.EdgyInt| is set in which case the patches are
of length \verb|ratio*DX+dx|. Reshape the grid to be 4D to
suit dimensions (micro,Vars,Ens,macro).
\begin{matlab}
%}
assert(patches.EdgyInt | mod(nSubP,2)==1, ...
    'configPatches1: nSubP must be odd')
i0=(nSubP+1)/2;
if ~patches.EdgyInt, dx = ratio*DX/(i0-1);  
else                 dx = ratio*DX/(nSubP-2);
end
patches.x = dx*(-i0+1:i0-1)'+X; % micro-grid
patches.x = reshape(patches.x,nSubP,1,1,nPatch);
%{
\end{matlab}




\subsection{Set ensemble inter-patch communication}
For \verb|EdgyInt| or centre interpolation respectively, the
right-edge\slash centre realisations \verb|1:nEnsem| are to
interpolate to left-edge~\verb|le|, and the left-edge\slash
centre realisations \verb|1:nEnsem| are to interpolate
to~\verb|re|. \verb|re| and \verb|li| are `transposes' of
each other as \verb|re(li)=le(ri)| are both \verb|1:nEnsem|.
Alternatively, one may use the statement
\begin{verbatim}
c=hankel(c(1:nSubP-1),c([nSubP 1:nSubP-2]));
\end{verbatim}
to \emph{correspondingly} generates all phase shifted copies
of microscale heterogeneity (see \verb|homoDiffEdgy1| of
\cref{sec:homoDiffEdgy1}).

The default is nothing shifty.  This setting reduces the
number of if-statements in function \verb|patchEdgeInt2()|.
\begin{matlab}
%}
nE = patches.nEnsem;
patches.le = 1:nE;  
patches.ri = 1:nE;  
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
  [mx,nc] = size(cs);
  nx = nSubP(1); 
  cs = repmat(cs,nSubP,1);
%{
\end{matlab}
If only one member of the ensemble is required, then
sub-sample to patch size, and store coefficients in
\verb|patches| as is.
\begin{matlab}
%}
  if nE==1, patches.cs = cs(1:nx-1,:); else
%{
\end{matlab}
But for $\verb|nEnsem|>1$ an ensemble of
\(m_x\)~phase-shifts of the coefficients is constructed
from the over-supply.  Here code phase-shifts over the
periods---the phase shifts are like Hankel-matrices.
\begin{matlab}
%}
    patches.nEnsem = mx;
    patches.cs = nan(nx-1,nc,mx);
    for i = 1:mx
        is = (i:i+nx-2);
        patches.cs(:,:,i) = cs(is,:);
    end
    patches.cs = reshape(patches.cs,nx-1,nc,[]);
%{
\end{matlab}
Further, set a cunning left\slash right
realisation of inter-patch coupling.  The aim is to
preserve symmetry in the system when also invoking
\verb|EdgyInt|.  What this coupling does without
\verb|EdgyInt| is unknown.  Use auto-replication.
\begin{matlab}
%}
    patches.le = mod((0:mx-1)'+mod(nx-2,mx),mx)+1;
    patches.ri = mod((0:mx-1)'-mod(nx-2,mx),mx)+1;
%{
\end{matlab}
Issue warning if the ensemble is likely to be affected by
lack of scale separation.  Need to justify this and the
arbitrary threshold more carefully??
\begin{matlab}
%}
if ratio*patches.nEnsem>0.9, warning( ...
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