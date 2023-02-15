% configPatches2() creates a data struct of the design of 2D
% patches for later use by the patch functions such as
% patchSys2().  AJR, Nov 2018 -- Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches2()}: configures spatial
patches in 2D}
\label{sec:configPatches2}
\localtableofcontents



Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative\slash step function
\verb|patchSys2()|. \cref{sec:configPatches2eg} lists an
example of its use.
\begin{matlab}
%}
function patches = configPatches2(fun,Xlim,Dom ...
    ,nPatch,ordCC,dx,nSubP,varargin)
%{
\end{matlab}



\paragraph{Input}
If invoked with no input arguments, then executes an example
of simulating a nonlinear diffusion \pde\ relevant to the
lubrication flow of a thin layer of fluid---see
\cref{sec:configPatches2eg} for an example code.
\begin{itemize}

\item \verb|fun| is the name of the user function,
\verb|fun(t,u,patches)| or \verb|fun(t,u)|, that computes
time-derivatives (or time-steps) of quantities on the 2D
micro-grid within all the 2D~patches.

\item \verb|Xlim| array/vector giving the rectangular
macro-space domain of the computation, namely
$[\verb|Xlim(1)|, \verb|Xlim(2)|] \times [\verb|Xlim(3)|,
\verb|Xlim(4)|]$. If \verb|Xlim| has two elements, then the
domain is the square domain of the same interval in both
directions.

\item \verb|Dom| sets the type of macroscale conditions for
the patches, and reflects the type of microscale boundary
conditions of the problem.   If \verb|Dom| is \verb|NaN| or
\verb|[]|, then the field~\verb|u| is doubly macro-periodic
in the 2D spatial domain, and resolved on equi-spaced
patches. If \verb|Dom| is a character string, then that
specifies the \verb|.type| of the following structure, with
\verb|.bcOffset| set to the default zero. Otherwise
\verb|Dom| is a structure with the following components.
\begin{itemize}

\item \verb|.type|, string, of either \verb|'periodic'| (the
default), \verb|'equispace'|, \verb|'chebyshev'|,
\verb|'usergiven'|.  For all cases except \verb|'periodic'|,
users \emph{must} code into \verb|fun| the micro-grid
boundary conditions that apply at the left\slash right\slash
bottom\slash top edges of the leftmost\slash rightmost\slash
bottommost\slash topmost patches, respectively.

\item \verb|.bcOffset|, optional one, two or four element
vector/array, in the cases of \verb|'equispace'| or
\verb|'chebyshev'| the patches are placed so the left\slash
right\slash top\slash bottom macroscale boundaries are
aligned to the left\slash right\slash top\slash bottom edges
of the corresponding extreme patches, but offset by
\verb|.bcOffset| of the sub-patch micro-grid spacing.  For
example, use \verb|bcOffset=0| when the micro-code applies
Dirichlet boundary values on the extreme edge micro-grid
points, whereas use \verb|bcOffset=0.5| when the microcode
applies Neumann boundary conditions halfway between the
extreme edge micro-grid points.  Similarly for the top and
bottom edges.

If \verb|.bcOffset| is a scalar, then apply the same offset
to all boundaries. If two elements, then apply the first
offset to both \(x\)-boundaries, and the second offset to
both \(y\)-boundaries. If four elements, then apply the
first two offsets to the respective \(x\)-boundaries, and
the last two offsets to the respective \(y\)-boundaries.

\item \verb|.X|, optional vector/array with \verb|nPatch(1)|
elements, in the case \verb|'usergiven'| it specifies the
\(x\)-locations of the centres of the patches---the user is
responsible the locations makes sense.

\item \verb|.Y|, optional vector/array with \verb|nPatch(2)|
elements, in the case \verb|'usergiven'| it specifies the
\(y\)-locations of the centres of the patches---the user is
responsible the locations makes sense.
\end{itemize}


\item \verb|nPatch| sets the number of equi-spaced spatial
patches: if scalar, then use the same number of patches in
both directions, otherwise \verb|nPatch(1:2)| gives the
number of patches~($\geq1$) in each direction.

\item \verb|ordCC| is the `order' of interpolation for
inter-patch coupling across empty space of the macroscale
patch values to the edge-values of the patches: currently
must be~$0,2,4,\ldots$; where $0$~gives spectral
interpolation.

\item \verb|dx| (real---scalar or two element) is usually
the sub-patch micro-grid spacing in~\(x\) and~\(y\).  If
scalar, then use the same \verb|dx| in both directions,
otherwise \verb|dx(1:2)| gives the spacing in each of the
two directions.

However, if \verb|Dom| is~\verb|NaN| (as for pre-2023), then
\verb|dx| actually is \verb|ratio| (scalar or two element),
namely the ratio of (depending upon \verb|EdgyInt|) either
the half-width or full-width of a patch to the equi-spacing
of the patch mid-points.  So either $\verb|ratio|=\tfrac12$
means the patches abut and $\verb|ratio|=1$ is overlapping
patches as in holistic discretisation, or $\verb|ratio|=1$
means the patches abut.  Small~\verb|ratio| should greatly
reduce computational time.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch: if scalar, then use the same
number in both directions, otherwise \verb|nSubP(1:2)| gives
the number in each direction. If not using \verb|EdgyInt|,
then must be odd so that there is/are centre-patch
micro-grid point\slash lines in each patch.

\item \verb|nEdge| (not yet implemented), \emph{optional},
default=1, for each patch, the number of edge values set by
interpolation at the edge regions of each patch.  The
default is one (suitable for microscale lattices with only
nearest neighbour interactions).

\item \verb|EdgyInt|, true/false, \emph{optional},
default=false.  If true, then interpolate to left\slash
right\slash top\slash bottom edge-values from right\slash
left\slash bottom\slash top next-to-edge values.  If false
or omitted, then interpolate from centre cross-patch lines.

\item \verb|nEnsem|,  \emph{optional-experimental},
default one, but if more, then an ensemble over this
number of realisations.

\item \verb|hetCoeffs|, \emph{optional}, default empty.
Supply a 2D or 3D array of microscale heterogeneous
coefficients to be used by the given microscale \verb|fun|
in each patch. Say the given array~\verb|cs| is of size
$m_x\times m_y\times n_c$, where $n_c$~is the number of
different sets of coefficients.  For example, in
heterogeneous diffusion, $n_c=2$ for the diffusivities in
the \emph{two} different spatial directions (or $n_c=3$ for
the diffusivity tensor). The coefficients are to be the same
for each and every patch; however, macroscale variations are
catered for by the $n_c$~coefficients being $n_c$~parameters
in some macroscale formula.
\begin{itemize}
\item If $\verb|nEnsem|=1$, then the array of coefficients
is just tiled across the patch size to fill up each patch,
starting from the $(1,1)$-point in each patch.

\item If $\verb|nEnsem|>1$ (value immaterial), then reset
$\verb|nEnsem|:=m_x\cdot m_y$ and construct an ensemble of
all $m_x\cdot m_y$ phase-shifts of the coefficients. In
this scenario, the inter-patch coupling couples different
members in the ensemble. When \verb|EdgyInt| is true, and
when the coefficients are diffusivities\slash elasticities
in~$x$ and~$y$ directions, respectively, then this
coupling cunningly preserves symmetry.

\end{itemize}

\item \verb|'parallel'|, true/false, \emph{optional},
default=false. If false, then all patch computations are on
the user's main \textsc{cpu}---although a user may well
separately invoke, say, a \textsc{gpu} to accelerate
sub-patch computations. 

If true, and it requires that you have \Matlab's Parallel
Computing Toolbox, then it will distribute the patches over
multiple \textsc{cpu}s\slash cores.  In \Matlab, only one
array dimension can be split in the distribution, so it
chooses the one space dimension~$x,y$ corresponding to the
highest~\verb|\nPatch| (if a tie, then chooses the rightmost
of~$x,y$).  A user may correspondingly distribute arrays
with property \verb|patches.codist|, or simply use formulas
invoking the preset distributed arrays \verb|patches.x|, and
\verb|patches.y|.  If a user has not yet established a
parallel pool, then a `local' pool is started.

\end{itemize}



\paragraph{Output} The struct \verb|patches| is created and
set with the following components.  If no output variable is
provided for \verb|patches|, then make the struct available
as a global variable.\footnote{When using \texttt{spmd}
parallel computing, it is generally best to avoid global
variables, and so instead prefer using an explicit output
variable.}
\begin{matlab}
%}
if nargout==0, global patches, end
%{
\end{matlab}
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,patches)| or \verb|fun(t,u)|, that computes
the time derivatives (or steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.periodic|: either true, for interpolation on
the macro-periodic domain; or false, for general
interpolation by divided differences over non-periodic
domain or unevenly distributed patches.

\item \verb|.stag| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling---not yet
implemented.

\item \verb|.Cwtsr| and \verb|.Cwtsl|, only for
macro-periodic conditions, are the
$\verb|ordCC|\times 2$-array of weights for the inter-patch
interpolation onto the right\slash top and left\slash bottom
edges (respectively) with patch:macroscale ratio as
specified or as derived from~\verb|dx|.

\item \verb|.x| (6D) is $\verb|nSubP(1)| \times1 \times1
\times1 \times \verb|nPatch(1)| \times1$ array of the
regular spatial locations~$x_{iI}$ of the microscale grid
points in every patch.  

\item \verb|.y| (6D) is $1 \times \verb|nSubP(2)| \times1
\times1 \times1 \times \verb|nPatch(2)|$ array of the
regular spatial locations~$y_{jJ}$ of the microscale grid
points in every patch.  

\item \verb|.ratio| $1\times 2$, only for
macro-periodic conditions, are the size ratios of
every patch.

\item \verb|.nEdge| is, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch.

\item \verb|.le|, \verb|.ri|, \verb|.bo|, \verb|.to|
determine inter-patch coupling of members in an ensemble.
Each a column vector of length~\verb|nEnsem|.

\item \verb|.cs| either
\begin{itemize}
\item \verb|[]| 0D, or 
\item if $\verb|nEnsem|=1$, $(\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times n_c$ 3D array of microscale
heterogeneous coefficients, or
\item if $\verb|nEnsem|>1$, $(\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times n_c\times m_xm_y$ 4D array of
$m_xm_y$~ensemble of phase-shifts of the microscale
heterogeneous coefficients.
\end{itemize}

\item \verb|.parallel|, logical: true if patches are
distributed over multiple \textsc{cpu}s\slash cores for the
Parallel Computing Toolbox, otherwise false (the default is
to activate the \emph{local} pool).

\item \verb|.codist|, \emph{optional}, describes the
particular parallel distribution of arrays over the active
parallel pool.  

\end{itemize}






\subsection{If no arguments, then execute an example}
\label{sec:configPatches2eg}
\begin{matlab}
%}
if nargin==0
disp('With no arguments, simulate example of nonlinear diffusion')
%{
\end{matlab}
The code here shows one way to get started: a user's script
may have the following three steps (``\into'' denotes
function recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches2 
\item ode23 integrator \into patchSys2 \into user's PDE
\item process results
\end{enumerate}

Establish global patch data struct to interface with a
function coding a nonlinear `diffusion' \pde: to be solved
on $6\times4$-periodic domain, with $9\times7$ patches,
spectral interpolation~($0$) couples the patches, with
$5\times5$ points forming the micro-grid in each patch, and
a sub-patch micro-grid spacing of~\(0.12\) (relatively large
for visualisation).  \cite{Roberts2011a} established that
this scheme is consistent with the \pde\ (as the patch
spacing decreases).
\begin{matlab}
%}
global patches
patches = configPatches2(@nonDiffPDE,[-3 3 -2 2], [] ...
    , [9 7], 0, 0.12, 5 ,'EdgyInt',false);
%{
\end{matlab}
Set an  initial condition of a perturbed-Gaussian using
auto-replication of the spatial grid.
\begin{matlab}
%}
u0 = exp(-patches.x.^2-patches.y.^2);
u0 = u0.*(0.9+0.1*rand(size(u0)));
%{
\end{matlab}
Initiate a plot of the simulation using only the microscale
values interior to the patches: optionally set $x$~and
$y$-edges to \verb|nan| to leave the gaps between patches.
\begin{matlab}
%}
figure(1), clf, colormap(0.8*hsv)
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
ifOurCf2eps([mfilename 'ic'])
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:configPatches2ic}initial
field~$u(x,y,t)$ at time $t=0$ of the patch scheme applied
to a nonlinear `diffusion'~\pde: \cref{fig:configPatches2t3}
plots the computed field at time $t=3$.}
\includegraphics[scale=0.9]{configPatches2ic}
\end{figure}
Integrate in time to $t=4$ using standard functions. In
\Matlab\ \verb|ode15s| would be natural as the patch scheme
is naturally stiff, but \verb|ode23| is quicker \cite
[Fig.~4] {Maclean2020a}.  Ask for output at non-uniform
times because the diffusion slows.
\begin{matlab}
%}
disp('Wait to simulate nonlinear diffusion h_t=(h^3)_xx+(h^3)_yy')
drawnow
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSys2,linspace(0,2).^2,u0(:));
else % octave version is quite slow for me
    lsode_options('absolute tolerance',1e-4);
    lsode_options('relative tolerance',1e-4);
    [ts,us] = odeOcts(@patchSys2,[0 1],u0(:));
end
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:configPatches2t3}.  Use \verb|patchEdgeInt2| to
interpolate patch-edge values.
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
\caption{\label{fig:configPatches2t3}field~$u(x,y,t)$ at
time $t=3$ of the patch scheme applied to a nonlinear
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


\IfFileExists{../Patch/nonDiffPDE.m}{\input{../Patch/nonDiffPDE.m}}{}



\begin{devMan}

\subsection{Parse input arguments and defaults}
\begin{matlab}
%}
p = inputParser;
fnValidation = @(f) isa(f, 'function_handle');%test for fn name
addRequired(p,'fun',fnValidation); 
addRequired(p,'Xlim',@isnumeric);
%addRequired(p,'Dom'); % nothing yet decided
addRequired(p,'nPatch',@isnumeric);
addRequired(p,'ordCC',@isnumeric);
addRequired(p,'dx',@isnumeric);
addRequired(p,'nSubP',@isnumeric);
addParameter(p,'nEdge',1,@isnumeric);
addParameter(p,'EdgyInt',false,@islogical);
addParameter(p,'nEnsem',1,@isnumeric);
addParameter(p,'hetCoeffs',[],@isnumeric);
addParameter(p,'parallel',false,@islogical);
%addParameter(p,'nCore',1,@isnumeric); % not yet implemented
parse(p,fun,Xlim,nPatch,ordCC,dx,nSubP,varargin{:});
%{
\end{matlab}
Set the optional parameters. 
\begin{matlab}
%}
patches.nEdge = p.Results.nEdge;
patches.EdgyInt = p.Results.EdgyInt;
patches.nEnsem = p.Results.nEnsem;
cs = p.Results.hetCoeffs;
patches.parallel = p.Results.parallel;
%patches.nCore = p.Results.nCore;
%{
\end{matlab}

Initially duplicate parameters for both space dimensions as
needed.
\begin{matlab}
%}
if numel(Xlim)==2,   Xlim = repmat(Xlim,1,2); end
if numel(nPatch)==1, nPatch = repmat(nPatch,1,2); end
if numel(dx)==1,     dx = repmat(dx,1,2); end
if numel(nSubP)==1,  nSubP = repmat(nSubP,1,2); end
%{
\end{matlab}

Check parameters.
\begin{matlab}
%}
assert(Xlim(1)<Xlim(2) ...
      ,'first pair of Xlim must be ordered increasing')
assert(Xlim(3)<Xlim(4) ...
      ,'second pair of Xlim must be ordered increasing')
assert(patches.nEdge==1 ...
      ,'multi-edge-value interp not yet implemented')
assert(all(2*patches.nEdge<nSubP) ...
      ,'too many edge values requested')
%if patches.nCore>1
%    warning('nCore>1 not yet tested in this version')
%    end
%{
\end{matlab}



For compatibility with pre-2023 functions, if parameter
\verb|Dom| is \verb|Nan|, then  we set the \verb|ratio| to
be the value of the so-called \verb|dx| vector.
\begin{matlab}
%}
if ~isstruct(Dom), pre2023=isnan(Dom);
else pre2023=false; end
if pre2023, ratio=dx; dx=nan; end
%{
\end{matlab}

Default macroscale conditions are periodic with evenly
spaced patches.
\begin{matlab}
%}
if isempty(Dom), Dom=struct('type','periodic'); end
if (~isstruct(Dom))&isnan(Dom), Dom=struct('type','periodic'); end
%{
\end{matlab}
If \verb|Dom| is a string, then just set type to that
string, and subsequently set corresponding defaults for
others fields.
\begin{matlab}
%}
if ischar(Dom), Dom=struct('type',Dom); end
%{
\end{matlab}
We allow different macroscale domain conditions in the
different directions. But for the moment do not allow
periodic to be mixed with the others (as the interpolation
mechanism is different code)---hence why we choose
\verb|periodic| be seven characters, whereas the others are
eight characters. The different conditions are coded in
different rows of  \verb|Dom.type|, so we duplicate the
string if only one row specified.
\begin{matlab}
%}
if size(Dom.type,1)==1, Dom.type=repmat(Dom.type,2,1); end
%{
\end{matlab}
Check what is and is not specified, and provide default of
zero 
(Dirichlet boundaries) if no \verb|bcOffset| specified when
needed.  Do so for both directions independently.
\begin{matlab}
%}
patches.periodic=false;
for p=1:2
switch Dom.type(p,:)
case 'periodic'
    patches.periodic=true;
    if isfield(Dom,'bcOffset')
    warning('bcOffset not available for Dom.type = periodic'), end
    msg=' not available for Dom.type = periodic';
    if isfield(Dom,'X'), warning(['X' msg]), end
    if isfield(Dom,'Y'), warning(['Y' msg]), end
case {'equispace','chebyshev'}
    if ~isfield(Dom,'bcOffset'), Dom.bcOffset=zeros(2,2); end
    % for mixed with usergiven, following should still work
    if numel(Dom.bcOffset)==1
        Dom.bcOffset=repmat(Dom.bcOffset,2,2); end
    if numel(Dom.bcOffset)==2
        Dom.bcOffset=repmat(Dom.bcOffset(:)',2,1); end
    msg=' not available for Dom.type = equispace or chebyshev';
    if (p==1)& isfield(Dom,'X'), warning(['X' msg]), end
    if (p==2)& isfield(Dom,'Y'), warning(['Y' msg]), end
case 'usergiven'
%    if isfield(Dom,'bcOffset')
%    warning('bcOffset not available for usergiven Dom.type'), end
    msg=' required for Dom.type = usergiven';
    if p==1, assert(isfield(Dom,'X'),['X' msg]), end
    if p==2, assert(isfield(Dom,'Y'),['Y' msg]), end
otherwise 
    error([Dom.type 'is unknown Dom.type'])
end%switch Dom.type
end%for p
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
coupling is \verb|ordCC| of~$0$ or (not yet??)~$-1$.
\todo{Perhaps implement staggered spectral coupling.}
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
if patches.stag, assert(all(mod(nPatch,2)==0), ...
    'Require an even number of patches for staggered grid')
end
%{
\end{matlab}



\paragraph{Set the macro-distribution of patches}
Third, set the centre of the patches in the macroscale grid
of patches.  Loop over the coordinate directions, setting
the distribution into~\verb|Q| and finally assigning to
array of corresponding direction.
\begin{matlab}
%}
for q=1:2
qq=2*q-1;
%{
\end{matlab}
Distribution depends upon \verb|Dom.type|:
\begin{matlab}
%}
switch Dom.type(q,:)
%{
\end{matlab}
%: case periodic
The periodic case is evenly spaced within the spatial domain.
Store the size ratio in \verb|patches|.
\begin{matlab}
%}
case 'periodic'
  Q=linspace(Xlim(qq),Xlim(qq+1),nPatch(q)+1);
  DQ=Q(2)-Q(1);
  Q=Q(1:nPatch(q))+diff(Q)/2;
  pEI=patches.EdgyInt;% abbreviation
  if pre2023, dx(q) = ratio(q)*DQ/(nSubP(q)-1-pEI)*(2-pEI);
  else        ratio(q) = dx(q)/DQ*(nSubP(q)-1-pEI)/(2-pEI);  
  end
  patches.ratio=ratio;
%{
\end{matlab}
%: case equispace
The equi-spaced case is also evenly spaced but with the
extreme edges aligned with the spatial domain boundaries,
modified by the offset.
\begin{matlab}
%}
case 'equispace'
  Q=linspace(Xlim(qq)+((nSubP(q)-1)/2-Dom.bcOffset(qq))*dx(q) ...
          ,Xlim(qq+1)-((nSubP(q)-1)/2-Dom.bcOffset(qq+1))*dx(q) ...
          ,nPatch(q));
  DQ=diff(Q(1:2));
  width=(1+patches.EdgyInt)/2*(nSubP(q)-1-patches.EdgyInt)*dx;
  if DQ<width*0.999999
     warning('too many equispace patches (double overlapping)')
     end
%{
\end{matlab}
%: case chebyshev
The Chebyshev case is spaced according to the Chebyshev
distribution in order to reduce macro-interpolation errors,
\(Q_i \propto -\cos(i\pi/N)\),  but with the extreme edges
aligned with the spatial domain boundaries, modified by the
offset, and modified by possible `boundary layers'.
\footnote{ However, maybe overlapping patches near a
boundary should be viewed as some sort of spatially analogue
of the `christmas tree' of projective integration and its
integration to a slow manifold.   Here maybe the overlapping
patches allow for a `christmas tree' approach to the
boundary layers.   Needs to be explored??}  
\begin{matlab}
%}
case 'chebyshev'
  halfWidth=dx(q)*(nSubP(q)-1)/2;
  Q1 = Xlim(1)+halfWidth-Dom.bcOffset(qq)*dx(q);
  Q2 = Xlim(2)-halfWidth+Dom.bcOffset(qq+1)*dx(q);
%  Q = (Q1+Q2)/2-(Q2-Q1)/2*cos(linspace(0,pi,nPatch));
%{
\end{matlab}
Search for total width of `boundary layers' so that in the
interior the patches are non-overlapping Chebyshev.   But
the width for assessing overlap of patches is the following
variable \verb|width|.
\begin{matlab}
%}
  width=(1+patches.EdgyInt)/2*(nSubP(q)-1-patches.EdgyInt)*dx(q);
  for b=0:2:nPatch(q)-2
    DQmin=(Q2-Q1-b*width)/2*( 1-cos(pi/(nPatch(q)-b-1)) );
    if DQmin>width, break, end
  end
  if DQmin<width*0.999999
     warning('too many Chebyshev patches (mid-domain overlap)')
     end
%{
\end{matlab}
Assign the centre-patch coordinates.
\begin{matlab}
%}
  Q =[ Q1+(0:b/2-1)*width ...
       (Q1+Q2)/2-(Q2-Q1-b*width)/2*cos(linspace(0,pi,nPatch(q)-b)) ...
       Q2+(1-b/2:0)*width ];
%{
\end{matlab}

%: case usergiven
The user-given case is entirely up to a user to specify, we just
ensure it has the correct shape of a row??.
\begin{matlab}
%}
case 'usergiven'
  if q==1, Q = reshape(Dom.X,1,[]);
  else     Q = reshape(Dom.Y,1,[]);
  end%if
end%switch Dom.type
%{
\end{matlab}
Assign \(Q\)-coordinates to the correct spatial direction.
At this stage they are all rows.
\begin{matlab}
%}
if q==1, X=Q; end
if q==2, Y=Q; end
end%for q
%{
\end{matlab}







\paragraph{Construct the micro-grids}
Fourth, construct the microscale grid in each patch. Reshape
the grid to be 6D to suit dimensions (micro,Vars,Ens,macro).
\begin{matlab}
%}
nSubP = reshape(nSubP,1,2); % force to be row vector
assert(patches.EdgyInt | all(mod(nSubP,2)==1), ...
    'configPatches2: nSubP must be odd')
i0 = (nSubP(1)+1)/2;
patches.x = reshape( dx(1)*(-i0+1:i0-1)'+X ...
                   ,nSubP(1),1,1,1,nPatch(1),1);
%{
\end{matlab}
Next the \(y\)-direction.
\begin{matlab}
%}
i0 = (nSubP(2)+1)/2;
patches.y = reshape( dx(2)*(-i0+1:i0-1)'+Y ...
                   ,1,nSubP(2),1,1,1,nPatch(2));
%{
\end{matlab}


\paragraph{Pre-compute weights for macro-periodic}
In the case of macro-periodicity, precompute the weightings
to interpolate field values for coupling. \todo{Might sometime
extend to coupling via derivative values.}   
\begin{matlab}
%}
if patches.periodic
  ratio = reshape(ratio,1,2); % force to be row vector
  patches.ratio=ratio; 
  if ordCC>0
    [Cwtsr,Cwtsl] = patchCwts(ratio,ordCC,patches.stag);
    patches.Cwtsr = Cwtsr;  patches.Cwtsl = Cwtsl;
  end%if
end%if patches.periodic
%{
\end{matlab}





\subsection{Set ensemble inter-patch communication}
For \verb|EdgyInt| or centre interpolation respectively, 
\begin{itemize}
\item the right-edge\slash centre realisations
\verb|1:nEnsem| are to interpolate to left-edge~\verb|le|,
and 
\item the left-edge\slash centre realisations
\verb|1:nEnsem| are to interpolate to~\verb|re|.
\end{itemize}
\verb|re| and \verb|li| are `transposes' of each other as
\verb|re(li)=le(ri)| are both \verb|1:nEnsem|. Similarly for
bottom-edge\slash centre interpolation to top-edge
via~\verb|to|, and top-edge\slash centre interpolation to
bottom-edge via~\verb|bo|.

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
many times in order to subsequently sub-sample: \verb|nSubP|
times should be enough. If \verb|cs| is more then 3D, then
the higher-dimensions are reshaped into the 3rd dimension.
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
$m_xm_y$~phase-shifts of the coefficients is constructed
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
lack of scale separation.  \todo{Maybe need to justify this
and the arbitrary threshold more carefully??}
\begin{matlab}
%}
if prod(ratio)*patches.nEnsem>0.9, warning( ...
'Probably poor scale separation in ensemble of coupled phase-shifts')
scaleSeparationParameter = ratio*patches.nEnsem
end
%{
\end{matlab}
End the two if-statements.
\begin{matlab}
%}
  end%if-else nEnsem>1  
end%if not-empty(cs)
%{
\end{matlab}



\paragraph{If parallel code} then first assume this is not
within an \verb|spmd|-environment, and so we invoke
\verb|spmd...end| (which starts a parallel pool if not
already started).  At this point, the global \verb|patches|
is copied for each worker processor and so it becomes
\emph{composite} when we distribute any one of the fields.
Hereafter, {\em all fields in the global variable
\verb|patches| must only be referenced within an
\verb|spmd|-environment.}%
\footnote{If subsequently outside spmd, then one must use
functions like \texttt{getfield(patches\{1\},'a')}.}
\begin{matlab}
%}
if patches.parallel
%  theparpool=gcp()
  spmd
%{
\end{matlab}
Second, decide which dimension is to be sliced among
parallel workers (for the moment, do not consider slicing
the ensemble). Choose the direction of most patches, biased
towards the last.
\begin{matlab}
%}
  [~,pari]=max(nPatch+0.01*(1:2));
  patches.codist=codistributor1d(4+pari);
%{
\end{matlab}
\verb|patches.codist.Dimension| is the index that is split
among workers.  Then distribute the appropriate coordinate
direction among the workers: the function must be invoked
inside an \verb|spmd|-group in order for this to work---so
we do not need \verb|parallel| in argument list.
\begin{matlab}
%}
  switch pari
    case 1, patches.x=codistributed(patches.x,patches.codist);
    case 2, patches.y=codistributed(patches.y,patches.codist);
  otherwise
    error('should never have bad index for parallel distribution')
  end%switch
  end%spmd
%{
\end{matlab}

If not parallel, then clean out \verb|patches.codist| if it exists.
May not need, but safer.
\begin{matlab}
%}
else% not parallel
  if isfield(patches,'codist'), rmfield(patches,'codist'); end
end%if-parallel
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
