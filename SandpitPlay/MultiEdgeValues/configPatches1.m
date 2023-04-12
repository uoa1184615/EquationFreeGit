% configPatches1() creates a data struct of the design of
% 1D patches for later use by the patch functions such as
% patchSys1(). AJR, Nov 2017 -- 23 Mar 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches1()}: configure spatial
patches in 1D}
\label{sec:configPatches1}
\localtableofcontents



Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative\slash step function
\verb|patchSys1()|. \cref{sec:configPatches1eg} lists an
example of its use.
\begin{matlab}
%}
function patches = configPatches1(fun,Xlim,Dom ...
    ,nPatch,ordCC,dx,nSubP,varargin)
version = '2023-03-23';
%{
\end{matlab}



\paragraph{Input}
If invoked with no input arguments, then executes an example
of simulating Burgers' \pde---see \cref{sec:configPatches1eg}
for the example code.
\begin{itemize}

\item \verb|fun| is the name of the user function,
\verb|fun(t,u,patches)| or \verb|fun(t,u)| or
\verb|fun(t,u,patches,...)|, that computes time derivatives
(or time-steps) of quantities on the 1D micro-grid within
all the 1D~patches.

\item \verb|Xlim| give the macro-space spatial domain of the
computation, namely the interval $[ \verb|Xlim(1)|,
\verb|Xlim(2)|]$.

\item \verb|Dom| sets the type of macroscale conditions for
the patches, and reflects the type of microscale boundary
conditions of the problem.   If \verb|Dom| is \verb|NaN| or
\verb|[]|, then the field~\verb|u| is macro-periodic in the
1D spatial domain, and resolved on equi-spaced patches. If
\verb|Dom| is a character string, then that specifies the
\verb|.type| of the following structure, with
\verb|.bcOffset| set to the default zero. Otherwise
\verb|Dom| is a structure with the following components.
\begin{itemize}

\item \verb|.type|, string, of either \verb|'periodic'| (the
default), \verb|'equispace'|, \verb|'chebyshev'|,
\verb|'usergiven'|.  For all cases except \verb|'periodic'|,
users \emph{must} code into \verb|fun| the micro-grid
boundary conditions that apply at the left(right) edge of
the leftmost(rightmost) patches.

\item \verb|.bcOffset|, optional one or two element array,
in the cases of \verb|'equispace'| or \verb|'chebyshev'|
the patches are placed so the left\slash right macroscale
boundaries are aligned to the left\slash right edges of the
corresponding extreme patches, but offset by \verb|bcOffset|
of the sub-patch micro-grid spacing.  For example, use
\verb|bcOffset=0| when applying Dirichlet boundary values on
the extreme edge micro-grid points, whereas use
\verb|bcOffset=0.5| when applying Neumann boundary conditions
halfway between the extreme edge micro-grid points.

\item \verb|.X|, optional array, in the case~\verb|'usergiven'|
it specifies the locations of the centres of the
\verb|nPatch| patches---the user is responsible it makes
sense.
\end{itemize}


\item \verb|nPatch| is the number of equi-spaced spatial
patches.

\item \verb|ordCC|, must be~$\geq -1$, is the `order' of
interpolation across empty space of the macroscale patch
values to the edge of the patches for inter-patch coupling:
where \verb|ordCC| of~$0$ or~$-1$ gives spectral
interpolation; and \verb|ordCC| being odd specifies
staggered spatial grids.

\item \verb|dx| (real) is usually the sub-patch micro-grid
spacing in~\(x\).

However, if \verb|Dom| is~\verb|NaN| (as for pre-2023), then
\verb|dx| actually is \verb|ratio|, namely the ratio of
(depending upon \verb|EdgyInt|) either the half-width or
full-width of a patch to the equi-spacing of the patch
mid-points---adjusted a little when $\verb|nEdge|>1$. So
either $\verb|ratio|=\tfrac12$ means the patches abut and
$\verb|ratio|=1$ is overlapping patches as in holistic
discretisation, or $\verb|ratio|=1$ means the patches abut. 
Small~\verb|ratio| should greatly reduce computational time.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch. If not using \verb|EdgyInt|,
then $\verb|nSubP/nEdge|$ must be odd integer so that there 
is/are centre-patch lattice point(s).  So for the defaults 
of $\verb|nEdge|=1$ and not \verb|EdgyInt|, then 
\verb|nSubP| must be odd.

\item \verb|'nEdge'|, \emph{optional}, default=1, the number
of edge values set by interpolation at the edge regions of
each patch.  The default is one (suitable for microscale
lattices with only nearest neighbour interactions).

\item \verb|EdgyInt|, true/false, \emph{optional},
default=false.  If true, then interpolate to left\slash
right edge-values from right\slash left next-to-edge values.
If false or omitted, then interpolate from centre-patch
values.

\item \verb|nEnsem|, \emph{optional-experimental},
default one, but if more, then an ensemble over this
number of realisations.

\item \verb|hetCoeffs|, \emph{optional}, default empty.
Supply a 1D or 2D array of microscale heterogeneous
coefficients to be used by the given microscale \verb|fun|
in each patch. Say the given array~\verb|cs| is of size
$m_x\times n_c$, where $n_c$~is the number of different sets
of coefficients. The coefficients are to be the same for
each and every patch; however, macroscale variations are
catered for by the $n_c$~coefficients being $n_c$~parameters
in some macroscale formula.
\begin{itemize}
\item If $\verb|nEnsem|=1$, then the array of coefficients
is just tiled across the patch size to fill up each patch,
starting from the first point in each patch.  Best accuracy 
usually obtained when the periodicity of the coefficients 
is a factor of \verb|nSubP-2*nEdge| for \verb|EdgyInt|, or 
a factor of \verb|(nSubP-nEdge)/2| for not \verb|EdgyInt|.

\item If $\verb|nEnsem|>1$ (value immaterial), then reset
$\verb|nEnsem|:=m_x$ and construct an ensemble of all
$m_x$~phase-shifts of the coefficients. In this scenario,
the inter-patch coupling couples different members in the
ensemble.  When \verb|EdgyInt| is true, and when the
coefficients are diffusivities\slash elasticities, then this
coupling cunningly preserves symmetry.
\end{itemize}

\item \verb|nCore|,  \emph{optional-experimental}, default
one, but if more, and only for non-EdgyInt, then
interpolates from an average over the core of a patch, a
core of size ??.  Then edge values are set according to
interpolation of the averages?? or so that average at edges
is the interpolant??

\item \verb|'parallel'|, true/false, \emph{optional},
default=false. If false, then all patch computations are on
the user's main \textsc{cpu}---although a user may well
separately invoke, say, a \textsc{gpu} to accelerate
sub-patch computations. 

If true, and it requires that you have \Matlab's Parallel
Computing Toolbox, then it will distribute the patches over
multiple \textsc{cpu}s\slash cores. In \Matlab, only one
array dimension can be split in the distribution, so it
chooses the one space dimension~$x$.  A user may
correspondingly distribute arrays with property
\verb|patches.codist|, or simply use formulas invoking the
preset distributed arrays \verb|patches.x|. If a user has
not yet established a parallel pool, then a `local' pool is
started.

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
patches.version = version;
%{
\end{matlab}
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,patches)| or \verb|fun(t,u)| or
\verb|fun(t,u,patches,...)|, that computes the time
derivatives (or steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.periodic|: either true, for interpolation on
the macro-periodic domain; or false, for general
interpolation by divided differences over non-periodic
domain or unevenly distributed patches.

\item \verb|.stag| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling.

\item \verb|.Cwtsr| and \verb|.Cwtsl|, only for
macro-periodic conditions, are the $\verb|ordCC|$-vector of
weights for the inter-patch interpolation onto the right and
left edges (respectively) with patch:macroscale ratio as
specified or as derived from~\verb|dx|.

\item \verb|.x| (4D) is $\verb|nSubP| \times1 \times1
\times \verb|nPatch|$ array of the regular spatial
locations~$x_{iI}$ of the $i$th~microscale grid point in
the $I$th~patch.  

\item \verb|.ratio|, only for macro-periodic conditions, is
the size ratio of every patch.

\item \verb|.nEdge| is, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch.

\item \verb|.le|, \verb|.ri| determine inter-patch coupling
of members in an ensemble. Each a column vector of
length~\verb|nEnsem|.

\item \verb|.cs| either
\begin{itemize}
\item \verb|[]| 0D, or 
\item if $\verb|nEnsem|=1$, $(\verb|nSubP(1)|-1)\times
n_c$ 2D array of microscale heterogeneous coefficients, or
\item if $\verb|nEnsem|>1$, $(\verb|nSubP(1)|-1) \times
n_c\times m_x$ 3D array of $m_x$~ensemble of phase-shifts
of the microscale
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
\label{sec:configPatches1eg}
\begin{matlab}
%}
if nargin==0
disp('With no arguments, simulate example of Burgers PDE')
%{
\end{matlab}
The code here shows one way to get started: a user's script
may have the following three steps (``\into'' denotes
function recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches1 
\item ode15s integrator \into patchSys1 \into user's PDE
\item process results
\end{enumerate}

Establish global patch data struct to point to and interface
with a function coding Burgers' \pde: to be solved on
$2\pi$-periodic domain, with eight patches, spectral
interpolation couples the patches, with micro-grid
spacing~$0.06$, and with seven microscale points forming
each patch.
\begin{matlab}
%}
global patches
patches = configPatches1(@BurgersPDE, [0 2*pi], ...
    'periodic', 8, 0, 0.06, 7);
%{
\end{matlab}
Set some initial condition, with some microscale randomness.
\begin{matlab}
%}
u0=0.3*(1+sin(patches.x))+0.1*randn(size(patches.x));
%{
\end{matlab}
Simulate in time using a standard stiff integrator and the
interface function \verb|patchSys1()|
(\cref{sec:patchSys1}).
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
[ts,us] = ode15s( @patchSys1,[0 0.5],u0(:));
else % octave version
[ts,us] = odeOcts(@patchSys1,[0 0.5],u0(:));
end
%{
\end{matlab}
Plot the simulation using only the microscale values
interior to the patches: either set $x$-edges to \verb|nan|
to leave the gaps; or use \verb|patchEdgyInt1| to
re-interpolate correct patch edge values and thereby join
the patches.  \cref{fig:config1Burgers} illustrates an
example simulation in time generated by the patch scheme
applied to Burgers'~\pde.
\begin{matlab}
%}
figure(1),clf
if 1, patches.x([1 end],:,:,:)=nan;  us=us.';
else us=reshape(patchEdgyInt1(us.'),[],length(ts));  
end
mesh(ts,patches.x(:),us)
view(60,40), colormap(0.7*hsv)
title('Burgers PDE: patches in space, continuous time')
xlabel('time t'), ylabel('space x'), zlabel('u(x,t)')
%{
\end{matlab}

\begin{figure}
\centering \caption{\label{fig:config1Burgers}field
$u(x,t)$ of the patch scheme applied to Burgers'~\pde.}
\includegraphics[scale=0.85]{configPatches1}
\end{figure}
Upon finishing execution of the example, optionally save 
the graph to be shown in \cref{fig:config1Burgers}, then 
exit this function.
\begin{matlab}
%}
ifOurCf2eps(mfilename)
return
end%if nargin==0
%{
\end{matlab}

\IfFileExists{../Patch/BurgersPDE.m}{\input{../Patch/BurgersPDE.m}}{}
\IfFileExists{../Patch/odeOcts.m}{\input{../Patch/odeOcts.m}}{}




\begin{devMan}

\subsection{Parse input arguments and defaults}
\begin{matlab}
%}
p = inputParser;
fnValidation = @(f) isa(f, 'function_handle'); %test for fn name
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
addParameter(p,'nCore',1,@isnumeric);
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
patches.nCore = p.Results.nCore;
%{
\end{matlab}

Check parameters.
\begin{matlab}
%}
assert(Xlim(1)<Xlim(2) ...
      ,'two entries of Xlim must be ordered increasing')
assert((mod(ordCC,2)==0)|(patches.nEdge==1) ...
      ,'Cannot yet have nEdge>1 and staggered patch grids')
assert(3*patches.nEdge<=nSubP ...
      ,'too many edge values requested')
assert(rem(nSubP,patches.nEdge)==0 ...
      ,'nSubP must be integer multiple of nEdge')
if ~patches.EdgyInt, assert(rem(nSubP/patches.nEdge,2)==1 ...
      ,'for non-edgyInt, nSubP/nEdge must be odd integer')
      end
if (patches.nEnsem>1)&(patches.nEdge>1)
      warning('not yet tested when both nEnsem and nEdge non-one')
      end
if patches.nCore>1
    warning('nCore>1 not yet tested in this version')
    end
%{
\end{matlab}


For compatibility with pre-2023 functions, if parameter
\verb|Dom| is \verb|Nan|, then  we set the \verb|ratio| to
be the value of the so-called \verb|dx| parameter.
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
string, and then get corresponding defaults for others
fields.
\begin{matlab}
%}
if ischar(Dom), Dom=struct('type',Dom); end
%{
\end{matlab}
Check what is and is not specified, and provide default of
Dirichlet boundaries if no \verb|bcOffset| specified when
needed.
\begin{matlab}
%}
patches.periodic=false;
switch Dom.type
case 'periodic'
    patches.periodic=true;
    if isfield(Dom,'bcOffset')
    warning('bcOffset not available for Dom.type = periodic'), end
    if isfield(Dom,'X')
    warning('X not available for Dom.type = periodic'), end
case {'equispace','chebyshev'}
    if ~isfield(Dom,'bcOffset'), Dom.bcOffset=[0;0]; end
    if length(Dom.bcOffset)==1
        Dom.bcOffset=repmat(Dom.bcOffset,2,1); end
    if isfield(Dom,'X')
    warning('X not available for Dom.type = equispace or chebyshev')
    end
case 'usergiven'
    if isfield(Dom,'bcOffset')
    warning('bcOffset not available for usergiven Dom.type'), end
    assert(isfield(Dom,'X'),'X required for Dom.type = usergiven')
otherwise 
    error([Dom.type ' is unknown Dom.type'])
end%switch Dom.type
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
coupling is \verb|ordCC| of~$0$ and~$-1$.
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


Third, set the centre of the patches in the macroscale grid
of patches, depending upon \verb|Dom.type|.
\begin{matlab}
%}
switch Dom.type
%{
\end{matlab}
%: case periodic
The periodic case is evenly spaced within the spatial domain.
Store the size ratio in \verb|patches|.
\begin{matlab}
%}
case 'periodic'
  X=linspace(Xlim(1),Xlim(2),nPatch+1);
  DX=X(2)-X(1);
  X=X(1:nPatch)+diff(X)/2;
  pEI=patches.EdgyInt;% abbreviation
  pnE=patches.nEdge;  % abbreviation
  if pre2023, dx = ratio*DX/(nSubP-pnE*(1+pEI))*(2-pEI);
  else        ratio = dx/DX*(nSubP-pnE*(1+pEI))/(2-pEI);  end
  patches.ratio=ratio;
%{
\end{matlab}
In the case of macro-periodicity, precompute the weightings
to interpolate field values for coupling. 
\todo{Might sometime extend to coupling via derivative values.}   
\begin{matlab}
%}
  if ordCC>0
    [Cwtsr,Cwtsl] = patchCwts(ratio,ordCC,patches.stag);
    patches.Cwtsr = Cwtsr;  patches.Cwtsl = Cwtsl;
  end
%{
\end{matlab}
%: case equispace
The equi-spaced case is also evenly spaced but with the
extreme edges aligned with the spatial domain boundaries,
modified by the offset.
%\todo{This warning needs refinement for multi-edges??}
\begin{matlab}
%}
case 'equispace'
  X=linspace(Xlim(1)+((nSubP-1)/2-Dom.bcOffset(1))*dx ...
            ,Xlim(2)-((nSubP-1)/2-Dom.bcOffset(2))*dx ,nPatch);
  DX=diff(X(1:2));
  width=(1+patches.EdgyInt)/2*(nSubP-1-patches.EdgyInt)*dx;
  if DX<width*0.999999
     warning('too many equispace patches (double overlapping)')
     end
%{
\end{matlab}
%: case chebyshev
The Chebyshev case is spaced according to the Chebyshev
distribution in order to reduce macro-interpolation errors,
\(X_i \propto -\cos(i\pi/N)\),  but with the extreme edges
aligned with the spatial domain boundaries, modified by the
offset, and modified by possible `boundary layers'.
\footnote{ However, maybe overlapping patches near a
boundary should be viewed as some sort of spatial analogue
of the `christmas tree' of projective integration and its
projection to a slow manifold.   Here maybe the overlapping
patches allow for a `christmas tree' approach to the
boundary layers.   Needs to be explored??}
\begin{matlab}
%}
case 'chebyshev'
  halfWidth=dx*(nSubP-1)/2;
  X1 = Xlim(1)+halfWidth-Dom.bcOffset(1)*dx;
  X2 = Xlim(2)-halfWidth+Dom.bcOffset(2)*dx;
%  X = (X1+X2)/2-(X2-X1)/2*cos(linspace(0,pi,nPatch));
%{
\end{matlab}
Search for total width of `boundary layers' so that in the
interior the patches are non-overlapping Chebyshev.   But
the width for assessing overlap of patches is the following
variable \verb|width|.  We need to find~\verb|b|, the number of patches `glued' together at the boundaries. 
\begin{matlab}
%}
  pEI=patches.EdgyInt;% abbreviation
  pnE=patches.nEdge;  % abbreviation
  width=(1+pEI)/2*(nSubP-pnE-pEI*pnE)*dx;
  for b=0:2:nPatch-2
    DXmin=(X2-X1-b*width)/2*( 1-cos(pi/(nPatch-b-1)) );
    if DXmin>width, break, end
  end%for
  if DXmin<width*0.999999
     warning('too many Chebyshev patches (mid-domain overlap)')
     end
%{
\end{matlab}
Assign the centre-patch coordinates.
\begin{matlab}
%}
  X = [ X1+(0:b/2-1)*width ...
        (X1+X2)/2-(X2-X1-b*width)/2*cos(linspace(0,pi,nPatch-b)) ...
        X2+(1-b/2:0)*width ];
%{
\end{matlab}

%: case usergiven
The user-given case is entirely up to a user to specify, we
just force it to have the correct shape of a row.
\begin{matlab}
%}
case 'usergiven'
  X = reshape(Dom.X,1,[]);
end%switch Dom.type
%{
\end{matlab}


Fourth, construct the microscale grid in each patch, centred
about the given mid-points~\verb|X|. Reshape the grid to be
4D to suit dimensions (micro,Vars,Ens,macro).
\begin{matlab}
%}
xs = dx*( (1:nSubP)-mean(1:nSubP) );
patches.x = reshape( xs'+X ,nSubP,1,1,nPatch);
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
\verb|re(li)=le(ri)| are both \verb|1:nEnsem|.
Alternatively, one may use the statement
\begin{verbatim}
  c=hankel(c(1:nSubP-1),c([nSubP 1:nSubP-2]));
\end{verbatim}
to \emph{correspondingly} generates all phase shifted copies
of microscale heterogeneity (see \verb|homoDiffEdgy1| of
\cref{sec:homoDiffEdgy1}).

The default is nothing shifty.  This setting reduces the
number of if-statements in function \verb|patchEdgeInt1()|.
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
many times in order to subsequently sub-sample: \verb|nSubP|
times should be enough. If \verb|cs| is more then 2D, then
the higher-dimensions are reshaped into the 2nd dimension.
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
$m_x$~phase-shifts of the coefficients is constructed from
the over-supply.  Here code phase-shifts over the
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
Further, set a cunning left\slash right realisation of
inter-patch coupling.  The aim is to preserve symmetry in
the system when also invoking \verb|EdgyInt|.  What this
coupling does without \verb|EdgyInt| is unknown.  Use
auto-replication.
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
Second, choose to slice parallel workers in the spatial
direction.
\begin{matlab}
%}
  pari = 1;
  patches.codist=codistributor1d(3+pari);
%{
\end{matlab}
\verb|patches.codist.Dimension| is the index that is split
among workers.  Then distribute the coordinate direction
among the workers: the function must be invoked inside an
\verb|spmd|-group in order for this to work---so we do not
need \verb|parallel| in argument list.
\begin{matlab}
%}
  switch pari
    case 1, patches.x=codistributed(patches.x,patches.codist);
  otherwise
    error('should never have bad index for parallel distribution')
  end%switch
  end%spmd
%{
\end{matlab}

If not parallel, then clean out \verb|patches.codist| if it
exists. May not need, but safer.
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


