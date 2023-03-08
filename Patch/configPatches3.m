% configPatches3() creates a data struct of the design of 3D
% patches for later use by the patch functions such as
% patchSys3().  AJR, Aug 2020 -- 2 Feb 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches3()}: configures spatial
patches in 3D}
\label{sec:configPatches3}
\localtableofcontents



Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative\slash step function
\verb|patchSys3()|, and possibly other patch functions.
\cref{sec:configPatches3eg,sec:homoDiffEdgy3} list examples
of its use.
\begin{matlab}
%}
function patches = configPatches3(fun,Xlim,Dom ...
    ,nPatch,ordCC,dx,nSubP,varargin)
%{
\end{matlab}



\paragraph{Input}
If invoked with no input arguments, then executes an example
of simulating a heterogeneous wave \pde---see
\cref{sec:configPatches3eg} for an example code.
\begin{itemize}

\item \verb|fun| is the name of the user function,
\verb|fun(t,u,patches)| or \verb|fun(t,u)| or
\verb|fun(t,u,patches,...)|, that computes time-derivatives
(or time-steps) of quantities on the 3D micro-grid within
all the 3D~patches.

\item \verb|Xlim| array/vector giving the rectangular-cuboid
macro-space domain of the computation: namely
$[\verb|Xlim(1)|, \verb|Xlim(2)|] \times [\verb|Xlim(3)|,
\verb|Xlim(4)| \times [\verb|Xlim(5)|, \verb|Xlim(6)|]$. If
\verb|Xlim| has two elements, then the domain is the cubic
domain of the same interval in all three directions.

\item \verb|Dom| sets the type of macroscale conditions for
the patches, and reflects the type of microscale boundary
conditions of the problem.   If \verb|Dom| is \verb|NaN| or
\verb|[]|, then the field~\verb|u| is triply macro-periodic
in the 3D spatial domain, and resolved on equi-spaced
patches. If \verb|Dom| is a character string, then that
specifies the \verb|.type| of the following structure, with
\verb|.bcOffset| set to the default zero.  Otherwise
\verb|Dom| is a structure with the following components.
\begin{itemize}

\item \verb|.type|, string, of either \verb|'periodic'| (the
default), \verb|'equispace'|, \verb|'chebyshev'|,
\verb|'usergiven'|.  For all cases except \verb|'periodic'|,
users \emph{must} code into \verb|fun| the micro-grid
boundary conditions that apply at the left\slash right\slash
bottom\slash top\slash back\slash front faces of the
leftmost\slash rightmost\slash bottommost\slash
topmost\slash backmost\slash frontmost patches,
respectively.

\item \verb|.bcOffset|, optional one, three or six element
vector/array, in the cases of \verb|'equispace'| or
\verb|'chebyshev'| the patches are placed so the left\slash
right macroscale boundaries are aligned to the left\slash
right faces of the corresponding extreme patches, but offset
by \verb|bcOffset| of the sub-patch micro-grid spacing.  For
example, use \verb|bcOffset=0| when the micro-code applies
Dirichlet boundary values on the extreme face micro-grid
points, whereas use \verb|bcOffset=0.5| when the microcode
applies Neumann boundary conditions halfway between the
extreme face micro-grid points.  Similarly for the top,
bottom, back, and front faces.

If \verb|.bcOffset| is a scalar, then apply the same offset
to all boundaries. If three elements, then apply the first
offset to both \(x\)-boundaries, the second offset to both
\(y\)-boundaries, and the third offset to both
\(z\)-boundaries. If six elements, then apply the first two
offsets to the respective \(x\)-boundaries, the middle two
offsets to the respective \(y\)-boundaries, and the last two
offsets to the respective \(z\)-boundaries.

\item \verb|.X|, optional vector/array with \verb|nPatch(1)|
elements, in the case \verb|'usergiven'| it specifies the
\(x\)-locations of the centres of the patches---the user is
responsible the locations makes sense.

\item \verb|.Y|, optional vector/array with \verb|nPatch(2)|
elements, in the case \verb|'usergiven'| it specifies the
\(y\)-locations of the centres of the patches---the user is
responsible the locations makes sense.

\item \verb|.Z|, optional vector/array with \verb|nPatch(3)|
elements, in the case \verb|'usergiven'| it specifies the
\(z\)-locations of the centres of the patches---the user is
responsible the locations makes sense.
\end{itemize}


\item \verb|nPatch| sets the number of equi-spaced spatial
patches: if scalar, then use the same number of patches in
all three directions, otherwise \verb|nPatch(1:3)| gives the
number~($\geq1$) of patches in each direction.

\item \verb|ordCC| is the `order' of interpolation for
inter-patch coupling across empty space of the macroscale
patch values to the face-values of the patches: currently
must be~$0,2,4,\ldots$; where $0$~gives spectral
interpolation.

\item \verb|dx| (real---scalar or three elements) is usually
the sub-patch micro-grid spacing in~\(x\), \(y\) and~\(z\). 
If scalar, then use the same \verb|dx| in all three
directions, otherwise \verb|dx(1:3)| gives the spacing in
each of the three directions.

However, if \verb|Dom| is~\verb|NaN| (as for pre-2023), then
\verb|dx| actually is \verb|ratio| (scalar or three elements),
namely the ratio of (depending upon \verb|EdgyInt|) either
the half-width or full-width of a patch to the equi-spacing
of the patch mid-points.  So either $\verb|ratio|=\tfrac12$
means the patches abut and $\verb|ratio|=1$ is overlapping
patches as in holistic discretisation, or $\verb|ratio|=1$
means the patches abut.  Small~\verb|ratio| should greatly
reduce computational time.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch: if scalar, then use the same
number in all three directions, otherwise \verb|nSubP(1:3)|
gives the number in each direction. If not using
\verb|EdgyInt|, then must be odd so that there is/are
centre-patch micro-grid point\slash planes in each patch.

\item \verb|'nEdge'| (not yet implemented), \emph{optional},
default=1, for each patch, the number of face values set by
interpolation at the face regions of each patch.  The
default is one (suitable for microscale lattices with only
nearest neighbour interactions).

\item \verb|'EdgyInt'|, true/false, \emph{optional},
default=false.  If true, then interpolate to left\slash
right\slash top\slash bottom\slash front\slash back
face-values from right\slash left\slash bottom\slash
top\slash back\slash front next-to-face values.  If false or
omitted, then interpolate from centre-patch planes.  

\item \verb|'nEnsem'|,  \emph{optional-experimental},
default one, but if more, then an ensemble over this number
of realisations.

\item \verb|'hetCoeffs'|, \emph{optional}, default empty.
Supply a 3D or 4D array of microscale heterogeneous
coefficients to be used by the given microscale \verb|fun|
in each patch. Say the given array~\verb|cs| is of size
$m_x\times m_y\times m_z\times n_c$, where $n_c$~is the
number of different arrays of coefficients.  For example, in
heterogeneous diffusion, $n_c=3$ for the diffusivities in
the \emph{three} different spatial directions (or $n_c=6$
for the diffusivity tensor).  The coefficients are to be the
same for each and every patch. However, macroscale
variations are catered for by the $n_c$~coefficients being
$n_c$~parameters in some macroscale formula.
\begin{itemize}
\item If $\verb|nEnsem|=1$, then the array of coefficients
is just tiled across the patch size to fill up each patch,
starting from the $(1,1,1)$-point in each patch.

\item If $\verb|nEnsem|>1$ (value immaterial), then reset
$\verb|nEnsem|:=m_x\cdot m_y\cdot m_z$ and construct an
ensemble of all $m_x\cdot m_y\cdot m_z$ phase-shifts of
the coefficients.  In this scenario, the inter-patch
coupling couples different members in the ensemble.  When
\verb|EdgyInt| is true, and when the coefficients are
diffusivities\slash elasticities in $x,y,z$-directions,
respectively, then this coupling cunningly preserves
symmetry.

\end{itemize}

\item \verb|'parallel'|, true/false, \emph{optional},
default=false. If false, then all patch computations are on
the user's main \textsc{cpu}---although a user may well
separately invoke, say, a \textsc{gpu} to accelerate
sub-patch computations. 

If true, and it requires that you have \Matlab's Parallel
Computing Toolbox, then it will distribute the patches over
multiple \textsc{cpu}s\slash cores. In \Matlab, only one
array dimension can be split in the distribution, so it
chooses the one space dimension~$x,y,z$ corresponding to
the highest~\verb|nPatch| (if a tie, then chooses the
rightmost of~$x,y,z$).  A user may correspondingly
distribute arrays with property \verb|patches.codist|, or
simply use formulas invoking the preset distributed arrays
\verb|patches.x|, \verb|patches.y|, and \verb|patches.z|. If
a user has not yet established a parallel pool, then a
`local' pool is started.

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
\verb|fun(t,u,patches)| or \verb|fun(t,u)| or
\verb|fun(t,u,patches,...)| that computes the time
derivatives (or steps) on the patchy lattice. 

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

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the
$\verb|ordCC|\times 3$-array of weights for the inter-patch
interpolation onto the right\slash top\slash front and
left\slash bottom\slash back faces (respectively) with
patch:macroscale ratio as specified or as derived
from~\verb|dx|.

\item \verb|.x| (8D) is $\verb|nSubP(1)| \times1 \times1
\times1 \times1 \times \verb|nPatch(1)| \times1 \times1$
array of the regular spatial locations~$x_{iI}$ of the
microscale grid points in every patch.  

\item \verb|.y| (8D) is $1 \times \verb|nSubP(2)| \times1
\times1 \times1 \times1 \times \verb|nPatch(2)| \times1$
array of the regular spatial locations~$y_{jJ}$ of the
microscale grid points in every patch.  

\item \verb|.z| (8D) is $1 \times1 \times \verb|nSubP(3)|
\times1 \times1 \times1 \times1 \times \verb|nPatch(3)|$
array of the regular spatial locations~$z_{kK}$ of the
microscale grid points in every patch.  

\item \verb|.ratio| $1\times 3$, only for macro-periodic
conditions, are the size ratios of every patch.

\item \verb|.nEdge| is, for each patch, the number of face
values set by interpolation at the face regions of each
patch.

\item \verb|.le|, \verb|.ri|, \verb|.bo|, \verb|.to|,
\verb|.ba|, \verb|.fr| determine inter-patch coupling of
members in an ensemble. Each a column vector of
length~\verb|nEnsem|.

\item \verb|.cs| either
\begin{itemize}
\item \verb|[]| 0D, or 
\item if $\verb|nEnsem|=1$, $(\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times (\verb|nSubP(3)|-1)\times n_c$ 4D
array of microscale heterogeneous coefficients, or
\item if $\verb|nEnsem|>1$, $(\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times (\verb|nSubP(3)|-1)\times
n_c\times m_xm_ym_z$ 5D array of $m_xm_ym_z$~ensemble of
phase-shifts of the microscale heterogeneous coefficients.
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
\label{sec:configPatches3eg}
\begin{matlab}
%}
if nargin==0
disp('With no arguments, simulate example of heterogeneous wave')
%{
\end{matlab}
The code here shows one way to get started: a user's script
may have the following three steps (``\into'' denotes function
recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches3 
\item ode23 integrator \into patchSys3 \into user's PDE
\item process results
\end{enumerate}

Set random heterogeneous
coefficients of period two in each of the three
directions. Crudely normalise by the harmonic mean so the
macro-wave time scale is roughly one. 
\begin{matlab}
%}
mPeriod = [2 2 2];
cHetr = exp(0.9*randn([mPeriod 3]));
cHetr = cHetr*mean(1./cHetr(:)) 
%{
\end{matlab}

Establish global patch data struct to interface with a
function coding a nonlinear `diffusion' \pde: to be solved
on $[-\pi,\pi]^3$-periodic domain, with $5^3$~patches,
spectral interpolation~($0$) couples the patches, each patch
with micro-grid spacing~$0.22$ (relatively large for
visualisation), and with $4^3$~points forming each patch.  
\begin{matlab}
%}
global patches
patches = configPatches3(@heteroWave3,[-pi pi], 'periodic' ...
    , 5, 0, 0.22, mPeriod+2 ,'EdgyInt',true ...
    ,'hetCoeffs',cHetr);
%{
\end{matlab}
Set a wave initial state using auto-replication of the
spatial grid, and as \cref{fig:configPatches3ic} shows. This
wave propagates diagonally across space. Concatenate the two
\(u,v\)-fields to be the two components of the fourth
dimension.
\begin{matlab}
%}
u0 = 0.5+0.5*sin(patches.x+patches.y+patches.z);
v0 =    -0.5*cos(patches.x+patches.y+patches.z)*sqrt(3);
uv0 = cat(4,u0,v0);
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:configPatches3ic}initial
field~$u(x,y,z,t)$ at time $t=0$ of the patch scheme
applied to a heterogeneous wave~\pde:
\cref{fig:configPatches3fin} plots the computed field at
time $t=6$.}
\includegraphics[scale=0.9]{configPatches3ic}
\end{figure}
Integrate in time to $t=6$ using standard functions. In
Matlab \verb|ode15s| would be natural as the patch scheme is
naturally stiff, but \verb|ode23| is much quicker
\cite[Fig.~4]{Maclean2020a}.
\begin{matlab}
%}
disp('Simulate heterogeneous wave u_tt=div[C*grad(u)]')
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSys3,linspace(0,6),uv0(:));
else %disp('octave version is very slow for me')
    lsode_options('absolute tolerance',1e-4);
    lsode_options('relative tolerance',1e-4);
    [ts,us] = odeOcts(@patchSys3,[0 1 2],uv0(:));
end
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:configPatches3fin}.  Use \verb|patchEdgeInt3| to
obtain patch-face values in order to most easily reconstruct
the array data structure.

Replicate $x$, $y$, and~$z$ arrays to get individual
spatial coordinates of every data point.  Then, optionally,
set faces to~\verb|nan| so the plot just shows
patch-interior data. 
\begin{matlab}
%}
figure(1), clf, colormap(0.8*jet)
xs = patches.x+0*patches.y+0*patches.z;
ys = patches.y+0*patches.x+0*patches.z;
zs = patches.z+0*patches.y+0*patches.x;
if 1, xs([1 end],:,:,:)=nan;  
      xs(:,[1 end],:,:)=nan;  
      xs(:,:,[1 end],:)=nan;  
end;%option
j=find(~isnan(xs));
%{
\end{matlab}
In the scatter plot, these functions \verb|pix()| and
\verb|col()| map the $u$-data values to the size of the
dots and to the colour of the dots, respectively.
\begin{matlab}
%}
pix = @(u) 15*abs(u)+7;
col = @(u) sign(u).*abs(u);
%{
\end{matlab}
Loop to plot at each and every time step.
\begin{matlab}
%}
for i = 1:length(ts)
  uv = patchEdgeInt3(us(i,:));
  u = uv(:,:,:,1,:);
  for p=1:2
    subplot(1,2,p)
    if (i==1)| exist('OCTAVE_VERSION','builtin')
      scat(p) = scatter3(xs(j),ys(j),zs(j),'filled');
      axis equal, caxis(col([0 1])), view(45-5*p,25)
      xlabel('x'), ylabel('y'), zlabel('z')
      title('view stereo pair cross-eyed')
    end % in matlab just update values
    set(scat(p),'CData',col(u(j)) ...
       ,'SizeData',pix((8+xs(j)-ys(j)+zs(j))/6+0*u(j)));
    legend(['time = ' num2str(ts(i),'%4.2f')],'Location','north')
  end
%{
\end{matlab}
Optionally save the initial condition to graphic file for
\cref{fig:configPatches2ic}, and optionally save the last
plot.
\begin{matlab}
%}
  if i==1,
    ifOurCf2eps([mfilename 'ic'])
    disp('Type space character to animate simulation')
    pause
  else pause(0.05)
  end
end% i-loop over all times
ifOurCf2eps([mfilename 'fin'])
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:configPatches3fin}field~$u(x,y,z,t)$
at time $t=6$ of the patch scheme applied to the
heterogeneous wave~\pde\ with initial condition in
\cref{fig:configPatches3ic}.}
\includegraphics[scale=0.9]{configPatches3fin}
\end{figure}

Upon finishing execution of the example, exit this function.
\begin{matlab}
%}
return
end%if no arguments
%{
\end{matlab}


\IfFileExists{../Patch/heteroWave3.m}{\input{../Patch/heteroWave3.m}}{}






\begin{devMan}

\subsection{Parse input arguments and defaults}
\begin{matlab}
%}
p = inputParser;
fnValidation = @(f) isa(f, 'function_handle'); %test for fn name
addRequired(p,'fun',fnValidation); 
addRequired(p,'Xlim',@isnumeric);
%addRequired(p,'Dom'); % too flexible
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

Initially duplicate parameters for three space dimensions as
needed.
\begin{matlab}
%}
if numel(Xlim)==2,   Xlim = repmat(Xlim,1,3); end
if numel(nPatch)==1, nPatch = repmat(nPatch,1,3); end
if numel(dx)==1,     dx = repmat(dx,1,3); end
if numel(nSubP)==1,  nSubP = repmat(nSubP,1,3); end
%{
\end{matlab}

Check parameters.
\begin{matlab}
%}
assert(Xlim(1)<Xlim(2) ...
      ,'first pair of Xlim must be ordered increasing')
assert(Xlim(3)<Xlim(4) ...
      ,'second pair of Xlim must be ordered increasing')
assert(Xlim(5)<Xlim(6) ...
      ,'third pair of Xlim must be ordered increasing')
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
if size(Dom.type,1)==1, Dom.type=repmat(Dom.type,3,1); end
%{
\end{matlab}
Check what is and is not specified, and provide default of
Dirichlet boundaries if no \verb|bcOffset| specified when
needed.  Do so for all three directions independently.
\begin{matlab}
%}
patches.periodic=false;
for p=1:3
switch Dom.type(p,:)
case 'periodic'
    patches.periodic=true;
    if isfield(Dom,'bcOffset')
    warning('bcOffset not available for Dom.type = periodic'), end
    msg=' not available for Dom.type = periodic';
    if isfield(Dom,'X'), warning(['X' msg]), end
    if isfield(Dom,'Y'), warning(['Y' msg]), end
    if isfield(Dom,'Z'), warning(['Z' msg]), end
case {'equispace','chebyshev'}
    if ~isfield(Dom,'bcOffset'), Dom.bcOffset=zeros(2,3); end
    % for mixed with usergiven, following should still work
    if numel(Dom.bcOffset)==1
        Dom.bcOffset=repmat(Dom.bcOffset,2,3); end
    if numel(Dom.bcOffset)==3
        Dom.bcOffset=repmat(Dom.bcOffset(:)',2,1); end
    msg=' not available for Dom.type = equispace or chebyshev';
    if (p==1)& isfield(Dom,'X'), warning(['X' msg]), end
    if (p==2)& isfield(Dom,'Y'), warning(['Y' msg]), end
    if (p==3)& isfield(Dom,'Z'), warning(['Z' msg]), end
case 'usergiven'
%    if isfield(Dom,'bcOffset')
%    warning('bcOffset not available for usergiven Dom.type'), end
    msg=' required for Dom.type = usergiven';
    if p==1, assert(isfield(Dom,'X'),['X' msg]), end
    if p==2, assert(isfield(Dom,'Y'),['Y' msg]), end
    if p==3, assert(isfield(Dom,'Z'),['Z' msg]), end
otherwise 
    error([Dom.type ' is unknown Dom.type'])
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
for q=1:3
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
The periodic case is evenly spaced within the spatial
domain. Store the size ratio in \verb|patches|.
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
  end%for
  if DQmin<width*0.999999
    warning('too many Chebyshev patches (mid-domain overlap)')
  end%if
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
The user-given case is entirely up to a user to specify, we
just ensure it has the correct shape of a row.
\begin{matlab}
%}
case 'usergiven'
  if q==1, Q = reshape(Dom.X,1,[]); end
  if q==2, Q = reshape(Dom.Y,1,[]); end
  if q==3, Q = reshape(Dom.Z,1,[]); end
end%switch Dom.type
%{
\end{matlab}
Assign \(Q\)-coordinates to the correct spatial direction.
At this stage they are all rows.
\begin{matlab}
%}
if q==1, X=Q; end
if q==2, Y=Q; end
if q==3, Z=Q; end
end%for q
%{
\end{matlab}







\paragraph{Construct the micro-grids}
Construct the microscale in each patch.  Reshape the grid
to be 8D to suit dimensions (micro,Vars,Ens,macro).
\begin{matlab}
%}
nSubP = reshape(nSubP,1,3); % force to be row vector
assert(patches.EdgyInt | all(mod(nSubP,2)==1), ...
    'configPatches3: nSubP must be odd')
i0 = (nSubP(1)+1)/2;
patches.x = reshape( dx(1)*(-i0+1:i0-1)'+X ...
                   ,nSubP(1),1,1,1,1,nPatch(1),1,1);
i0 = (nSubP(2)+1)/2;
patches.y = reshape( dx(2)*(-i0+1:i0-1)'+Y ...
                   ,1,nSubP(2),1,1,1,1,nPatch(2),1);
i0 = (nSubP(3)+1)/2;
patches.z = reshape( dx(3)*(-i0+1:i0-1)'+Z ...
                   ,1,1,nSubP(3),1,1,1,1,nPatch(3));
%{
\end{matlab}



\paragraph{Pre-compute weights for macro-periodic} In the
case of macro-periodicity, precompute the weightings to
interpolate field values for coupling. \todo{Might sometime
extend to coupling via derivative values.}
\begin{matlab}
%}
if patches.periodic
  ratio = reshape(ratio,1,3); % force to be row vector
  patches.ratio = ratio; 
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
\item the right-face\slash centre realisations
\verb|1:nEnsem| are to interpolate to left-face~\verb|le|,
and 
\item the left-face\slash centre realisations
\verb|1:nEnsem| are to interpolate to~\verb|re|.
\end{itemize}
\verb|re| and \verb|li| are `transposes' of each other as
\verb|re(li)=le(ri)| are both \verb|1:nEnsem|. Similarly for
bottom-face\slash centre interpolation to top-face
via~\verb|to|, top-face\slash centre interpolation to
bottom-face via~\verb|bo|, back-face\slash centre
interpolation to front-face via~\verb|fr|, and
front-face\slash centre interpolation to back-face
via~\verb|ba|.

The default is nothing shifty.  This setting reduces the
number of if-statements in function \verb|patchEdgeInt3()|.
\begin{matlab}
%}
nE = patches.nEnsem;
patches.le = 1:nE;  patches.ri = 1:nE;  
patches.bo = 1:nE;  patches.to = 1:nE; 
patches.ba = 1:nE;  patches.fr = 1:nE; 
%{
\end{matlab}

However, if heterogeneous coefficients are supplied via
\verb|hetCoeffs|, then do some non-trivial replications.
First, get microscale periods, patch size, and replicate
many times in order to subsequently sub-sample: \verb|nSubP|
times should be enough. If \verb|cs| is more then 4D, then
the higher-dimensions are reshaped into the 4th dimension.
\begin{matlab}
%}
if ~isempty(cs)
  [mx,my,mz,nc] = size(cs);
  nx = nSubP(1); ny = nSubP(2); nz = nSubP(3);
  cs = repmat(cs,nSubP);
%{
\end{matlab}
If only one member of the ensemble is required, then
sub-sample to patch size, and store coefficients in
\verb|patches| as is.
\begin{matlab}
%}
  if nE==1, patches.cs = cs(1:nx-1,1:ny-1,1:nz-1,:); else
%{
\end{matlab}
But for $\verb|nEnsem|>1$ an ensemble of
$m_xm_ym_z$~phase-shifts of the coefficients is
constructed from the over-supply.  Here code phase-shifts
over the periods---the phase shifts are like
Hankel-matrices.
\begin{matlab}
%}
    patches.nEnsem = mx*my*mz;
    patches.cs = nan(nx-1,ny-1,nz-1,nc,mx,my,mz);
    for k = 1:mz
      ks = (k:k+nz-2);
      for j = 1:my
        js = (j:j+ny-2);
        for i = 1:mx
          is = (i:i+nx-2);
          patches.cs(:,:,:,:,i,j,k) = cs(is,js,ks,:);
        end
      end
    end
    patches.cs = reshape(patches.cs,nx-1,ny-1,nz-1,nc,[]);
%{
\end{matlab}
Further, set a cunning left\slash right\slash bottom\slash
top\slash front\slash back realisation of inter-patch
coupling.  The aim is to preserve symmetry in the system
when also invoking \verb|EdgyInt|.  What this coupling does
without \verb|EdgyInt| is unknown.  Use auto-replication.
\begin{matlab}
%}
    mmx=(0:mx-1)'; mmy=0:my-1; mmz=shiftdim(0:mz-1,-1);
    le = mod(mmx+mod(nx-2,mx),mx)+1;
    patches.le = reshape(  le+mx*(mmy+my*mmz)  ,[],1);
    ri = mod(mmx-mod(nx-2,mx),mx)+1;
    patches.ri = reshape(  ri+mx*(mmy+my*mmz)  ,[],1);
    bo = mod(mmy+mod(ny-2,my),my)+1;
    patches.bo = reshape( 1+mmx+mx*(bo-1+my*mmz) ,[],1);
    to = mod(mmy-mod(ny-2,my),my)+1;
    patches.to = reshape( 1+mmx+mx*(to-1+my*mmz) ,[],1);
    ba = mod(mmz+mod(nz-2,mz),mz)+1;
    patches.ba = reshape( 1+mmx+mx*(mmy+my*(ba-1)) ,[],1);
    fr = mod(mmz-mod(nz-2,mz),mz)+1;
    patches.fr = reshape( 1+mmx+mx*(mmy+my*(fr-1)) ,[],1);
%{
\end{matlab}
Issue warning if the ensemble is likely to be affected by
lack of scale separation.  \todo{Need to justify this and 
the arbitrary threshold more carefully??}
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
  spmd
%{
\end{matlab}
Second, decide which dimension is to be sliced among
parallel workers (for the moment, do not consider slicing
the ensemble). Choose the direction of most patches, biased
towards the last.
\begin{matlab}
%}
  [~,pari]=max(nPatch+0.01*(1:3));
  patches.codist=codistributor1d(5+pari);
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
    case 3, patches.z=codistributed(patches.z,patches.codist);
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
