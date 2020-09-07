% Creates a data struct of the design of 3D patches for
% later use by the patch functions such as patchSmooth3() 
% AJR, Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{configPatches3()}: configures spatial
patches in 3D}
\label{sec:configPatches3}
\localtableofcontents





\subsection{Introduction}

Makes the struct~\verb|patches| for use by the patch\slash
gap-tooth time derivative\slash step function
\verb|patchSmooth3()|. 
\cref{sec:configPatches3eg,sec:homoDiffEdgy3} lists an
example of its use.
\begin{matlab}
%}
function configPatches3(fun,Xlim,BCs,nPatch,ordCC,ratio,nSubP ...
    ,varargin)
global patches
%{
\end{matlab}

\paragraph{Input}
If invoked with no input arguments, then executes an example
of simulating a heterogeneous wave \pde.
\begin{itemize}

\item \verb|fun| is the name of the user function,
\verb|fun(t,u,x,y,z)|, that computes time-derivatives (or
time-steps) of quantities on the 3D micro-grid within all
the 3D patches.

\item \verb|Xlim| array/vector giving the macro-space domain
of the computation: patches are distributed equi-spaced over
the interior of the rectangle \([\verb|Xlim(1)|,
\verb|Xlim(2)|] \times [\verb|Xlim(3)|, \verb|Xlim(4)|
\times [\verb|Xlim(5)|, \verb|Xlim(6)|]\): if \verb|Xlim| is
of length two, then the domain is the cubic domain of the
same interval in all three directions.

\item \verb|BCs| eventually and somehow will define the
macroscale boundary conditions.  Currently, \verb|BCs| is
ignored and the system is assumed macro-periodic in the
domain.

\item \verb|nPatch| determines the number of equi-spaced
spaced patches: if scalar, then use the same number of
patches in both directions, otherwise \verb|nPatch(1:3)|
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
both directions, otherwise \verb|ratio(1:3)| gives the ratio
in each direction.

\item \verb|nSubP| is the number of equi-spaced microscale
lattice points in each patch: if scalar, then use the same
number in both directions, otherwise \verb|nSubP(1:3)| gives
the number in each direction. If not using \verb|EdgyInt|,
then must be odd so that there is a central micro-grid 
point/planes in each patch.

\item \verb|'nEdge'| (not yet implemented), \emph{optional},
default=1, for each patch, the number of edge values set by
interpolation at the edge regions of each patch.  The
default is one (suitable for microscale lattices with only
nearest neighbour interactions).

\item \verb|'EdgyInt'|, true/false, \emph{optional},
default=false.  If true, then interpolate to left\slash
right\slash top\slash bottom\slash front\slash back
face-values from right\slash left\slash bottom\slash
top\slash back\slash front next-to-face values.  

\item \verb|'nEnsem'|,  \emph{optional-experimental},
default one, but if more, then an ensemble over this
number of realisations.

\item \verb|'hetCoeffs'|, \emph{optional}, default empty.
Supply a 3/4D array of microscale heterogeneous coefficients
to be used by the given microscale \verb|fun| in each patch.
Say the given array~\verb|cs| is of size \(m_x\times
m_y\times m_z\times n_c\), where \(n_c\)~is the number of
different arrays of coefficients. For example, in
heterogeneous diffusion, \(n_c=3\) for the diffusivities in
the \emph{three} different spatial directions.  The
coefficients are to be the same for each and every patch. 
However, macroscale variations are catered for by the
\(n_c\)~coefficients being \(n_c\)~parameters in some
macroscale formula.
\begin{itemize}

\item If \(\verb|nEnsem|=1\), then the array of coefficients
is just tiled across the patch size to fill up each patch,
starting from the \((1,1)\)-element.

\item If \(\verb|nEnsem|>1\) (value immaterial), then reset
\(\verb|nEnsem|:=m_x\cdot m_y\cdot m_z\) and construct an
ensemble of all \(m_x\cdot m_y\cdot m_z\) phase-shifts of
the coefficients. In this scenario, the inter-patch coupling
couples different members in the ensemble.  When
\verb|EdgyInt| is true, and when the coefficients are
diffusivities\slash elasticities in~\(x,y,z\) directions,
respectively, then this coupling cunningly preserves
symmetry.

\end{itemize}

\end{itemize}



\paragraph{Output} The \emph{global} struct \verb|patches|
is created and set with the following components.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x,y,z)| that computes the time derivatives (or
steps) on the patchy lattice. 

\item \verb|.ordCC| is the specified order of inter-patch
coupling. 

\item \verb|.stag| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling---not yet
implemented.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the
\(\verb|ordCC|\times 3\)-array of weights for the
inter-patch interpolation onto the right\slash top and
left\slash bottom edges (respectively) with patch:macroscale
ratio as specified.

\item \verb|.x| (8D) is \(\verb|nSubP(1)| \times1 \times1
\times1 \times1 \times \verb|nPatch(1)| \times1 \times1\)
array of the regular spatial locations~\(x_{ijk}\) of the
microscale grid points in every patch.  

\item \verb|.y| (8D) is \(1 \times \verb|nSubP(2)| \times1
\times1 \times1 \times1 \times \verb|nPatch(2)| \times1\)
array of the regular spatial locations~\(y_{ijk}\) of the
microscale grid points in every patch.  

\item \verb|.z| (8D) is \(1 \times1 \times \verb|nSubP(3)|
\times1 \times1 \times1 \times1 \times \verb|nPatch(3)|\)
array of the regular spatial locations~\(z_{ijk}\) of the
microscale grid points in every patch.  

\item \verb|.ratio| \(1\times 3\), are the size ratios of
every patch.

\item \verb|.nEdge| is, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch.

\item \verb|.le|, \verb|.ri|, \verb|.bo|, \verb|.to|,
\verb|.ba|, \verb|.fr| determine inter-patch coupling of
members in an ensemble. Each a column vector of
length~\verb|nEnsem|.

\item \verb|.cs| either
\begin{itemize}
\item \verb|[]| 0D, or 
\item if \(\verb|nEnsem|=1\), \((\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times (\verb|nSubP(3)|-1)\times n_c\) 4D
array of microscale heterogeneous coefficients, or
\item if \(\verb|nEnsem|>1\), \((\verb|nSubP(1)|-1)\times
(\verb|nSubP(2)|-1)\times (\verb|nSubP(3)|-1)\times
n_c\times m_xm_y\) 5D array of \(m_xm_y\)~ensemble of
phase-shifts of the microscale heterogeneous coefficients.
\end{itemize}

\item \verb|.codist|, \emph{optional}, describes the particular parallel distribution of arrays over the active parallel pool (the default is to activate the \emph{local} pool).  Can use \verb|isfield(patches,'codist')| to test whether parallel distribution has been invoked by a user.

\end{itemize}






\subsection{If no arguments, then execute an example}
\label{sec:configPatches3eg}
\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
The code here shows one way to get started: a user's script
may have the following three steps (arrows indicate function
recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches3 
\item ode23 integrator \into patchSmooth3 \into user's PDE
\item process results
\end{enumerate}

Set random heterogeneous
coefficients of period two in each of the three
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
mPeriod = [2 2 2];
cHetr = exp(0.3*randn([mPeriod 3]));
cHetr = cHetr*mean(1./cHetr(:)) 
%{
\end{matlab}
Establish global patch data struct to interface with a
function coding a nonlinear `diffusion' \pde: to be solved
on \([-\pi,\pi]^3\)-periodic domain, with \(5^3\) patches,
spectral interpolation~(\(0\)) couples the patches, each
patch of half-size ratio~\(0.4\) (relatively large for
visualisation), and with \(4^3\) points forming each
patch. 
\begin{matlab}
%}
configPatches3(@heteroWave3,[-pi pi], nan, 5 ...
    , 0, 0.35, mPeriod+2 ,'EdgyInt',true ...
    ,'hetCoeffs',cHetr);
%{
\end{matlab}
Set a  wave initial state using auto-replication of the
spatial grid, and as \cref{fig:configPatches3ic} shows.
\begin{matlab}
%}
u0 = 0.5+0.5*sin(patches.x+patches.y+patches.z);
v0 =    -0.5*cos(patches.x+patches.y+patches.z)*sqrt(3);
uv0 = cat(4,u0,v0);
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:configPatches3ic}initial
field~\(u(x,y,z,t)\) at time \(t=0\) of the patch scheme
applied to a heterogeneous wave~\pde:
\cref{fig:configPatches3fin} plots the computed field at
time \(t=6\).}
\includegraphics[scale=0.9]{configPatches3ic}
\end{figure}
Integrate in time using standard functions. In Matlab
\verb|ode15s| would be natural as the patch scheme is
naturally stiff, but \verb|ode23| is much quicker.
\begin{matlab}
%}
disp('Simulate heterogeneous wave u_tt=div[C*grad(u)]')
if ~exist('OCTAVE_VERSION','builtin')
    [ts,us] = ode23(@patchSmooth3,linspace(0,6),uv0(:));
else %disp('octave version is very slow for me')
    lsode_options('absolute tolerance',1e-4);
    lsode_options('relative tolerance',1e-4);
    [ts,us] = odeOcts(@patchSmooth3,[0 1 2],uv0(:));
end
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:configPatches3fin}.  Use \verb|patchEdgeInt3| to
interpolate patch-edge values (even if not drawn) in order
to most easily reconstruct the array data structure.

Replicate \(x\), \(y\), and \(z\) arrays to get individual
spatial coordinates of every data point.  Then, optionally,
set faces to \verb|nan| so the plot just shows
patch-interior data. 
\begin{matlab}
%}
figure(1), clf, colormap(0.8*jet)
xs = patches.x+0*patches.y+0*patches.z;
ys = patches.y+0*patches.x+0*patches.z;
zs = patches.z+0*patches.y+0*patches.x;
xs([1 end],:,:,:)=nan;  
xs(:,[1 end],:,:)=nan;  
xs(:,:,[1 end],:)=nan;  
j=find(~isnan(xs));
%{
\end{matlab}
The functions \verb|pix()| and \verb|col()| map the
\(u\)-data values to the size of the dots and to the colour
of the dots, respectively.
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
      scat(p) = scatter3(xs(j),ys(j),zs(j) ...
              ,pix(u(j)),col(u(j)),'filled');
      axis equal, caxis(col([0 1])), view(45-5*p,25)
      xlabel('x'), ylabel('y'), zlabel('z')
      title('view stereo pair cross-eyed')
    else % in matlab just update values
      set(scat(p),'CData',col(u(j)),'SizeData',pix(u(j)));
    end
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
\caption{\label{fig:configPatches3fin}field~\(u(x,y,z,t)\)
at time \(t=6\) of the patch scheme applied to the
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


\input{../Patch/heteroWave3.m}






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
addParameter(p,'hetCoeffs',[],@isnumeric);
addParameter(p,'parallel',false,@islogical);
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
parallel = p.Results.parallel;
%patches.nCore = p.Results.nCore;
%{
\end{matlab}

Initially duplicate parameters for both space dimensions as
needed.
\begin{matlab}
%}
if numel(Xlim)==2, Xlim = repmat(Xlim,1,3); end
if numel(nPatch)==1, nPatch = repmat(nPatch,1,3); end
if numel(ratio)==1, ratio = repmat(ratio,1,3); end
if numel(nSubP)==1, nSubP = repmat(nSubP,1,3); end
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
values.)  Store the size ratio in \verb|patches|.
\begin{matlab}
%}
ratio = reshape(ratio,1,3); % force to be row vector
patches.ratio=ratio; 
if ordCC>0, patchCwts(ratio,ordCC,patches.stag), end
%{
\end{matlab}


Third, set the centre of the patches in a the macroscale
grid of patches assuming periodic macroscale domain.
\begin{matlab}
%}
X = linspace(Xlim(1),Xlim(2),nPatch(1)+1);
DX = X(2)-X(1);
X = X(1:nPatch(1))+diff(X)/2;
Y = linspace(Xlim(3),Xlim(4),nPatch(2)+1);
DY = Y(2)-Y(1);
Y = Y(1:nPatch(2))+diff(Y)/2;
Z = linspace(Xlim(5),Xlim(6),nPatch(3)+1);
DZ = Z(2)-Z(1);
Z = Z(1:nPatch(3))+diff(Z)/2;
%{
\end{matlab}
Construct the microscale in each patch, assuming Dirichlet
patch edges, and a half-patch length of~\(\verb|ratio(1)|
\cdot \verb|DX|\), \(\verb|ratio(2)| \cdot \verb|DY|\)
and~\(\verb|ratio(3)| \cdot \verb|DZ|\), unless
\verb|patches.EdgyInt| is set in which case the patches are
of length \verb|ratio*DX+dx|, \verb|ratio*DY+dy| and
\verb|ratio*DZ+dz|.
\begin{matlab}
%}
nSubP = reshape(nSubP,1,3); % force to be row vector
assert(patches.EdgyInt | all(mod(nSubP,2)==1), ...
    'configPatches3: nSubP must be odd')
i0 = (nSubP(1)+1)/2;
if ~patches.EdgyInt, dx = ratio(1)*DX/(i0-1);
else                 dx = ratio(1)*DX/(nSubP(1)-2);
end
patches.x = dx*(-i0+1:i0-1)'+X;  % micro-grid
patches.x = reshape(patches.x,nSubP(1),1,1,1,1,nPatch(1),1,1);
i0 = (nSubP(2)+1)/2;
if ~patches.EdgyInt, dy = ratio(2)*DY/(i0-1);
else                 dy = ratio(2)*DY/(nSubP(2)-2);
end
patches.y = dy*(-i0+1:i0-1)'+Y;  % micro-grid
patches.y = reshape(patches.y,1,nSubP(2),1,1,1,1,nPatch(2),1);
i0 = (nSubP(3)+1)/2;
if ~patches.EdgyInt, dz = ratio(3)*DZ/(i0-1);
else                 dz = ratio(3)*DZ/(nSubP(3)-2);
end
patches.z = dz*(-i0+1:i0-1)'+Z;  % micro-grid
patches.z = reshape(patches.z,1,1,nSubP(3),1,1,1,1,nPatch(3));
%{
\end{matlab}
If parallel code, then first check a pool has started, or start one (presumably the `local' one). 
Second, decide which dimension is to be sliced among parallel workers.  Choose the direction of most patches, biased towards the last.
\begin{matlab}
%}
if parallel
  theparpool=gcp();
  [~,pari]=max(nPatch+0.01*(1:3));
  patches.codist=codistributor1d(5+pari);
else pari=0;
     if isfield(patches,'codist')
        rmfield(patches,'codist'); end
end
%{
\end{matlab}
\verb|patches.codist.Dimension| is the index that is split among workers.
Then distribute the appropriate coordinate direction among the workers.
\begin{matlab}
%}
switch pari
  case 0 % do nothing
  case 1, patches.x=codistributed(patches.x,patches.codist);
  case 2, patches.y=codistributed(patches.y,patches.codist);
  case 3, patches.z=codistributed(patches.z,patches.codist);
  otherwise
    error('should never have bad index for parallel distribution')
end
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
many times in order to subsequently sub-sample: 
\verb|nSubP| times should be enough.
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
\(m_xm_ym_z\)~phase-shifts of the coefficients is
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
top realisation of inter-patch coupling.  The aim is to
preserve symmetry in the system when also invoking
\verb|EdgyInt|.  What this coupling does without
\verb|EdgyInt| is unknown.  Use auto-replication.
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
