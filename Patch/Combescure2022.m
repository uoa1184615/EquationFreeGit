% For an example nonlinear elasticity in 1D, simulate and
% use MatCont to continue parametrised equilibria. An
% example of working via patches in space. Adapted from the
% example Figure 3(a) of Combescure(2022). AJR Nov 2022 --
% 13 Mar 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{Combescure2022}: simulation and
continuation of a 1D example nonlinear elasticity, via
patches}
\label{sec:Combescure2022}


Here we explore a nonlinear 1D elasticity problem with
complicated microstructure.  Executes a simulation. 
\emph{But the main aim is to show how one may use the 
MatCont continuation toolbox \cite[]{Govaerts2019} together 
with the Patch Scheme toolbox} \cite[]{Maclean2020a} in 
order to explore parameter space by continuing branches of
equilibria, etc. 

\begin{figure}
\centering
\caption{\label{fig:toyElas}1D arrangement of non-linear
springs with connections to (a) next-to-neighbour node 
\protect\cite[Fig.~3(a)]{Combescure2022}.  The blue box is
one micro-cell of one period, width~\(2b\), containing an 
odd and an even~\(i\).}
\setlength{\unitlength}{0.01\linewidth}
\begin{picture}(100,31)
\put(0,0){\framebox(100,31){}}
\put(0,0){\includegraphics[width=\linewidth]{Figs/toyElas}}
\put(36,4){\color{blue}\framebox(27,23){cell}}
\end{picture}
\end{figure}

\cref{fig:toyElas} shows the microscale elasticity---adapted
from Fig.~3(a) by \cite{Combescure2022}. Let the spatial
microscale lattice be at rest at points~\(x_i\), with
constant spacing~\(b\).  With displacement
variables~\(u_i(t)\), simulate the microscale lattice toy
elasticity system with 2-periodicity: for \(p=1,2\)
(respectively black and red in \cref{fig:toyElas}) and for
every~\(i\),
\begin{align}
&\epsilon^p_i:=\frac1{pb}(u_{i+p/2}-u_{i-p/2}),
&&\sigma^p_i:=w'_p(\epsilon^p_i),
\nonumber\\
&\DD t{u_{i}}= \sum_{p=1}^2\frac1{pb}(\sigma^p_{i+p/2}-\sigma^p_{i-p/2}),
&&w'_p(\epsilon):=\epsilon-M_p\epsilon^3+\epsilon^5.
\label{eq:heteroNLE}
\end{align}
The system has a microscale heterogeneity via the two
different functions~\(w'_p(\epsilon)\)
\cite[\S4]{Combescure2022}:
\begin{itemize}
\item microscale `instability' (structure) arises with
\(M_1:=2\) and \(M_2:=1\)
(\cref{fig:Comb22diffuSvis2b,fig:Comb22cpl}(b)); and
\item large scale `instability' (structure) arises with
\(M_1:=-1\) and \(M_2:=3\)
(\cref{fig:Comb22diffuLvis1,fig:Comb22cpl}(a)).
\end{itemize}

\paragraph{Microscale case} Set \(M_1:=2\) and \(M_2:=1\)\,.
 We fix the boundary conditions \(u(0)=0\) and parametrise
solutions by~\(u(L)\). There are equilibria \(u\approx
u(L)x/L\), but under large compression (large
negative~\(u(L)\)) interesting structures develop.
\cref{fig:Comb22diffuSvis2b} shows boundary layers with
microscale variations develop for \(u(L)<-13\). This figure
plots a strain~\(\epsilon\) as the strain is nearly constant
across the interior, so the boundary layers show up clearly.
As~\(u(L)\) decreases further, \cref{fig:Comb22diffuSvis2b}
shows the family of equilibria form complicated folds.
\cref{tblMicro} lists that MatCont also reports some branch
points and neutral saddle equilibria in this same regime
(see \cref{fig:Comb22cpl}(b)). I have not yet followed any
of the branches.

\begin{SCfigure}
\centering
\caption{\label{fig:Comb22diffuSvis2b}the case of microscale
`instability' appears as fluctuations close to both
boundaries.  As the system is physically compressed, the
equilibrium curve has complicated folds, as shown here (and
\cref{fig:Comb22cpl}(b)).}
\includegraphics[scale=0.8]{Comb22diffuSvis2b}
\end{SCfigure}


\begin{SCtable}
\centering\caption{\label{tblMicro}Interesting equilibria
for the cases of small scale instability:  \(M_1:=2\),
\(M_2:=1\) (\cref{fig:Comb22diffuSvis2b,fig:Comb22cpl}(b)). 
The rightmost column gives the \(-u(L)\)~parameter values
for corresponding critical points in the three-patch code
(\cref{fig:Comb22diffuSvis2N3}).}
\begin{tabular}{@{}rp{12.1em}r@{}}
\hline
$-u(L)$&MatCont description &\text{Patch}\\\hline
14.684 & Branch point &14.599\\
14.702 & Limit point &14.610\\
14.612 & Neutral Saddle Equilibrium &-\\
14.063 & Neutral Saddle Equilibrium &-\\
13.972 & Limit point &13.817\\
13.988 & Branch point &13.828\\
17.184 & Branch point &17.197\\
- & Limit point &17.227\\
17.183 & Neutral Saddle Equilibrium &17.211\\
%15.034 & Neutral Saddle Equilibrium \\
%15.024 & Limit point \\
%15.032 & Branch point \\
%17.987 & Branch point \\
%17.993 & Limit point \\
%17.987 & Neutral Saddle Equilibrium \\
\hline
\end{tabular}
\end{SCtable}

The previous paragraph's discussion is for a full domain
simulation, albeit done through an imposed computational
framework of physically abutting patches.
\cref{fig:Comb22diffuSvis2N3} shows the corresponding
MatCont continuation for the patch scheme with \(N=3\)
patches in the domain. Just three patches may well be
reasonable as the structures in this problem are the two
boundary layers, and a constant interior.
\cref{fig:Comb22diffuSvis2N3} shows the patch scheme
reasonably resolves these. \cref{tblMicro} also lists the
special points, as reported by MatCont, in the equilibria of
the patch scheme. The locations of these special points
reasonably match those found by the full domain simulation.

Importantly, MatCont is about \emph{ten times quicker to
execute on the patches} than on the full domain code. This
speed-up indicates that on larger scale problems the patch
scheme could be very useful in continuation explorations.

\begin{SCfigure}
\centering
\caption{\label{fig:Comb22diffuSvis2N3}using just three
patches, the case of microscale instability appears as
fluctuations close to both boundaries.  As the system is
physically compressed, the equilibrium curve has complicated
folds, as shown, and that approximately match
\cref{fig:Comb22diffuSvis2b}.  But it is computed ten times
quicker.}
\includegraphics[scale=0.8]{Comb22diffuSvis2N3}
\end{SCfigure}


\paragraph{Large scale case} Set \(M_1:=-1\) and
\(M_2:=3\)\,.  We fix the boundary conditions \(u(0)=0\) and
parametrise solutions by~\(u(L)\). There are equilibria
\(u\approx u(L)x/L\), but under large compression (large
negative~\(u(L)\)) interesting structures develop.
\cref{fig:Comb22diffuLvis1} shows an interior region of
higher magnitude strain develops. Again, this figure plots a
strain~\(\epsilon\) as the strain is nearly constant across
the domain, so the interior structure shows up clearly.
As~\(u(L)\) decreases further, \cref{fig:Comb22diffuLvis1}
shows the family of equilibria form complicated folds.
\cref{tblLarge} lists that MatCont also reports some branch
points and neutral saddle equilibria in this regime (see
\cref{fig:Comb22cpl}(a)). I have not yet followed any of the
branches.

\begin{SCfigure}
\centering
\caption{\label{fig:Comb22diffuLvis1}the case of large scale
`instability'.  Spatial structure appears in the middle of
the domain.  As the system is physically compressed, the
equilibrium curve has complicated folds, as shown here and
in \cref{fig:Comb22cpl}(a).}
\includegraphics[scale=0.8]{Comb22diffuLvis1}
\end{SCfigure}

\begin{SCtable}
\centering\caption{\label{tblLarge}Interesting equilibria
for the cases of large scale instability:  \(M_1:=-1\),
\(M_2:=3\) (\cref{fig:Comb22diffuLvis1,fig:Comb22cpl}(a)).}
\begin{tabular}{@{}rp{12.1em}r@{}}
\hline
$-u(L)$&MatCont description \\\hline
21.295 & Limit point \\
18.783 & Branch point \\
18.762 & Neutral Saddle Equilibrium \\
18.761 & Neutral Saddle Equilibrium \\
18.761 & Limit point \\
18.934 & Branch point \\
19.393 & Branch point \\
19.928 & Branch point \\
20.490 & Branch point \\
21.055 & Branch point \\
21.627 & Branch point \\
\hline
\end{tabular}
% these are from N=3 patches 
%21.469 & Branch point \\
%23.342 & Neutral Saddle Equilibrium \\
%23.462 & Branch point \\
%29.95 & Hopf \\
\end{SCtable}

The patch scheme with \(N=3\) patches does not make
reasonable predictions here.  I suspect this failure is
because the nontrivial interior structure here occupies too
much of the domain to fit into one `small' patch.  Here the
patch scheme may be useful if the physical domain is larger.



\subsection{Configure heterogeneous toy elasticity systems}
\label{sec:chtes}


Set some physical parameters.  Each cell is of
width~\(dx:=2b\) as I choose to store~\(u_i\) for odd~\(i\)
in \verb|u((i+1)/2,1,:)| and for even~\(i\) in
\verb|u(i/2,2,:)|, that is, the the physical displacements
form the  array
\begin{equation*}
\verb|u|=\begin{bmatrix} u_1&u_2\\ u_3&u_4\\ u_5&u_6\\ \vdots&\vdots \end{bmatrix}.
\end{equation*}
Then corresponding velocities are adjoined as 3rd and 4th column.
\begin{matlab}
%}
clear all
global b M vis
b = 1   % separation of lattice points
N = 42  % # lattice steps in L
L = b*N % length of domain
%{
\end{matlab}
The nonlinear coefficients of stress-strain are in
array~\verb|M|, chosen by~\verb|theCase|.
\begin{matlab}
%}
theCase = 2
switch theCase
case 1, M = [0 0 0 0]  % linear spring coefficients
case 2, M = [ 2 1 1 1] % micro scale instability??
case 3, M = [-1 3 1 1] % large scale instability??
end% switch
vis = 0.1 % does not appear to affect the equilibria
tEnd = 25
%{
\end{matlab}
Patch parameters: here \verb|nSubP| is the number of cells.
\begin{matlab}
%}
edgyInt = true
nSubP = 6, nPatch = 5 % gives full-domain on N=42, dx=2
%nSubP = 6, nPatch = 3 % patches for some crude comparison
%{
\end{matlab}


Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice elasticity
system~\cref{eq:heteroNLE}. Solved with \verb|nPatch|
patches, and interpolation (as high-order as possible) to
provide the edge-values of the inter-patch coupling
conditions.   
\begin{matlab}
%}
global patches
configPatches1(@heteroNLE,[0 L],'equispace',nPatch ...
    ,0,2*b,nSubP,'EdgyInt',edgyInt);
xx = patches.x+[-1 1]*b/2; % staggered sub-cell positions
%{
\end{matlab}





\subsection{Simulate in time}
Set the initial displacement and velocity of a simulation.
Integrate some time using standard integrator.
\begin{matlab}
%}
u0 = [ sin(pi/L*xx)  -0*0.14*cos(pi/L*xx) ];
tic
[ts,ust] = ode23(@patchSys1, tEnd*linspace(0,1,41), u0(:) ...
                ,[],patches,0);
cpuIntegrateTime = toc
%{
\end{matlab}

\paragraph{Plot space-time surface of the simulation} To see
the edge values of the patches, interpolate and then adjoin
a row of \verb|nan|s between patches. Because of the
odd/even storage we need to do a lot of permuting and
reshaping.   First, array of sub-cell coordinates in a
column for each patch, separating patches also by an extra
row of nans.
\begin{matlab}
%}
xs = reshape( permute( xx ,[2 1 3 4]), 2*nSubP,nPatch);  
xs(end+1,:) = nan;  
%{
\end{matlab}
Interpolate patch edge values, at all times simultaneously
by including time data into the 2nd dimension, and 2nd
reshaping it into the 3rd dimension.
\begin{matlab}
%}
uvs = reshape( permute( reshape(ust ...
      ,length(ts),nSubP,4,1,nPatch) ,[2 3 1 4 5]) ,nSubP,[],1,nPatch);
uvs = reshape( patchEdgeInt1(uvs) ,nSubP,4,[],nPatch);
%{
\end{matlab}
Extract displacement field, merge the 1st two columns,
permute the time variations to the 3rd, separate patches by
NaNs, and merge spatial data into the 1st column.
\begin{matlab}
%}
us = reshape( permute( uvs(:,1:2,:,:) ...
     ,[2 1 4 3]) ,2*nSubP,nPatch,[]);
us(end+1,:,:) = nan;
us = reshape(us,[],length(ts));
%{
\end{matlab}
Plot space-time surface of displacements over the macroscale
duration of the simulation.
\begin{matlab}
%}
  figure(1), clf()
  mesh(ts,xs(:),us) 
  view(60,40), colormap(0.8*jet), axis tight
  xlabel('time t'), ylabel('space x'), zlabel('u(x,t)') 
%{
\end{matlab}
Ditto for the velocity.
\begin{matlab}
%}
vs = reshape( permute( uvs(:,3:4,:,:) ...
     ,[2 1 4 3]) ,2*nSubP,nPatch,[]);
vs(end+1,:,:) = nan;
vs = reshape(vs,[],length(ts));
  figure(2), clf()
  mesh(ts,xs(:),vs) 
  view(60,40), colormap(0.8*jet), axis tight
  xlabel('time t'), ylabel('space x'), zlabel('v(x,t)') 
  drawnow
%{
\end{matlab}



\subsection{MatCont continuation}
First, use \verb|fsolve| to find an equilibrium at some
starting compressive displacement---a compression that
differs depending upon the case of nonlinearity.
\begin{matlab}
%}
muL0 = 12+6*(theCase==3)
u0 = [ -muL0*xx/L 0*xx ];
u0([1 end],:,:,:)=nan;
patches.i = find(~isnan(u0));
nVars=length(patches.i)
ueq=fsolve(@(v) dudtSys(0,v,muL0),u0(patches.i));
%{
\end{matlab}
Start search for equilibria at other compression parameters.
 Starting from zero, need 1000+ to find both the large-scale
and small-scale instability cases.   But need less points
when starting from parameter~\(12\) or so.
\begin{matlab}
%}
disp('Searching for equilibria, may take 1000+ secs')
[uv0,vec0]=init_EP_EP(@matContSys,ueq,muL0,[1]);
opt=contset; % initialise MatCont options
opt=contset(opt,'Singularities',true); %to report branch points, p.24
opt=contset(opt,'MaxNumPoints',400); % restricts how far matcont goes
opt=contset(opt,'Backward',true); % strangely, needs to go backwards??
[uv,vec,s,h,f]=cont(@equilibrium, uv0, [], opt); %MatCont continuation
%{
\end{matlab}

\paragraph{Post-process the report}
\begin{matlab}
%}
disp('List of interesting critical points')
muLs=uv(nVars+1,:);
for j=1:numel(s)
   disp([num2str(muLs(s(j).index),5) ' & ' s(j).msg ' \\'])
end
%{
\end{matlab}
Find a range of parameter and corresponding indices where
all the critical points occur.
\begin{matlab}
%}
p1=muLs(end); pe=muLs(1);
if numel(s)>3, for j=2:numel(s)-1
    p1=min(p1,muLs(s(j).index));
    pe=max(pe,muLs(s(j).index));
end, end
pMid=(p1+pe)/2, pWid=abs(pe-p1)
iPars=find(abs(muLs(:)-pMid)<pWid);%include some either side
%{
\end{matlab}
Choose an `evenly spaced' subset of the range so we only
plot up to  sixty of the parameter values reported in the
range.
\begin{matlab}
%}
nPars=numel(iPars)
dP=ceil((nPars-1)/60)
iP=1:dP:nPars;
muLP=muLs(iPars(iP));
%{
\end{matlab}
Interpolate patch edge values, at all parameters
simultaneously by including parameter-wise data into the 2nd
dimension, and 2nd reshaping it into the 3rd dimension.
\begin{matlab}
%}
uvs=nan(numel(iP),numel(u0));
uvs(:,patches.i)=uv(1:nVars,iPars(iP))';
uvs = reshape( permute( reshape(uvs ...
      ,length(muLP),nSubP,4,1,nPatch) ,[2 3 1 4 5]) ,nSubP,[],1,nPatch);
uvs = reshape( patchEdgeInt1(uvs) ,nSubP,4,[],nPatch);
%{
\end{matlab}
Extract displacement field, merge the 1st two columns,
permute the parameter variations to the 3rd, separate
patches by NaNs, and merge spatial data into the 1st column.
\begin{matlab}
%}
us = reshape( permute( uvs(:,1:2,:,:) ...
     ,[2 1 4 3]) ,2*nSubP,nPatch,[]);
us(end+1,:,:) = nan;
us = reshape(us,[],length(muLP));
%{
\end{matlab}
Plot space-time surface of displacements over the macroscale
duration of the simulation.
\begin{matlab}
%}
  figure(4), clf()
  mesh(muLP,xs(:),us) 
  view(60,40), colormap(0.8*jet), axis tight
  xlabel('-u(L)'), ylabel('space x'), zlabel('u(x)') 
%{
\end{matlab}
Plot space-time surface of strain, differences in
displacements, over the parameter variation.
\begin{matlab}
%}
  figure(5), clf()
  mesh(muLP,xs(1:end-1),diff(us)) 
  view(45,20), colormap(0.8*jet), axis tight
  xlabel('-u(L)'), ylabel('space x'), zlabel('strain \delta u(x)') 
  ifOurCf2eps(['Comb22diffu' num2str(theCase)],[12 9])%optionally save
%{
\end{matlab}


\begin{figure}
\centering
\caption{\label{fig:Comb22cpl}cross-sections through
\cref{fig:Comb22diffuLvis1,fig:Comb22diffuSvis2b}: (a)~large
scale case, at the mid-point in space of
\cref{fig:Comb22diffuLvis1}; (b)~microscale case, in a
boundary layer of \cref{fig:Comb22diffuSvis2b}.  These
cross-sections are labelled with the various critical
points.}
\begin{tabular}{@{}cc@{}}
(a) large scale case & (b) microscale case\\
\includegraphics[scale=0.75]{Figs/Comb22cpl3}&
\includegraphics[scale=0.75]{Figs/Comb22cpl2}
\end{tabular}
\end{figure}

\paragraph{Labelled parameter plot} Get the labelled 2D
plots of \cref{fig:Comb22cpl} via MatCont's \verb|cpl|
function.  In high-D problems it is unlikely that any one
variable is a good thing to plot, so I show how to plot
something else, here a strain. I use all the computed points
so reform~\verb|uvs| (possibly better to have merged the
critical points into the list of plotted parameters??).  
\begin{matlab}
%}
uvs = nan(numel(muLs),numel(u0));
uvs(:,patches.i) = uv(1:nVars,:)';
uvs = reshape( uvs ,[],nSubP,4,nPatch);
%{
\end{matlab}
As a function of  the parameter, plot the strain in the
middle of the domain (the middle of the middle patch),
unless it is the microscale case when we plot a strain near
the middle of the left boundary layer.
\begin{matlab}
%}
if theCase==2, thePatch=1;
else thePatch=(nPatch+1)/2;
end%if
  figure(7),clf
  du = diff( uvs(:,nSubP/2,1:2,thePatch) ,1,3);
  cpl([muLs;du'],[],s);
  xlabel('-u(L)')
  if thePatch==1,  ylabel('boundary layer strain')
  else ylabel('mid-domain strain')
  end
  ifOurCf2eps(['Comb22cpl' num2str(theCase)],[9 7])%optionally save
%{
\end{matlab}





\subsection{\texttt{matContSys}: basic function for MatCont analysis}
This is the simple `odefile' of the patch scheme wrapped
around the microcode.
\begin{matlab}
%}
function out = matContSys%(t,coordinates,flag,y,z)
out{1} = [];%@init;
out{2} = @dudtSys;
out{3} = [];%@jacobian;
out{4} = [];%@jacobianp;
out{5} = [];%@hessians;
out{6} = [];%@hessiansp;
out{7} = [];
out{8} = [];
out{9} = [];
end% function matContSys
%{
\end{matlab}



\subsection{\texttt{dudtSys()}: wraps around the patch wrapper}
This function adjoins \verb|patches| to the argument list,
places the variables within the patch structure, and then
extracts their time derivatives to return. Used by both
MatCont and \verb|fsolve|.
\begin{matlab}
%}
function ut = dudtSys(t,u,p)
global patches
%{
\end{matlab}
The 4 here is the number of variables in each micro-cell, 
that is, notionally `at' each \(x\)-grid point.
\begin{matlab}
%}
U=nan(1,4)+patches.x; 
U(patches.i)=u(:);
Ut=patchSys1(t,U,patches,p);
ut=Ut(patches.i);
end
%{
\end{matlab}




\subsection{\texttt{heteroNLE()}: forced heterogeneous elasticity}
\label{sec:heteroNLE}

This function codes the lattice heterogeneous example
elasticity inside the patches.  Computes the time derivative
at each point in the interior of a patch, output
in~\verb|uvt|.
\begin{matlab}
%}
function uvt = heteroNLE(t,uv,patches,muL)
  if nargin<4, muL=0; end% default end displacement is zero
  global b M vis 
%{
\end{matlab}
Separate state vector into displacement and velocity fields:
\(u_{ijI}\)~is the displacement at the \(j\)th~point in the
\(i\)th 2-cell in the \(I\)th~patch; similarly for
velocity~\(v_{ijI}\).  That is, physically neighbouring
points have different~\(j\), whereas physical
next-to-neighbours have \(i\)~different by one.
\begin{matlab}
%}
  u=uv(:,1:2,:,:); v=uv(:,3:4,:,:); % separate u and v=du/dt
%{
\end{matlab}

Provide boundary conditions, here fixed displacement and
velocity in the left/right sub-cells of the
leftmost/rightmost patches.
\begin{matlab}
%}
u(1,:,:,1)=0;
v(1,:,:,1)=0;
u(end,:,:,end)=-muL;
v(end,:,:,end)=0;
%{
\end{matlab}

Compute the two different strain fields, and also a first
derivative for some optional viscosity.
\begin{matlab}
%}
  eps2 = diff(u)/(2*b);
  eps1 = [u(:,2,:,:)-u(:,1,:,:) u([2:end 1],1,:,:)-u(:,2,:,:)]/b;
  eps1(end,2,:,:)=nan; % as this value is fake
  vx1  = [v(:,2,:,:)-v(:,1,:,:) v([2:end 1],1,:,:)-v(:,2,:,:)]/b;
  vx1(end,2,:,:)=nan; % as this value is fake
%{
\end{matlab}
Set corresponding nonlinear stresses
\begin{matlab}
%}
  sig2 = eps2-M(2)*eps2.^3+M(4)*eps2.^5;
  sig1 = eps1-M(1)*eps1.^3+M(3)*eps1.^5;
%{
\end{matlab}
Preallocate output array, and fill in time derivatives of
displacement and velocity, from velocity and gradient of
stresses, respectively.
\begin{matlab}
%}
  uvt = nan+uv;          % preallocate output array
  i=2:size(uv,1)-1;
  % rate of change of position
  uvt(i,1:2,:,:) = v(i,:,:,:);  
  % rate of change of velocity +some artificial viscosity??
  uvt(i,3:4,:,:) = diff(sig2) ...
    +[ sig1(i,1,:,:)-sig1(i-1,2,:,:)  diff(sig1(i,:,:,:),1,2)] ... 
  +vis*[ vx1(i,1,:,:)-vx1(i-1,2,:,:)  diff(vx1(i,:,:,:),1,2) ]; 

end% function heteroNLE
%{
\end{matlab}
%}
