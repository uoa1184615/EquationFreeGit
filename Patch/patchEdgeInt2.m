% patchEdgeInt2() provides the interpolation across 2D space
% for 2D patches of simulations of a smooth lattice system
% such as PDE discretisations.  AJR, Nov 2018 -- Dec 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt2()}: sets 2D patch
edge values from 2D macroscale interpolation}
\label{sec:patchEdgeInt2}


Couples 2D patches across 2D space by computing their edge
values via macroscale interpolation.  Research
\cite[]{Roberts2011a, Bunder2019c} indicates the patch
centre-values are sensible macroscale variables, and
macroscale interpolation of these determine patch-edge
values. However, for computational homogenisation in
multi-D, interpolating patch next-to-edge values appears
better \cite[]{Bunder2020a}. This function is primarily used
by \verb|patchSmooth2()| but is also useful for user
graphics. 

Script \verb|patchEdgeInt2test.m| verifies this code.
 
Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|) or otherwise via the global struct
\verb|patches|.
\begin{matlab}
%}
function u = patchEdgeInt2(u,patches)
if nargin<2, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
\(\verb|prod(nSubP)|  \cdot \verb|nVars| \cdot \verb|nEnsem|
\cdot \verb|prod(nPatch)|\) where there are \(\verb|nVars|
\cdot \verb|nEnsem|\) field values at each of the points in
the \(\verb|nSubP1| \cdot \verb|nSubP2| \cdot \verb|nPatch1|
\cdot \verb|nPatch2|\) grid on the \(\verb|nPatch1| \cdot
\verb|nPatch2|\) array of patches.

\item \verb|patches| a struct set by \verb|configPatches2()|
which includes the following information.
\begin{itemize}

\item \verb|.x| is \(\verb|nSubP1| \times1 \times1 \times1
\times \verb|nPatch1| \times1 \) array of the spatial
locations~\(x_{iI}\) of the microscale grid points in every
patch. Currently it \emph{must} be an equi-spaced lattice on
both macro- and microscales.

\item \verb|.y| is similarly \(1 \times \verb|nSubP2| \times1
\times1 \times1 \times \verb|nPatch2|\) array of the spatial
locations~\(y_{jJ}\) of the microscale grid points in every
patch. Currently it \emph{must} be an equi-spaced lattice on
both macro- and microscales.

\item \verb|.ordCC| is order of interpolation, currently
(Nov 2020) only \(\{0,2,4,\ldots\}\)

\item \verb|.stag| in \(\{0,1\}\) is one for staggered grid
(alternating) interpolation.

\item \verb|.Cwtsr| and \verb|.Cwtsl| define the coupling
coefficients for finite width interpolation in both the
\(x,y\)-directions.

\item \verb|.EdgyInt| true/false is true for interpolating
patch-edge values from opposite next-to-edge values (often
preserves symmetry).

\item \verb|.nEnsem| the number of realisations in the ensemble.

\item \verb|.parallel| whether serial or parallel.

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 6D array, \(\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nVars| \cdot \verb|nEnsem| \cdot
\verb|nPatch1| \cdot \verb|nPatch2|\), of the fields with
edge values set by interpolation (and corner vales set
to~\verb|NaN|).
\end{itemize}







\begin{devMan}

Test for reality of the field values, and define a function
accordingly.  Could be problematic if some variables are
real and some are complex, or if variables are of quite
different sizes. 
\begin{matlab}
%}
  if max(abs(imag(u(:))))<1e-9*max(abs(u(:)))
       uclean=@(u) real(u);
  else uclean=@(u) u; 
  end
%{
\end{matlab}

Determine the sizes of things. Any error arising in the
reshape indicates~\verb|u| has the wrong size.
\begin{matlab}
%}
[~,ny,~,~,~,Ny] = size(patches.y);
[nx,~,~,~,Nx,~] = size(patches.x);
nEnsem = patches.nEnsem;
nVars = round(numel(u)/numel(patches.x)/numel(patches.y)/nEnsem);
assert(numel(u) == nx*ny*Nx*Ny*nVars*nEnsem ...
  ,'patchEdgeInt2: input u has wrong size for parameters')
u = reshape(u,[nx ny nVars nEnsem Nx Ny ]);
%{
\end{matlab}
Get the size ratios of the patches in each direction.
\begin{matlab}
%}
rx = patches.ratio(1);
ry = patches.ratio(2);
%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, Dirichlet, Neumann, Robin?? These index
vectors point to patches and their four immediate
neighbours. 
\begin{matlab}
%}
I=1:Nx; Ip=mod(I,Nx)+1; Im=mod(I-2,Nx)+1;
J=1:Ny; Jp=mod(J,Ny)+1; Jm=mod(J-2,Ny)+1;
%{
\end{matlab}
The centre of each patch (as \verb|nx| and~\verb|ny| are
odd for centre-patch interpolation) is at indices
\begin{matlab}
%}
i0 = round((nx+1)/2);
j0 = round((ny+1)/2);
%{
\end{matlab}



\paragraph{Lagrange interpolation gives patch-edge values}
Compute centred differences of the mid-patch values for
the macro-interpolation, of all fields. Assumes the domain
is macro-periodic.
\begin{matlab}
%}
ordCC=patches.ordCC;
if ordCC>0 % then finite-width polynomial interpolation
%{
\end{matlab}
The patch-edge values are either interpolated from the next-to-edge values, or from the centre-cross values (not the patch-centre value itself as that seems to have worse properties in general).  Have not yet implemented core averages.
\begin{matlab}
%}
  if patches.EdgyInt % next-to-edge values    
    Ux = u([2 nx-1],2:(ny-1),:,:,I,J);
    Uy = u(2:(nx-1),[2 ny-1],:,:,I,J);
  else 
    Ux = u(i0,2:(ny-1),:,:,I,J);
    Uy = u(2:(nx-1),j0,:,:,I,J);
  end;
%{
\end{matlab}
Use finite difference formulas for the interpolation, so store finite differences (\(\mu\delta,\delta^2,\mu\delta^3,\delta^4,\ldots\)) in these arrays.  When parallel, in order to preserve the distributed array structure we use an index at the end for the differences.
\begin{matlab}
%}
  if patches.parallel
    dmux = zeros([size(Ux),ordCC],patches.codist); % 7D
    dmuy = zeros([size(Uy),ordCC],patches.codist); % 7D
  else
    dmux = zeros([size(Ux),ordCC]); % 7D
    dmuy = zeros([size(Uy),ordCC]); % 7D
  end
%{
\end{matlab}
First compute differences \(\mu\delta\) and \(\delta^2\) in both space directions.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
  else %disp('starting standard interpolation')   
    dmux(:,:,:,:,I,:,1) = (Ux(:,:,:,:,Ip,:) ...
                          -Ux(:,:,:,:,Im,:))/2; %\mu\delta 
    dmux(:,:,:,:,I,:,2) =  Ux(:,:,:,:,Ip,:) ...
       -2*Ux(:,:,:,:,I,:) +Ux(:,:,:,:,Im,:);    % \delta^2    
    dmuy(:,:,:,:,:,J,1) = (Uy(:,:,:,:,:,Jp) ...
                          -Uy(:,:,:,:,:,Jm))/2; %\mu\delta 
    dmuy(:,:,:,:,:,J,2) =  Uy(:,:,:,:,:,Jp) ...
       -2*Uy(:,:,:,:,:,J) +Uy(:,:,:,:,:,Jm);    % \delta^2
  end% if odd/even
%{
\end{matlab}
Recursively take \(\delta^2\) of these to form higher order
centred differences in both space directions.
\begin{matlab}
%}
   for k = 3:ordCC    
    dmux(:,:,:,:,I,:,k)   =     dmux(:,:,:,:,Ip,:,k-2) ...
      -2*dmux(:,:,:,:,I,:,k-2) +dmux(:,:,:,:,Im,:,k-2);    
    dmuy(:,:,:,:,:,J,k)   =     dmuy(:,:,:,:,:,Jp,k-2) ...
      -2*dmuy(:,:,:,:,:,J,k-2) +dmuy(:,:,:,:,:,Jm,k-2);
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for
each patch \cite[]{Roberts06d, Bunder2013b}, using weights
computed in \verb|configPatches2()|. Here interpolate to
specified order.

For the case where next-to-edge values interpolate to the
opposite edge-values: when we have an ensemble of
configurations, different configurations might be coupled to
each other, as specified by \verb|patches.le|,
\verb|patches.ri|, \verb|patches.to| and \verb|patches.bo|.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two edges
u(nx,2:(ny-1),:,patches.ri,I,:) ...
  = Ux(1,:,:,:,:,:)*(1-patches.stag) ...
    +sum( shiftdim(patches.Cwtsr(:,1),-6).*dmux(1,:,:,:,:,:,:) ,7);  
u(1 ,2:(ny-1),:,patches.le,I,:) ...
  = Ux(k,:,:,:,:,:)*(1-patches.stag) ...      
    +sum( shiftdim(patches.Cwtsl(:,1),-6).*dmux(k,:,:,:,:,:,:) ,7);
u(2:(nx-1),ny,:,patches.to,:,J) ...
  = Uy(:,1,:,:,:,:)*(1-patches.stag) ...
    +sum( shiftdim(patches.Cwtsr(:,2),-6).*dmuy(:,1,:,:,:,:,:) ,7);
u(2:(nx-1),1 ,:,patches.bo,:,J) ...
  = Uy(:,k,:,:,:,:)*(1-patches.stag) ...
    +sum( shiftdim(patches.Cwtsl(:,2),-6).*dmuy(:,k,:,:,:,:,:) ,7);
u([1 nx],[1 ny],:,:,:,:)=nan; % remove corner values
%{
\end{matlab}



\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic.
\begin{matlab}
%}
else% spectral interpolation
%{
\end{matlab}
We interpolate in terms of the patch index, \(j\)~say, not
directly in space. As the macroscale fields are
\(N\)-periodic in the patch index~\(j\), the macroscale
Fourier transform writes the centre-patch values as
\(U_j=\sum_{k}C_ke^{ik2\pi j/N}\). Then the edge-patch
values \(U_{j\pm r} =\sum_{k}C_ke^{ik2\pi/N(j\pm r)}
=\sum_{k}C'_ke^{ik2\pi j/N}\) where
\(C'_k=C_ke^{ikr2\pi/N}\). For \(N\)~patches we resolve
`wavenumbers' \(|k|<N/2\), so set row vector
\(\verb|ks|=k2\pi/N\) for `wavenumbers' \(\mathcode`\,="213B
k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)\) for
odd~\(N\), and \(\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, \pm(k_{\max}+1) -k_{\max}, \ldots, -1)\) for
even~\(N\).

Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches2|
tests there are an even number of patches). Then the
patch-ratio is effectively halved. The patch edges are near
the middle of the gaps and swapped.
\begin{matlab}
%}
 if patches.stag % transform by doubling the number of fields
 error('staggered grid not yet implemented??')
   v=nan(size(u)); % currently to restore the shape of u
   u=cat(3,u(:,1:2:nPatch,:),u(:,2:2:nPatch,:));
   stagShift=reshape(0.5*[ones(nVars,1);-ones(nVars,1)],1,1,[]);
   iV=[nVars+1:2*nVars 1:nVars]; % scatter interp to alternate field
   r=r/2;           % ratio effectively halved
   nPatch=nPatch/2; % halve the number of patches
   nVars=nVars*2;   % double the number of fields
 else % the values for standard spectral
    stagShift = 0;  
    iV = 1:nVars;
 end
%{
\end{matlab}
Now set wavenumbers in the two directions into two vectors
at the correct dimension.  In the case of even~\(N\) these
compute the \(+\)-case for the highest wavenumber zig-zag
mode, \(\mathcode`\,="213B k=(0,1, \ldots, k_{\max},
+(k_{\max}+1) -k_{\max}, \ldots, -1)\).
\begin{matlab}
%}
  kMax = floor((Nx-1)/2); 
  krx = shiftdim( rx*2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax) ,-3);
  kMay = floor((Ny-1)/2); 
  kry = shiftdim( ry*2*pi/Ny*(mod((0:Ny-1)+kMay,Ny)-kMay) ,-4);
%{
\end{matlab}

Compute the Fourier transform of the centre-cross values.
Unless doing patch-edgy interpolation when FT the
next-to-edge values.  If there are an even number of points,
then if complex, treat as positive wavenumber, but if real,
treat as cosine. When using an ensemble of configurations,
different configurations might be coupled to each other, as
specified by \verb|patches.le|, \verb|patches.ri|,
\verb|patches.to| and \verb|patches.bo|.
\begin{matlab}
%}
ix=(2:nx-1)';  iy=2:ny-1; % indices of interior
if ~patches.EdgyInt
     % here try central cross interpolation
     Cle = fft(fft(u(i0,iy,:,:,:,:),[],5),[],6);
     Cbo = fft(fft(u(ix,j0,:,:,:,:),[],5),[],6);
     Cri=Cle; Cto=Cbo;
else % edgyInt uses next-to-edge values
     Cle = fft(fft(u(   2,iy ,:,patches.le,:,:),[],5),[],6);
     Cri = fft(fft(u(nx-1,iy ,:,patches.ri,:,:),[],5),[],6);
     Cbo = fft(fft(u(ix,2    ,:,patches.bo,:,:),[],5),[],6);
     Cto = fft(fft(u(ix,ny-1 ,:,patches.to,:,:),[],5),[],6);
end     
%{
\end{matlab}
Now invert the double Fourier transforms to complete
interpolation. (Should stagShift be multiplied by rx/ry??)
Enforce reality when appropriate. 
\begin{matlab}
%}
u(nx,iy,:,:,:,:) = uclean( ifft(ifft( ...
    Cle.*exp(1i*(stagShift+krx))  ,[],5),[],6) );
u( 1,iy,:,:,:,:) = uclean( ifft(ifft( ...
    Cri.*exp(1i*(stagShift-krx))  ,[],5),[],6) );
u(ix,ny,:,:,:,:) = uclean( ifft(ifft( ...
    Cbo.*exp(1i*(stagShift+kry))  ,[],5),[],6) );
u(ix, 1,:,:,:,:) = uclean( ifft(ifft( ...
    Cto.*exp(1i*(stagShift-kry))  ,[],5),[],6) );
end% if spectral 
%{
\end{matlab}
Nan the corner values of every 2D patch.
\begin{matlab}
%}
u([1 nx],[1 ny],:,:,:,:)=nan; 
end% function patchEdgeInt2
%{
\end{matlab}
Fin, returning the 6D array of field values with
interpolated edges. 
\end{devMan} 
%}
