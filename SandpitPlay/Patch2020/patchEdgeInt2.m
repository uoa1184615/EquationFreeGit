% Provides the interpolation across 2D space for 2D patches
% of simulations of a smooth lattice system such as PDE
% discretisations.  AJR, Nov 2018 -- 15 Apr 2020 -- July 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section[\texttt{patchEdgeInt2()}: 2D patch edge values from
2D interpolation] {\texttt{patchEdgeInt2()}: sets 2D patch
edge values from 2D macroscale interpolation}
\label{sec:patchEdgeInt2}


Couples 2D patches across 2D space by computing their edge
values via macroscale interpolation.  Assumes that the patch
centre-values are sensible macroscale variables, and patch
edge values are determined by macroscale interpolation of
the patch-centre values??  Communicate patch-design variables
via the global struct~\verb|patches|.
\begin{matlab}
%}
function u = patchEdgeInt2(u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector of length \(\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nVars|\cdot \verb|nEnsem| \cdot
\verb|nPatch1| \cdot \verb|nPatch2|\) where there are
\verb|nVars| field values at each of the points in the
\(\verb|nSubP1| \cdot \verb|nSubP2| \cdot \verb|nPatch1|
\cdot \verb|nPatch2|\) grid on the \(\verb|nPatch1| \cdot
\verb|nPatch2|\) array of patches.

\item \verb|patches| a struct set by \verb|configPatches2()|
which includes the following information.
\begin{itemize}

\item \verb|.x| is \(\verb|nSubP1| \times1 \times1 \times1
\times \verb|nPatch1| \times1 \) array of the spatial
locations~\(x_{ij}\) of the microscale grid points in every
patch. Currently it \emph{must} be an equi-spaced lattice on
both macro- and microscales.

\item \verb|.y| is similarly \(\times1 \verb|nSubP1| \times1
\times1 \times1 \times \verb|nPatch1|\) array of the spatial
locations~\(y_{ij}\) of the microscale grid points in every
patch. Currently it \emph{must} be an equi-spaced lattice on
both macro- and microscales.

\item \verb|.ordCC| is order of interpolation, currently
(July 2020) only \(\{0,2,4,\ldots\}\)

\item \verb|.stag| in \(\{0,1\}\) is one for staggered grid
(alternating) interpolation.

\item \verb|.Cwtsr| and \verb|.Cwtsl| define the coupling
coefficients for finite width interpolation.

\item \verb|.EdgyInt| true/false is true for interpolating
patch-edge values from opposite next-to-edge values (often
preserves symmetry).

\end{itemize}
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item \verb|u| is \(\verb|nSubP1| \cdot \verb|nSubP2| \cdot
\verb|nVars|\cdot \verb|nEnsem| \cdot \verb|nPatch1| \cdot
\verb|nPatch2|\) 6D array of the fields with edge values set
by interpolation.
\end{itemize}







\begin{devMan}

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
%  nSubP=[nx ny], nPatch=[Nx Ny], nVars=nVars, sizeu=size(u)
u = reshape(u,[nx ny nVars nEnsem Nx Ny ]);
%{
\end{matlab}
With Dirichlet patches, the half-length of a patch is
\(h=dx(n_\mu-1)/2\) (or~\(-2\) for specified flux), unless
we are interpolating from next-to-edge values, and the
ratio needed for interpolation is then \(r=h/\Delta X\).
Compute lattice sizes from inside the patches as the edge
values may be \nan{}s, etc.
\begin{matlab}
%}
dx = patches.x(3,1,1,1,1,1)-patches.x(2,1,1,1,1,1);
DX = patches.x(2,1,1,1,2,1)-patches.x(2,1,1,1,1,1);
if patches.EdgyInt==0, rx = dx*(nx-1)/2/DX;
    else               rx = dx*(nx-2)/DX;
    end
dy = patches.y(1,3,1,1,1,1)-patches.y(1,2,1,1,1,1);
DY = patches.y(1,2,1,1,1,2)-patches.y(1,2,1,1,1,1);
if patches.EdgyInt==0, ry = dy*(ny-1)/2/DY;
    else               ry = dy*(ny-2)/DY;
    end

%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, Dirichlet, Neumann, Robin?? These index vectors
point to patches and their two immediate neighbours---currently 
not needed. 
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
Currently, only next-to-edge interpolation is implemented.
\begin{matlab}
%}
ordCC=patches.ordCC;
if ordCC>0 % then finite-width polynomial interpolation
  if patches.EdgyInt % next-to-edge values    
    uCorex = u([2 nx-1],2:(ny-1),:,:,I,J);
    uCorey = u(2:(nx-1),[2 ny-1],:,:,I,J);
  else 
    %disp('polynomial interpolation, currently couple from the centre-cross??')
    uCorex = u(i0,2:(ny-1),:,:,I,J);
    uCorey = u(2:(nx-1),j0,:,:,I,J);
  end;
  dmux = zeros([ordCC,size(uCorex)]); % 7D
  dmuy = zeros([ordCC,size(uCorey)]); % 7D
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet implemented for alternative (odd) patch coupling')
  else %disp('starting standard interpolation')   
    dmux(1,:,:,:,:,I,:) = (uCorex(:,:,:,:,Ip,:) -uCorex(:,:,:,:,Im,:))/2; % \mu\delta 
    dmux(2,:,:,:,:,I,:) = (uCorex(:,:,:,:,Ip,:) -2*uCorex(:,:,:,:,I,:) +uCorex(:,:,:,:,Im,:)); % \delta^2    
    dmuy(1,:,:,:,:,:,J) = (uCorey(:,:,:,:,:,Jp) -uCorey(:,:,:,:,:,Jm))/2; % \mu\delta 
    dmuy(2,:,:,:,:,:,J) = (uCorey(:,:,:,:,:,Jp) -2*uCorey(:,:,:,:,:,J) +uCorey(:,:,:,:,:,Jm)); % \delta^2
  end% if odd/even
%{
\end{matlab}
Recursively take \(\delta^2\) of these to form higher order
centred differences.
\begin{matlab}
%}
   for k = 3:ordCC    
    dmux(k,:,:,:,:,I,:) = dmux(k-2,:,:,:,:,Ip,:) ...
    -2*dmux(k-2,:,:,:,:,I,:) +dmux(k-2,:,:,:,:,Im,:);    
    dmuy(k,:,:,:,:,:,J) = dmuy(k-2,:,:,:,:,:,Jp) ...
    -2*dmuy(k-2,:,:,:,:,:,J) +dmuy(k-2,:,:,:,:,:,Jm);
  end
%dmux33=dmux(:,:,:,:,:,3,3),sizedmux33=size(dmux33)
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for
each patch \cite[]{Roberts06d, Bunder2013b}, using weights
computed in \verb|configPatches2()|. Here interpolate to
specified order.

Where next-to-edge values interpolate to the opposite
edge-values. When we have an ensemble of configurations,
different configurations might be coupled to each other, as
specified by \verb|patches.le|, \verb|patches.ri|,
\verb|patches.to| and \verb|patches.bo|.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two edges
u(nx,2:(ny-1),:,patches.ri,I,:) ...
  = uCorex(1,:,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsr(:,1).*dmux(:,1,:,:,:,:,:) ),1);  
u(1 ,2:(ny-1),:,patches.le,I,:) ...
  = uCorex(k,:,:,:,:,:)*(1-patches.stag) ...      
    +shiftdim(sum( patches.Cwtsl(:,1).*dmux(:,k,:,:,:,:,:) ),1);
u(2:(nx-1),ny,:,patches.to,:,J) ...
  = uCorey(:,1,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsr(:,2).*dmuy(:,:,1,:,:,:,:) ),1);
u(2:(nx-1),1 ,:,patches.bo,:,J) ...
  = uCorey(:,k,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsl(:,2).*dmuy(:,:,k,:,:,:,:) ),1);
u([1 nx],[1 ny],:,:,:,:)=nan; % remove corner values
%{
\end{matlab}



\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic. We interpolate in
terms of the patch index~\(j\), say, not directly in space. 
As the macroscale fields are \(N\)-periodic in the patch
index~\(j\), the macroscale Fourier transform writes the
centre-patch values as \(U_j=\sum_{k}C_ke^{ik2\pi j/N}\).
Then the edge-patch values \(U_{j\pm r}
=\sum_{k}C_ke^{ik2\pi/N(j\pm r)} =\sum_{k}C'_ke^{ik2\pi
j/N}\) where \(C'_k=C_ke^{ikr2\pi/N}\). For \(N\)~patches we
resolve `wavenumbers' \(|k|<N/2\), so set row vector
\(\verb|ks|=k2\pi/N\) for `wavenumbers' \(\mathcode`\,="213B
k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)\) for
odd~\(N\), and \(\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, \pm(k_{\max}+1) -k_{\max}, \ldots, -1)\) for
even~\(N\).
\begin{matlab}
%}
else% spectral interpolation
%{
\end{matlab}
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
Compute the Fourier transform of the patch centre-values for
all the fields.  Unless doing patch-edgy interpolation when
FT the next-to-edge values.  If there are an even number of
points, then if complex, treat as positive wavenumber, but
if real, treat as cosine. When using an ensemble of
configurations, different configurations might be coupled to
each other, as specified by \verb|patches.le|,
\verb|patches.ri|, \verb|patches.to| and \verb|patches.bo|.
\begin{matlab}
%}
ix=(2:nx-1)';  iy=2:ny-1; % indices of interior
if ~patches.EdgyInt
     Cle = fft(fft(u(i0,j0,:,:,:,:),[],5),[],6); 
     Cbo = Cle;
else 
     Cle = fft(fft(u(   2,iy ,:,patches.le,:,:),[],5),[],6);
     Cri = fft(fft(u(nx-1,iy ,:,patches.ri,:,:),[],5),[],6);
     Cbo = fft(fft(u(ix,2    ,:,patches.bo,:,:),[],5),[],6);
     Cto = fft(fft(u(ix,ny-1 ,:,patches.to,:,:),[],5),[],6);
end     
% fill in the cross of Fourier-shifted mid-values
% something strange here as both krx and kry are row vectors
if ~patches.EdgyInt 
  % y-fraction of kry along left/right edges
  ks = (iy-j0)*2/(ny-1).*kry; 
  Cle = Cle.*exp(1i*ks); 
  Cri = Cle;
  % x-fraction of krx along bottom/top edges
  ks = (ix-i0)*2/(nx-1).*krx;
  Cbo = Cbo.*exp(1i*ks); 
  Cto = Cbo;
end
u(nx,iy,:,:,:,:) = uclean( ifft(ifft( ...
    Cle.*exp(1i*(stagShift+krx))  ,[],5),[],6) );
u( 1,iy,:,:,:,:) = uclean( ifft(ifft( ...
    Cri.*exp(1i*(stagShift-krx))  ,[],5),[],6) );
u(ix,ny,:,:,:,:) = uclean( ifft(ifft( ...
    Cbo.*exp(1i*(stagShift+kry))  ,[],5),[],6) );
u(ix, 1,:,:,:,:) = uclean( ifft(ifft( ...
    Cto.*exp(1i*(stagShift-kry))  ,[],5),[],6) );
u([1 nx],[1 ny],:,:,:)=nan; % remove corner values

end% if spectral 
end% function patchEdgeInt2
%{
\end{matlab}
Fin, returning the 6D array of field values with
interpolated edges. 
\end{devMan} 
%}
