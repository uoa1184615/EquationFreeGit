% Provides the interpolation across 3D space for 3D patches
% of simulations of a smooth lattice system such as PDE
% discretisations.  AJR, Aug--Sep 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section[\texttt{spmdPatchEdgeInt3()}: 3D patch face values from
3D interpolation] {\texttt{spmdPatchEdgeInt3()}: sets 3D patch
face values from 3D macroscale interpolation}
\label{sec:spmdPatchEdgeInt3}

Couples 3D patches across 3D space by computing their face
values via macroscale interpolation.  Assumes that the patch
centre-values are sensible macroscale variables, and patch
face values are determined by macroscale interpolation of
the patch-centre values??  Communicate patch-design variables
via the global struct~\verb|patches|.
\begin{matlab}
%}
function u = spmdPatchEdgeInt3(u,patches)
%global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector of length \(\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nSubP3| \cdot \verb|nVars|\cdot
\verb|nEnsem| \cdot \verb|nPatch1| \cdot \verb|nPatch2|
\cdot \verb|nPatch3|\) where there are \verb|nVars| field
values at each of the points in the \(\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nSubP3| \cdot \verb|nPatch1| \cdot
\verb|nPatch2| \cdot \verb|nPatch3|\) grid on the
\(\verb|nPatch1| \cdot \verb|nPatch2| \cdot \verb|nPatch3|\)
array of patches.

\item \verb|patches| a struct set by \verb|configPatches3()|
which includes the following information.
\begin{itemize}

\item \verb|.x| is \(\verb|nSubP1| \times1 \times1 \times1
\times1 \times \verb|nPatch1| \times1 \times1 \) array of
the spatial locations~\(x_{ijk}\) of the microscale grid
points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.y| is similarly \(\times1 \verb|nSubP2| \times1
\times1 \times1 \times1 \times \verb|nPatch2| \times1\)
array of the spatial locations~\(y_{ijk}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.z| is similarly \(\times1 \times1 \verb|nSubP3|
\times1 \times1 \times1 \times1 \times \verb|nPatch3|\)
array of the spatial locations~\(z_{ijk}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.ordCC| is order of interpolation, currently
(Aug 2020) only \(\{0,2,4,\ldots\}\)

\item \verb|.stag| in \(\{0,1\}\) is one for staggered grid
(alternating) interpolation.

\item \verb|.Cwtsr| and \verb|.Cwtsl| define the coupling
coefficients for finite width interpolation.

\item \verb|.EdgyInt| true/false is true for interpolating
patch-face values from opposite next-to-face values (often
preserves symmetry).

\end{itemize}
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item \verb|u| is \(\verb|nSubP1| \cdot \verb|nSubP2| \cdot
\verb|nSubP3| \cdot \verb|nVars|\cdot \verb|nEnsem| \cdot
\verb|nPatch1| \cdot \verb|nPatch2| \cdot \verb|nPatch3|\)
8D array of the fields with face values set by
interpolation.
\end{itemize}







\begin{devMan}

Test for reality of the field values, and define a function
accordingly.  Could be problematic if some variables are
real and some are complex, or if variables are of quite
different sizes.
Have to do such function definition outside of \verb|spmd|-block.
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
%spmd%%%%%%%%%%
%warning('spmdPatchEdgeInt3 before patches.z')
%patch=patches
[~,~,nz,~,~,~,~,Nz] = size(patches.z);
%warning('spmdPatchEdgeInt3 before patches.y')
[~,ny,~,~,~,~,Ny,~] = size(patches.y);
%warning('spmdPatchEdgeInt3 before patches.x')
[nx,~,~,~,~,Nx,~,~] = size(patches.x);
%warning('spmdPatchEdgeInt3 before nEnsem')
nEnsem = patches.nEnsem;
%warning('spmdPatchEdgeInt3 before nVars')
nVars = round( numel(u)/numel(patches.x) ...
    /numel(patches.y)/numel(patches.z)/nEnsem );
numelu=numel(u);
assert(numel(u) == nx*ny*nz*Nx*Ny*Nz*nVars*nEnsem ...
  ,'spmdPatchEdgeInt3: input u has wrong size for parameters')
u = reshape(u,[nx ny nz nVars nEnsem Nx Ny Nz]);
%{
\end{matlab}
Get the size ratios of the patches in each direction.
\begin{matlab}
%}
rx = patches.ratio(1);
ry = patches.ratio(2);
rz = patches.ratio(3);
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
K=1:Nz; Kp=mod(K,Nz)+1; Km=mod(K-2,Nz)+1;
%{
\end{matlab}
The centre of each patch (as \verb|nx| and~\verb|ny| are
odd for centre-patch interpolation) is at indices
\begin{matlab}
%}
i0 = round((nx+1)/2);
j0 = round((ny+1)/2);
k0 = round((nz+1)/2);
%{
\end{matlab}


\paragraph{Lagrange interpolation gives patch-face values}
Compute centred differences of the mid-patch values for
the macro-interpolation, of all fields. Assumes the domain
is macro-periodic.
Currently, only next-to-face interpolation is implemented.
\begin{matlab}
%}
ordCC=patches.ordCC;
if ordCC>0 % then finite-width polynomial interpolation
  if patches.EdgyInt % next-to-face values    
    uCorex = u([2 nx-1],2:(ny-1),2:(nz-1),:,:,I,J,K);
    uCorey = u(2:(nx-1),[2 ny-1],2:(nz-1),:,:,I,J,K);
    uCorez = u(2:(nx-1),2:(ny-1),[2 nz-1],:,:,I,J,K);
  else 
    %disp('currently couple from the mid-planes??')
    uCorex = u(i0,2:(ny-1),2:(nz-1),:,:,I,J,K);
    uCorey = u(2:(nx-1),j0,2:(nz-1),:,:,I,J,K);
    uCorez = u(2:(nx-1),2:(ny-1),k0,:,:,I,J,K);
  end;
  dmux = zeros([ordCC,size(uCorex)]); % 9D
  dmuy = zeros([ordCC,size(uCorey)]); % 9D
  dmuz = zeros([ordCC,size(uCorez)]); % 9D
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
  else %disp('starting standard interpolation')   
    dmux(1,:,:,:,:,:,I,:,:) = (uCorex(:,:,:,:,:,Ip,:,:) ...
    -uCorex(:,:,:,:,:,Im,:,:))/2; % \mu\delta 
    dmux(2,:,:,:,:,:,I,:,:) = (uCorex(:,:,:,:,:,Ip,:,:) ...
    -2*uCorex(:,:,:,:,:,I,:,:) +uCorex(:,:,:,:,:,Im,:,:)); % \delta^2    
    dmuy(1,:,:,:,:,:,:,J,:) = (uCorey(:,:,:,:,:,:,Jp,:) ...
    -uCorey(:,:,:,:,:,:,Jm,:))/2; % \mu\delta 
    dmuy(2,:,:,:,:,:,:,J,:) = (uCorey(:,:,:,:,:,:,Jp,:) ...
    -2*uCorey(:,:,:,:,:,:,J,:) +uCorey(:,:,:,:,:,:,Jm,:)); % \delta^2
    dmuz(1,:,:,:,:,:,:,:,K) = (uCorez(:,:,:,:,:,:,:,Kp) ...
    -uCorez(:,:,:,:,:,:,:,Km))/2; % \mu\delta 
    dmuz(2,:,:,:,:,:,:,:,K) = (uCorez(:,:,:,:,:,:,:,Kp) ...
    -2*uCorez(:,:,:,:,:,:,:,K) +uCorez(:,:,:,:,:,:,:,Km)); % \delta^2
  end% if odd/even
%{
\end{matlab}
Recursively take \(\delta^2\) of these to form higher order
centred differences.
\begin{matlab}
%}
   for k = 3:ordCC    
    dmux(k,:,:,:,:,:,I,:,:) =     dmux(k-2,:,:,:,:,:,Ip,:,:) ...
    -2*dmux(k-2,:,:,:,:,:,I,:,:) +dmux(k-2,:,:,:,:,:,Im,:,:);    
    dmuy(k,:,:,:,:,:,:,J,:) =     dmuy(k-2,:,:,:,:,:,:,Jp,:) ...
    -2*dmuy(k-2,:,:,:,:,:,:,J,:) +dmuy(k-2,:,:,:,:,:,:,Jm,:);
    dmuz(k,:,:,:,:,:,:,:,K) =     dmuz(k-2,:,:,:,:,:,:,:,Kp) ...
    -2*dmuz(k-2,:,:,:,:,:,:,:,K) +dmuz(k-2,:,:,:,:,:,:,:,Km);
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet face values for
each patch \cite[]{Roberts06d, Bunder2013b}, using weights
computed in \verb|configPatches2()|. Here interpolate to
specified order.

Where next-to-face values interpolate to the opposite
face-values. When we have an ensemble of configurations,
different configurations might be coupled to each other, as
specified by \verb|patches.le|, \verb|patches.ri|,
\verb|patches.to|, \verb|patches.bo|,
\verb|patches.fr| and \verb|patches.ba|.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two faces
u(nx,2:(ny-1),2:(nz-1),:,patches.ri,I,:,:) ...
  = uCorex(1,:,:,:,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsr(:,1).*dmux(:,1,:,:,:,:,:,:,:) ),1);  
u(1 ,2:(ny-1),2:(nz-1),:,patches.le,I,:,:) ...
  = uCorex(k,:,:,:,:,:,:,:)*(1-patches.stag) ...      
    +shiftdim(sum( patches.Cwtsl(:,1).*dmux(:,k,:,:,:,:,:,:,:) ),1);
u(2:(nx-1),ny,2:(nz-1),:,patches.to,:,J,:) ...
  = uCorey(:,1,:,:,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsr(:,2).*dmuy(:,:,1,:,:,:,:,:,:) ),1);
u(2:(nx-1),1 ,2:(nz-1),:,patches.bo,:,J,:) ...
  = uCorey(:,k,:,:,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsl(:,2).*dmuy(:,:,k,:,:,:,:,:,:) ),1);
u(2:(nx-1),2:(ny-1),nz,:,patches.fr,:,:,K) ...
  = uCorez(:,:,1,:,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsr(:,3).*dmuz(:,:,:,1,:,:,:,:,:) ),1);
u(2:(nx-1),2:(ny-1),1 ,:,patches.ba,:,:,K) ...
  = uCorez(:,:,k,:,:,:,:,:)*(1-patches.stag) ...
    +shiftdim(sum( patches.Cwtsl(:,3).*dmuz(:,:,:,k,:,:,:,:,:) ),1);
%{
\end{matlab}



\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic. We interpolate in
terms of the patch index~\(j\), say, not directly in space. 
As the macroscale fields are \(N\)-periodic in the patch
index~\(j\), the macroscale Fourier transform writes the
centre-patch values as \(U_j=\sum_{k}C_ke^{ik2\pi j/N}\).
Then the face-patch values \(U_{j\pm r}
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
patch-ratio is effectively halved. The patch faces are near
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
  krx = shiftdim( rx*2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax) ,-4);
  kMay = floor((Ny-1)/2); 
  kry = shiftdim( ry*2*pi/Ny*(mod((0:Ny-1)+kMay,Ny)-kMay) ,-5);
  kMaz = floor((Nz-1)/2); 
  krz = shiftdim( rz*2*pi/Nz*(mod((0:Nz-1)+kMaz,Nz)-kMaz) ,-6);
%{
\end{matlab}
%Test for reality of the field values, and define a function
%accordingly.  Could be problematic if some variables are
%real and some are complex, or if variables are of quite
%different sizes.
%\begin{matlab}
%%}
%  if max(abs(imag(u(:))))<1e-9*max(abs(u(:)))
%       uclean=@(u) real(u);
%  else uclean=@(u) u; 
%  end
%%{
%\end{matlab}
Compute the Fourier transform of the patch centre-values for
all the fields.  Unless doing patch-edgy interpolation when
FT the next-to-face values.  If there are an even number of
points, then if complex, treat as positive wavenumber, but
if real, treat as cosine. When using an ensemble of
configurations, different configurations might be coupled to
each other, as specified by \verb|patches.le|,
\verb|patches.ri|, \verb|patches.to|, \verb|patches.bo|,
\verb|patches.fr| and \verb|patches.ba|.
\begin{matlab}
%}
% indices of interior
ix=(2:nx-1)';  iy=2:ny-1;  iz=shiftdim(2:nz-1,-1); 
if ~patches.EdgyInt
     Cle = fft(fft(fft( u(i0,j0,k0,:,:,:,:,:) ...
         ,[],6),[],7),[],8); 
     Cbo = Cle;  Cba = Cle;
else 
     Cle = fft(fft(fft( u(   2,iy,iz ,:,patches.le,:,:,:) ...
         ,[],6),[],7),[],8);
     Cri = fft(fft(fft( u(nx-1,iy,iz ,:,patches.ri,:,:,:) ...
         ,[],6),[],7),[],8);
     Cbo = fft(fft(fft( u(ix,2   ,iz ,:,patches.bo,:,:,:) ...
         ,[],6),[],7),[],8);
     Cto = fft(fft(fft( u(ix,ny-1,iz ,:,patches.to,:,:,:) ...
         ,[],6),[],7),[],8);
     Cba = fft(fft(fft( u(ix,iy,2    ,:,patches.ba,:,:,:) ...
         ,[],6),[],7),[],8);
     Cfr = fft(fft(fft( u(ix,iy,nz-1 ,:,patches.fr,:,:,:) ...
         ,[],6),[],7),[],8);
end     
%{
\end{matlab}
Fill in the cross of Fourier-shifted mid-values
\begin{matlab}
%}
if ~patches.EdgyInt 
  % yz-fraction of kry/z along left/right faces
  ks = (iy-j0)*2/(ny-1).*kry+(iz-k0)*2/(nz-1).*krz; 
  Cle = Cle.*exp(1i*ks); 
  Cri = Cle;
  % xz-fraction of krx/z along bottom/top faces
  ks = (ix-i0)*2/(nx-1).*krx+(iz-k0)*2/(nz-1).*krz;
  Cbo = Cbo.*exp(1i*ks); 
  Cto = Cbo;
  % xy-fraction of krx/y along back/front faces
  ks = (ix-i0)*2/(nx-1).*krx+(iy-j0)*2/(ny-1).*kry;
  Cba = Cba.*exp(1i*ks); 
  Cfr = Cba;
end
%{
\end{matlab}
Now invert the triple Fourier transforms to complete
interpolation.
\begin{matlab}
%}
u(nx,iy,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cle.*exp(1i*(stagShift+krx))  ,[],6),[],7),[],8) );
u( 1,iy,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cri.*exp(1i*(stagShift-krx))  ,[],6),[],7),[],8) );
u(ix,ny,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cbo.*exp(1i*(stagShift+kry))  ,[],6),[],7),[],8) );
u(ix, 1,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cto.*exp(1i*(stagShift-kry))  ,[],6),[],7),[],8) );
u(ix,iy,nz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cba.*exp(1i*(stagShift+krz))  ,[],6),[],7),[],8) );
u(ix,iy, 1,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cfr.*exp(1i*(stagShift-krz))  ,[],6),[],7),[],8) );
end% if spectral 
%{
\end{matlab}
Nan the values in corners and edges, of every 3D patch.
\begin{matlab}
%}
u(:,[1 ny],[1 nz],:,:,:,:,:)=nan; 
u([1 nx],:,[1 nz],:,:,:,:,:)=nan; 
u([1 nx],[1 ny],:,:,:,:,:,:)=nan; 
%end%spmd %?????????????????
end% function spmdPatchEdgeInt3
%{
\end{matlab}
Fin, returning the 8D array of field values with
interpolated faces. 
\end{devMan} 
%}
