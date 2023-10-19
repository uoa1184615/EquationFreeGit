% patchEdgeIntCore() provides the interpolation across 1D of
% space for multi-D patches of a lattice system such as PDE
% discretisations.  AJR, 19 Oct 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeIntCore()}: sets multi-D patch
face values from 1D macroscale interpolation in one direction}
\label{sec:patchEdgeIntCore}


Couples multi-D patches across multi-D space by computing
their face values via macroscale interpolation.  Assumes
patch face values are determined by macroscale interpolation
of the patch centre-plane values \cite[]{Roberts2011a,
Bunder2019d}, or patch next-to-face values which appears
better \cite[]{Bunder2020a}.  This function is used by
\verb|patchEdgeInt1()|, \verb|patchEdgeInt2()| and 
\verb|patchEdgeInt3()|.

Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global
struct~\verb|patches|.
\begin{matlab}
%}
function u = patchEdgeIntCore(l,u,q,patches,stagShift ...
            ,n1,nl,n2,nEnsem,N1,Nl,N2,lLo,lHi)
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|l| the direction in space of this interpolation.

\item \verb|u| is a vector\slash array with number of elements
$\verb|n1| \cdot (\verb|nl|-2) \cdot \verb|n2| \cdot \verb|nEnsem|
\cdot \verb|N1| \cdot \verb|Nl| \cdot \verb|N2|$.
 
\item \verb|q| spatial coordinates in the direction with
number of elements $\verb|nl| \cdot \verb|Nl|$.

\item \verb|patches| a struct set by \verb|configPatches.()|.

\item Get \verb|stagShift| from \verb|patches|??

\item \ldots

\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 7D array, $\verb|n1| \cdot
\verb|nl| \cdot \verb|n2| \cdot
\verb|nEnsem| \cdot \verb|N1| \cdot \verb|Nl|
\cdot \verb|N2|$, of the fields with face values set by
interpolation.
\end{itemize}








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

Reshape the field values accordingly, fattened for the new
edge values, reshape the coordinate values, and set
macroscale indices and periodic neighbours.
\begin{matlab}
%}
u = cat(2,nan(n1,  1 ,n2,nEnsem,N1,Nl,N2) ...
   ,reshape(u,n1,nl-2,n2,nEnsem,N1,Nl,N2) ...
         ,nan(n1,  1 ,n2,nEnsem,N1,Nl,N2) );
q = reshape(q,1,nl,1,1,1,Nl,1);
Il=1:Nl; Ip=mod(Il,Nl)+1; Im=mod(Il-2,Nl)+1;
l0 = round((nl+1)/2);
%{
\end{matlab}



\subsection{Periodic macroscale interpolation schemes}
First get the size ratios of the patches in the direction.
\begin{matlab}
%}
if patches.periodic(l)
rl = patches.ratio(l);
%{
\end{matlab}

\subsubsection{Lagrange interpolation gives patch-face values}
Compute centred differences of the mid-patch values for the
macro-interpolation, of all fields.  Here the domain is
macro-periodic.
\begin{matlab}
%}
ordCC = patches.ordCC;
if ordCC>0 % then finite-width polynomial interpolation
%{
\end{matlab}


\paragraph{\(l\)-normal face values} The patch-edge values
are either interpolated from the next-to-edge-face values,
or from the centre-cross-plane values (not the patch-centre
value itself as that seems to have worse properties in
general).  Have not yet implemented core averages.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u(:,[2 nl-1],:,:,:,:,:);
  else % interpolate centre-cross values
    U = u(:,   l0   ,:,:,:,:,:);
  end;%if patches.EdgyInt
%{
\end{matlab}
Just in case any last array dimension(s) are one, we force a
padding of the sizes, then adjoin the extra dimension for
the subseluent array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,7-length(szUO)) ordCC];
%{
\end{matlab}
Use finite difference formulas for the interpolation, so
store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in these arrays.  When parallel, in order
to preserve the distributed array structure we use an index
at the end for the differences.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 8D
  else   dmu = zeros(szUO,patches.codist); % 8D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  %disp('starting standard interpolation')  
  dmu(:,:,:,:,:,Il,:,1) = (U(:,:,:,:,:,Ip,:) ...
                          -U(:,:,:,:,:,Im,:))/2; %\mu\delta 
  dmu(:,:,:,:,:,Il,:,2) = (U(:,:,:,:,:,Ip,:) ...
     -2*U(:,:,:,:,:,Il,:) +U(:,:,:,:,:,Im,:));   %\delta^2    
%{
\end{matlab}
Recursively take $\delta^2$ of these to form successively
higher order centred differences in space.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,:,Il,:,k) =     dmu(:,:,:,:,:,Ip,:,k-2) ...
    -2*dmu(:,:,:,:,:,Il,:,k-2) +dmu(:,:,:,:,:,Im,:,k-2);    
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet face values for
each patch \cite[]{Roberts06d, Bunder2013b}, using the weights
pre-computed by \verb|configPatches3()|. Here interpolate to
specified order.

For the case where next-to-face values interpolate to the
opposite face-values: when we have an ensemble of
configurations, different configurations might be coupled to
each other, as specified in \verb|lLo| and~\verb|lHi| (from
\verb|patches.le|, \verb|patches.ri|, \verb|patches.bo|,
\verb|patches.to|, \verb|patches.ba|, and
\verb|patches.fr|).
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two faces
u(:,nl,:,lHi,:,:,:) ...
  = U(:,1,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,l),-7).*dmu(:,1,:,:,:,:,:,:) ,8);  
u(:,1 ,:,lLo,:,:,:) ...
  = U(:,k,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,l),-7).*dmu(:,k,:,:,:,:,:,:) ,8);
%{
\end{matlab}




\subsubsection{Case of spectral interpolation}
Assumes the domain is macro-periodic.
\begin{matlab}
%}
else% patches.ordCC<=0, spectral interpolation
%{
\end{matlab}
We interpolate in terms of the patch index, $I$~say, not
directly in space. As the macroscale fields are $N$-periodic
in the patch index~$I$, the macroscale Fourier transform
writes the centre-patch values as $U_I=\sum_{k}C_ke^{ik2\pi
I/N}$. Then the face-patch values $U_{I\pm r}
=\sum_{k}C_ke^{ik2\pi/N(I\pm r)} =\sum_{k}C'_ke^{ik2\pi
I/N}$ where $C'_k =C_ke^{ikr2\pi/N}$. For $N$~patches we
resolve `wavenumbers' $|k|<N/2$, so set row vector
$\verb|ks| =k2\pi/N$ for `wavenumbers' $\mathcode`\,="213B
k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)$ for
odd~$N$, and $\mathcode`\,="213B k=(0,1, \ldots, k_{\max},
\pm(k_{\max}+1) -k_{\max}, \ldots, -1)$ for even~$N$.


\paragraph{\(l\)-normal face values} Now set wavenumbers
into a vector at the correct dimension.  In the case of
even~$N$ these compute the $+$-case for the highest
wavenumber zig-zag mode, $\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, +(k_{\max}+1) -k_{\max}, \ldots, -1)$.
\begin{matlab}
%}
  kMax = floor((Nl-1)/2); 
  kr = shiftdim( rl*2*pi/Nl*(mod((0:Nl-1)+kMax,Nl)-kMax) ,-4);
%{
\end{matlab}
Compute the Fourier transform of the patch values on the
centre-planes for all the fields.  Unless doing patch-edgy
interpolation when FT the next-to-face values.  If there are
an even number of points, then if complex, treat as positive
wavenumber, but if real, treat as cosine. When using an
ensemble of configurations, different configurations might
be coupled to each other, as specified by \verb|lLo|
and~\verb|lHi| (from \verb|patches.le|, \verb|patches.ri|,
\verb|patches.to|, \verb|patches.bo|, \verb|patches.fr| and
\verb|patches.ba|).
\begin{matlab}
%}
if ~patches.EdgyInt
     Cm = fft( u(:,l0,:,:,:,:,:) ,[],6); 
     Cp = Cm;
else 
     Cm = fft( u(:,   2,:,lLo,:,:,:) ,[],6);
     Cp = fft( u(:,nl-1,:,lHi,:,:,:) ,[],6);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Now invert the Fourier transforms to complete interpolation.
Enforce reality when appropriate with \verb|uclean()|. 
\begin{matlab}
%}
u(:,nl,:,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],6) );
u(:, 1,:,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],6) );
%{
\end{matlab}

End the two periodic cases of spectral or finite-difference
interpolation.
\begin{matlab}
%}
end% if ordCC>0  
%{
\end{matlab}



\subsection{Non-periodic macroscale interpolation}
\begin{matlab}
%}
else% patches.periodic false
assert(~patches.stag, ...
'not implementing staggered grids for non-periodic')
%{
\end{matlab}
Determine the order of interpolation~\verb|pl| (potentially
different in the different directions!), and hence size of
the (forward) divided difference tables in~\verb|F|~(8D) for
interpolating to left/right, top/bottom, and front/back
faces. Because of the product-form of the patch grid, and
because we are doing \emph{only} either edgy interpolation
or cross-patch interpolation (\emph{not} just the centre
patch value), the interpolations are all 1D interpolations.
\begin{matlab}
%}
if patches.ordCC<1,  pl = Nl-1;
  else pl = min(patches.ordCC,Nl-1); 
  end
%{
\end{matlab}

\subsubsection{\(l\)-direction values}
Set function values in second column of the tables for every
variable and across ensemble.  For~\verb|EdgyInt|, the
`reversal' of the next-to-face values are because their
values are to interpolate to the opposite face of each
patch. \todo{Have no plans to implement core averaging as
yet.}
\begin{matlab}
%}
  F = nan(n1,patches.EdgyInt+1,n2,nEnsem,N1,Nl,N2,pl+1);
  if patches.EdgyInt % interpolate next-to-face values
    F(:,:,:,:,:,:,:,1) = u(:,[nl-1 2],:,:,:,:,:);
    Q = q(:,[nl-1 2],:,:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,:,1) = u(:,   l0   ,:,:,:,:,:);
    Q = q(:,   l0   ,:,:,:,:,:);
  end%if patches.EdgyInt
%{
\end{matlab}

\paragraph{Form tables of divided differences} 
Compute tables of (forward) divided differences
\cite[e.g.,][]{DividedDifferences} for every variable, and
across ensemble, and in both directions, and for all three
types of faces (left/right, top/bottom, and front/back).
Recursively find all divided differences in the respective
direction.
\begin{matlab}
%}
for m = 1:pl
  i = 1:Nl-m;
  F(:,:,:,:,:,i,:,m+1) ...
  = ( F(:,:,:,:,:,i+1,:,m)-F(:,:,:,:,:,i,:,m)) ...
   ./(Q(:,:,:,:,:,i+m,:)  -Q(:,:,:,:,:,i,:));
end
%{
\end{matlab}

\paragraph{Interpolate with divided differences} Now
interpolate to find the face-values on left/right faces
at~\verb|Qface| for every interior in other coordinates.
\begin{matlab}
%}
Qface = q(:,[1 nl],:,:,:,:,:);
%{
\end{matlab}
Code Horner's recursive evaluation of the interpolation
polynomials.  Indices~\verb|i| are those of the left face of
each interpolation stencil, because the table is of forward
differences.  This alternative: the case of order~\(p_x\),
\(p_y\) and~\(p_z\) interpolation across the domain,
asymmetric near the boundaries of the rectangular domain.
\begin{matlab}
%}
  i = max(1,min(1:Nl,Nl-ceil(pl/2))-floor(pl/2));
  Uface = F(:,:,:,:,:,i,:,pl+1);
  for m = pl:-1:1
    Uface = F(:,:,:,:,:,i,:,m) ...
    +(Qface-Q(:,:,:,:,:,i+m-1,:)).*Uface;
  end
%{
\end{matlab}

Finally, insert face values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(:,1 ,:,lLo,:,:,:) = Uface(:,1,:,:,:,:,:);
u(:,nl,:,lHi,:,:,:) = Uface(:,2,:,:,:,:,:);
%{
\end{matlab}



\subsubsection{Optional NaNs for safety}
We want a user to set outer face values on the extreme
patches according to the microscale boundary conditions that
hold at the extremes of the domain.  Consequently, unless
testing, override their computed interpolation values
with~\verb|NaN|.
\begin{matlab}
%}
if isfield(patches,'intTest')&&patches.intTest
else %disp('usual case requires user to set bdry values')
    u(:, 1,:,:,:, 1,:) = nan;
    u(:,nl,:,:,:,Nl,:) = nan;
end%if
%{
\end{matlab}

End of the non-periodic interpolation code.
\begin{matlab}
%}
end%if patches.periodic else
%{
\end{matlab}


Fin, returning the 7D array of field values with
interpolated face. 
\begin{matlab}
%}
end% function patchEdgeIntCore
%{
\end{matlab}
%}
