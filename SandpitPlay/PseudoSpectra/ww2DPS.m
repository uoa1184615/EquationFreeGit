%{ 
compute pseudospectra of ww2D to for damped linear ideal
wave system on staggered patches of staggered micro-grid in
2D. Minimal case takes a minute.  Pseudospectrum seems to
show some sub-patch fast-waves are a little sensitive in
that there is a 'spike' of the zero-contour to real-part
just positive. AJR, May 2023
%}
clear all
global p N r dX micro cDrag cVisc ordCoup
%{
Linear drag on velocity: find all but one zero eigenvalue
moves to be -cDrag; and that all the imaginary eigenvalues
have real part of -cDrag/2.
%}
cDrag=0
%{
Viscous drag in momentum PDEs: 
%}
cVisc=1e-2

% p = patch micro-grid paramater, size is (4p+1)^2 p=1,2,...
p=2
% N = even number of patches in each dirn of periodic domain
% for spectral interpolation need N/2 to be odd
N=6
% r = fraction of distance a patch inner-edge is to centre of next patch
r=0.1
% ordCoup = 0 for spectral and 1 for low-order interpatch coupling
ordCoup=0
% dX = space step of the inter-patch distance of patch centres
dX=2*pi/N
X=(0.5:N)*dX; Y=X; % mid-point patches
% micro = index into locations of the PDEs on the micro-grid
% basename for graphics files
basename='ww2DXeig';
if ordCoup, basename=[basename 'Lo']; end

% dx = space step of the micro-grid
dx=r*dX/(2*p-1)
% locations on axes of micro-grid points in all patches 
[x,y]=ndgrid(dx*(-2*p:2*p),X);
x=x+y; y=x;
% Maybe use double letters for 2D grids
% mesh of all locations, used for ICs
[xx,XX,yy,YY]=ndgrid(dx*(-2*p:2*p),X,dx*(-2*p:2*p),X);
xx=xx+XX; yy=yy+YY;

% micro-grid indexes inside a pxp patch where
n=4*p+1 % patch microgrid is nxn
iMid=2*p+1; % mid-point index
iie=3:2:n-2; % interior even points (relative to centre-patch)
iio=4:2:n-3; % interior odd points (empty if p=1)
io=2:2:n-1; % odd points including edges.
ie=1:2:n; % even points including edges.

% index of patch types: IH,IH are H-patches; IQ,IH are
% U-patches; IH,IQ are V-patches; and IQ,IQ are empty patches
IH=1:2:N; IQ=2:2:N;

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
u=nan(n,N,n,N);
% Must set some values to find pattern of data in the 4D array
% Here set to equilibrium.
h0=@(x,y) 0; %sin(x+y);
u0=@(x,y) 0; %sin(x);
v0=@(x,y) 0; %-cos(y);
% say h (e,H,e,H) (e,H,o,Q) (o,Q,e,H)
u(iie,IH,iie,IH)=h0(xx(iie,IH,iie,IH),yy(iie,IH,iie,IH)); 
u(iie,IH,iio,IQ)=h0(xx(iie,IH,iio,IQ),yy(iie,IH,iio,IQ));
u(iio,IQ,iie,IH)=h0(xx(iio,IQ,iie,IH),yy(iio,IQ,iie,IH));
% u field (o,H,e,H) (o,H,o,Q) (e,Q,e,H)
u(iio,IH,iie,IH)=u0(xx(iio,IH,iie,IH),yy(iio,IH,iie,IH)); 
u(iio,IH,iio,IQ)=u0(xx(iio,IH,iio,IQ),yy(iio,IH,iio,IQ));
u(iie,IQ,iie,IH)=u0(xx(iie,IQ,iie,IH),yy(iie,IQ,iie,IH));
% v field (e,H,o,H) (e,H,e,Q) (o,Q,o,H)
u(iie,IH,iio,IH)=v0(xx(iie,IH,iio,IH),yy(iie,IH,iio,IH)); 
u(iie,IH,iie,IQ)=v0(xx(iie,IH,iie,IQ),yy(iie,IH,iie,IQ));
u(iio,IQ,iio,IH)=v0(xx(iio,IQ,iio,IH),yy(iio,IQ,iio,IH));
% finally get the gather/scatter indices of micro-grid PDE points 
micro=find(~isnan(u));
nMicro=length(micro)
huv=u(micro); 


%% Compute a time derivative as quasi-check
ut=nan(n,N,n,N);
ut(micro)=ww2DXdudt(0,huv);


%% Compute eigenvalues then plot in complex plane
small=1; % perturb the surface by this much, here linear problem
Jac=nan(nMicro);
for j=1:nMicro
    Jac(:,j)=ww2DXdudt(j,huv+small*((1:nMicro)==j)')/small;
end
figure(1)
figSc = inf % <- set figure scaling
multiscalePseudoSpectra(Jac,figSc,1,200)
figName = [mfilename num2str(round(abs(log10(figSc))))]
set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 14] ...
      ,'renderer','Painters')
print('-depsc2',figName)
matlab2tikz([figName '.tex'],'showInfo',false ...
    ,'noSize',true,'parseStrings',false,'showWarnings',false ...
    ,'extraCode',['\tikzsetnextfilename{' figName '}'] ...
    ,'extraAxisOptions','\extraAxisOptions' ...
    ,'checkForUpdates',false)




function dhuvdt=ww2DXdudt(t,huv)
% Compute time derivatives of shallow water wave system on
% staggered patches of staggered micro-grid, with extra
% outer-edge perimeter. 
% AJR, 1 Mar 2019
% I/O: t = time, ignored as autonomous
% huv = h,u,v values on staggered patches of staggered micro-grid
% dhuvdt = array of time derivatives of huv values
% Parameters (coupling has other global params)
% micro = index into locations of the PDEs on the micro-grid
% cDrag = linear drag on velocity
% cVisc = linear viscosity coeff
global micro cDrag cVisc ordCoup

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
% Specify the function that couples patches by
% providing patch edge values through  whatever
% mechanism
if ordCoup, [ut,IH,IQ,iio,iie,dx]=stag2DXcouple(huv);
else [ut,IH,IQ,iio,iie,dx]=stag2DXspectral(huv);
end
io=[2 iio iio(end)+2];
ie=[1 iie iie(end)+2];
% scatter data fields into separate arrays for safety
h=nan(size(ut)); u=h; v=h;
h(ie,IH,ie,IH)=ut(ie,IH,ie,IH);
h(ie,IH,io,IQ)=ut(ie,IH,io,IQ);
h(io,IQ,ie,IH)=ut(io,IQ,ie,IH);
u(io,IH,ie,IH)=ut(io,IH,ie,IH); 
u(io,IH,io,IQ)=ut(io,IH,io,IQ);
u(ie,IQ,ie,IH)=ut(ie,IQ,ie,IH);
v(ie,IH,io,IH)=ut(ie,IH,io,IH); 
v(ie,IH,ie,IQ)=ut(ie,IH,ie,IQ);
v(io,IQ,io,IH)=ut(io,IQ,io,IH);

%% Compute time derivatives in 3 different patch types
ut=nan(size(ut));
% dx = space step of the micro-grid
dx2=2*dx; dy2=dx2; dx2sq=dx2^2; 

% h_t=-u_x-v_y at (e,H,e,H) (e,H,o,Q) (o,Q,e,H)
ut(iie,IH,iie,IH)=-(u(iie+1,IH,iie,IH)-u(iie-1,IH,iie,IH))/dx2 ...
                  -(v(iie,IH,iie+1,IH)-v(iie,IH,iie-1,IH))/dy2 ; 
ut(iie,IH,iio,IQ)=-(u(iie+1,IH,iio,IQ)-u(iie-1,IH,iio,IQ))/dx2 ...
                  -(v(iie,IH,iio+1,IQ)-v(iie,IH,iio-1,IQ))/dy2 ;
ut(iio,IQ,iie,IH)=-(u(iio+1,IQ,iie,IH)-u(iio-1,IQ,iie,IH))/dx2 ...
                  -(v(iio,IQ,iie+1,IH)-v(iio,IQ,iie-1,IH))/dy2 ;
                  
% u_t=-h_x-?u+?del2u  at (o,H,e,H) (o,H,o,Q) (e,Q,e,H)
dels=(u(iio-2,IH,iie,IH)-2*u(iio,IH,iie,IH)+u(iio+2,IH,iie,IH))/dx2sq ...
    +(u(iio,IH,iie-2,IH)-2*u(iio,IH,iie,IH)+u(iio,IH,iie+2,IH))/dx2sq ;
ut(iio,IH,iie,IH)=-(h(iio+1,IH,iie,IH)-h(iio-1,IH,iie,IH))/dx2 ...
    -cDrag*u(iio,IH,iie,IH) ...
    +cVisc*dels; 
dels=(u(iio-2,IH,iio,IQ)-2*u(iio,IH,iio,IQ)+u(iio+2,IH,iio,IQ))/dx2sq ...
    +(u(iio,IH,iio-2,IQ)-2*u(iio,IH,iio,IQ)+u(iio,IH,iio+2,IQ))/dx2sq ;
ut(iio,IH,iio,IQ)=-(h(iio+1,IH,iio,IQ)-h(iio-1,IH,iio,IQ))/dx2 ...
    -cDrag*u(iio,IH,iio,IQ) ...
    +cVisc*dels ;
dels=(u(iie-2,IQ,iie,IH)-2*u(iie,IQ,iie,IH)+u(iie+2,IQ,iie,IH))/dx2sq ...
    +(u(iie,IQ,iie-2,IH)-2*u(iie,IQ,iie,IH)+u(iie,IQ,iie+2,IH))/dx2sq ;
ut(iie,IQ,iie,IH)=-(h(iie+1,IQ,iie,IH)-h(iie-1,IQ,iie,IH))/dx2 ...
    -cDrag*u(iie,IQ,iie,IH) ...
    +cVisc*dels ;
    
% v_t=-h_y-?v+?del2v  at (e,H,o,H) (e,H,e,Q) (o,Q,o,H)
dels=(v(iie-2,IH,iio,IH)-2*v(iie,IH,iio,IH)+v(iie+2,IH,iio,IH))/dx2sq ...
    +(v(iie,IH,iio-2,IH)-2*v(iie,IH,iio,IH)+v(iie,IH,iio+2,IH))/dx2sq ;
ut(iie,IH,iio,IH)=-(h(iie,IH,iio+1,IH)-h(iie,IH,iio-1,IH))/dy2 ...
    -cDrag*v(iie,IH,iio,IH)  ...
    +cVisc*dels ; 
dels=(v(iie-2,IH,iie,IQ)-2*v(iie,IH,iie,IQ)+v(iie+2,IH,iie,IQ))/dx2sq ...
    +(v(iie,IH,iie-2,IQ)-2*v(iie,IH,iie,IQ)+v(iie,IH,iie+2,IQ))/dx2sq ;
ut(iie,IH,iie,IQ)=-(h(iie,IH,iie+1,IQ)-h(iie,IH,iie-1,IQ))/dy2 ...
    -cDrag*v(iie,IH,iie,IQ)  ...
    +cVisc*dels ;
dels=(v(iio-2,IQ,iio,IH)-2*v(iio,IQ,iio,IH)+v(iio+2,IQ,iio,IH))/dx2sq ...
    +(v(iio,IQ,iio-2,IH)-2*v(iio,IQ,iio,IH)+v(iio,IQ,iio+2,IH))/dx2sq ;
ut(iio,IQ,iio,IH)=-(h(iio,IQ,iio+1,IH)-h(iio,IQ,iio-1,IH))/dy2 ...
    -cDrag*v(iio,IQ,iio,IH)  ...
    +cVisc*dels ;
    
% extract and reshape
dhuvdt=ut(micro);
end





function [u,IH,IQ,iio,iie,dx]=stag2DXspectral(huv)
% Compute patch-edge values on staggered patches of
% staggered micro-grid when patches have edges that are two
% grid points thick.  AJR, 1 Mar 2019
% I/O: huv = vector of h,u,v values on staggered
%     patches of staggered micro-grid
% u = 4D array of huv values, with patch-edge values interpolated
% IH = H-patch indices
% IQ = U,V-patch indices, in conjunction with IH
% iio = indices of internal lattice, odd-points
% iie = indices of internal lattice, even-points
% dx = microscale lattice spacing
% Parameters
% p = patch micro-grid paramater, size is (4p+1)^2 p=1,2,...
% N = even number of patches in each dirn of periodic domain
% r = fraction of distance a patch inner-edge is to centre of next patch
% dX = space step of the inter-patch distance of patch centres
% micro = index into locations of the PDEs on the micro-grid
global p N r dX micro 

% dx = space step of the micro-grid
dx=r*dX/(2*p-1);
ro=2*p*dx/dX;  % = fraction of distance to outer-edge

% micro-grid indexes inside a nxn patch where
n=4*p+1; % patch microgrid is nxn
iMid=2*p+1; % mid-point index
iie=3:2:n-2; % interior even points (relative to the mid-patch)
iio=4:2:n-3; % interior odd points (empty when p=1)

% index of patch types: IH,IH are H-patches; IQ,IH are
% U-patches; IH,IQ are V-patches; and IQ,IQ are empty patches
IH=1:2:N; IQ=2:2:N;

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
u=nan(n,N,n,N);
u(micro)=huv;

%% Interpolate spectrally to set necessary patch edge values
% Omits the corner points on all patches
% N/2, 2D Fourier coeffs of the three types of patches
cHk=fft2(squeeze(u(iMid,IH,iMid,IH))); 
cUk=fft2(squeeze(u(iMid,IQ,iMid,IH))); 
cVk=fft2(squeeze(u(iMid,IH,iMid,IQ))); 
% do we need some dXs in here? no, they cancel
kMax=(N/2-1)/2; %max wavenumber for interpolation, N/2 must be odd
% shift one patch, one macro-grid, to the right/up
k1=pi/(N/2)*(rem(kMax+(0:N/2-1),N/2)-kMax);
kr = r*k1;    kro = ro*k1;    % shift patch-half-width to the right/up
krp=(1+r)*k1; krop=(1+ro)*k1; % shift one + patch-half-width to the right/up
krm=(1-r)*k1; krom=(1-ro)*k1; % shift one - patch-half-width to the right/up

% First compute all the inner-edge gridpoints
% a dash means shift in x-dirn, undashed is in y-dirn
% a plus = shift right/up, a minus = shift left/down --- in a cell
for j=iie
  ks=(j-iMid)*dx/dX*k1;
  % interpolate U and V values to H-patches
  u(n-1,IH,j,IH)=ifft2(cUk.*exp(1i*( -krm'+ks )));
  u(2  ,IH,j,IH)=ifft2(cUk.*exp(1i*( -krp'+ks )));
  u(j,IH,n-1,IH)=ifft2(cVk.*exp(1i*( -krm +ks')));
  u(j,IH,2  ,IH)=ifft2(cVk.*exp(1i*( -krp +ks')));
  % interpolate H to U-patches
  u(n-1,IQ,j,IH)=ifft2(cHk.*exp(1i*( +krp'+ks )));
  u(2  ,IQ,j,IH)=ifft2(cHk.*exp(1i*( +krm'+ks )));
  % interpolate H to V-patches
  u(j,IH,n-1,IQ)=ifft2(cHk.*exp(1i*( +krp +ks')));
  u(j,IH,2  ,IQ)=ifft2(cHk.*exp(1i*( +krm +ks')));
end
for j=2:2:n-1
  ks=(j-iMid)*dx/dX*k1;
  % interpolate V to U-patches
  u(j,IQ,n-1,IH)=ifft2(cVk.*exp(1i*( k1'+ks'-krm )));
  u(j,IQ,2  ,IH)=ifft2(cVk.*exp(1i*( k1'+ks'-krp )));
  % interpolate U to V-patches
  u(n-1,IH,j,IQ)=ifft2(cUk.*exp(1i*( k1 +ks -krm')));
  u(2  ,IH,j,IQ)=ifft2(cUk.*exp(1i*( k1 +ks -krp')));
  % also do top/bot u of V-patches
  u(j,IH,n-1,IQ)=ifft2(cUk.*exp(1i*( krp +ks' -k1')));
  u(j,IH,2  ,IQ)=ifft2(cUk.*exp(1i*( krm +ks' -k1')));
  % and left/right v of U-patches
  u(n-1,IQ,j,IH)=ifft2(cVk.*exp(1i*( krp' +ks -k1 )));
  u(2  ,IQ,j,IH)=ifft2(cVk.*exp(1i*( krm' +ks -k1 )));
end

% Second compute all the outer-edge gridpoints
% a dash means shift in x-dirn, undashed is in y-dirn
% a plus = shift right/up, a minus = shift left/down --- in a cell
for j=2:2:n-1
  ks=(j-iMid)*dx/dX*k1;
  % interpolate U and V values to H-patches
  u(n,IH,j,IH)=ifft2(cVk.*exp(1i*( +kro'-k1+ks )));
  u(1,IH,j,IH)=ifft2(cVk.*exp(1i*( -kro'-k1+ks )));
  u(j,IH,n,IH)=ifft2(cUk.*exp(1i*( +kro-k1'+ks')));
  u(j,IH,1,IH)=ifft2(cUk.*exp(1i*( -kro-k1'+ks')));
  % interpolate H to U-patches
  u(j,IQ,n,IH)=ifft2(cHk.*exp(1i*( +k1'+ks'+kro )));
  u(j,IQ,1,IH)=ifft2(cHk.*exp(1i*( +k1'+ks'-kro )));
  % third interpolate H to V-patches
  u(n,IH,j,IQ)=ifft2(cHk.*exp(1i*( +k1+ks+kro' )));
  u(1,IH,j,IQ)=ifft2(cHk.*exp(1i*( +k1+ks-kro' )));
end
for j=1:2:n-2
  ks=(j-iMid)*dx/dX*k1;
  % interpolate top/bot H to H-patches
  u(    j,IH,n,IH)=ifft2(cHk.*exp(1i*( +ks'+kro )));
  u(n+1-j,IH,1,IH)=ifft2(cHk.*exp(1i*( -ks'-kro )));
  % and left/right H of H-patches
  u(n,IH,n+1-j,IH)=ifft2(cHk.*exp(1i*( +kro'-ks )));
  u(1,IH,    j,IH)=ifft2(cHk.*exp(1i*( -kro'+ks )));
  % interpolate top/bot U to U-patches
  u(    j,IQ,n,IH)=ifft2(cUk.*exp(1i*( +ks'+kro )));
  u(n+1-j,IQ,1,IH)=ifft2(cUk.*exp(1i*( -ks'-kro )));
  % and left/right U of U-patches
  u(n,IQ,n+1-j,IH)=ifft2(cUk.*exp(1i*( +kro'-ks )));
  u(1,IQ,    j,IH)=ifft2(cUk.*exp(1i*( -kro'+ks )));
  % interpolate top/bot V to V-patches
  u(    j,IH,n,IQ)=ifft2(cVk.*exp(1i*( +ks'+kro )));
  u(n+1-j,IH,1,IQ)=ifft2(cVk.*exp(1i*( -ks'-kro )));
  % and left/right V of V-patches
  u(n,IH,n+1-j,IQ)=ifft2(cVk.*exp(1i*( +kro'-ks )));
  u(1,IH,    j,IQ)=ifft2(cVk.*exp(1i*( -kro'+ks )));
end

end





