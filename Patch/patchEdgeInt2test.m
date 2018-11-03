%{
Testing is passed by all this randomness. Nov 2018
%}
clear all, close all
global patches
for realisation=1:99
Lx=1+3*rand, Ly=1+3*rand
nSubP=1+2*randi(3,1,2)
ratios=rand(1,2)/2
nPatch=3+randi(4,1,2)
%nPatch=[4 6], nSubP=3
configPatches2(@sin,[0 Lx 0 Ly],nan,nPatch,0,ratios,nSubP)
kx=randi([0 ceil((nPatch(1)-1)/2)])
ky=randi([0 ceil((nPatch(2)-1)/2)])
phix=pi*rand*(2*kx~=nPatch(1))
phiy=pi*rand*(2*ky~=nPatch(2))
% generate 2D array via auto-replication
u0=sin(2*pi*kx*patches.x(:)/Lx+phix).*sin(2*pi*ky*patches.y(:)'/Ly+phiy);
% reshape into 4D array
[nx,Nx]=size(patches.x);
[ny,Ny]=size(patches.y);
u0=reshape(u0,[nx Nx ny Ny]);
u0=permute(u0,[1 3 2 4]);
% copy and NaN the edges
u=u0; u([1 end],:,:,:)=nan; u(:,[1 end],:,:)=nan;
u=patchEdgeInt2(u);
err=u-u0;
normerr=norm(err(:))
if normerr>1e-12, error('2D interpolation failed'), end
end
