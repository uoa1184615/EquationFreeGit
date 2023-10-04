function h = slices(X,Y,Z,U,ix,iy,iz,fn)
% Analogous to slice() but better for contorted data grids
% as no interpolation and allows NaNs in the data grid. 
% If some elements of X,Y,Z are ever NaN, then those points
% are not plotted. (It is not enough to NaN elements of U.)
% AJR, 4 Oct 2023

% optional eighth argument is either @surf (default), or
% @mesh, or @contour4
if nargin<8, fn=@surf; end

% cater for extra argument to contour4
if func2str(fn)=="contour4"
  cU = linspace(min(U(:)),max(U(:)),21);
  fn = @(x,y,z,u) contour4(x,y,z,u,cU);
end%if

% Cater for grids to be input as 1D arrays of any 1D shape
if sum(size(X)>1)>1, nX=size(X,1);
else nX=length(X); X=reshape(X,[],1,1); end
if sum(size(Y)>1)>1, nY=size(Y,2);
else nY=length(Y); Y=reshape(Y,1,[],1); end
if sum(size(Z)>1)>1, nZ=size(Z,3);
else nZ=length(Z); Z=reshape(Z,1,1,[]); end
% expand any to requisite full shape
X=X+zeros(nX,nY,nZ);
Y=Y+zeros(nX,nY,nZ);
Z=Z+zeros(nX,nY,nZ);

% x surfaces
for i=ix(:)'
  fn(squeeze(X(i,:,:)),squeeze(Y(i,:,:)),squeeze(Z(i,:,:)) ...
    ,squeeze(U(i,:,:)) )
  hold on
end
% y surfaces
for i=iy(:)'
  fn(squeeze(X(:,i,:)),squeeze(Y(:,i,:)),squeeze(Z(:,i,:)) ...
    ,squeeze(U(:,i,:)) )
  hold on
end
% z surfaces
for i=iz(:)'
  fn(squeeze(X(:,:,i)),squeeze(Y(:,:,i)),squeeze(Z(:,:,i)) ...
    ,squeeze(U(:,:,i)) )
  hold on
end
hold off

end%function slices

