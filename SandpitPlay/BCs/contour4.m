function hand = contour4(X,Y,Z,C,nC)
% Draws contours of a quantity C on a surface in 3D-space
% from data sampled on some deformed rectangular grid giving
% the surface shape: X,Y,Z,C must all be the same sized 2D
% arrays.  Uses the current colormap for contour colours. 
% Optional nC, default=11, is the number of contour levels,
% or if vector then specifies the specific contour levels.
% The X,Y,Z arrays should be set sensibly via ndgrid (or
% equivalent), not meshgrid.   AJR, 5 Oct 2023

% if no arguments then draw an example
if nargin==0
  [X,Y]=ndgrid(-0.8:0.1:1,-1:0.2:1);
  Z=0.5*(X.^2-Y.^2);
  clf, colormap(0.8*jet)
  contour4(X,Y,Z, 0.5*exp(-X.^2-Y.^2));
  axis equal, colorbar
  xlabel('x'), ylabel('y'), zlabel('z')
  return
end%if nargin

[ni,nj]=size(C);
assert(all(size(X)==[ni nj] & size(Y)==[ni nj] & size(Z)==[ni nj]) ...
     ,'X,Y,Z,C must all be the same size of 2D array')
i=1:ni-1; j=1:nj-1;

% decide on contours and the colour mapping
if nargin<5, nC=11; end
if length(nC)==1
  cRng = quantile(C(:),[.1 .9]);
  cRng = cRng+[-1 1]*diff(cRng)/8;
  cs = linspace(cRng(1),cRng(2),nC);
  cLim = cRng+[-.5 .5]*diff(cs(1:2));
else%length(nC)>1
  cs=sort(nC(:)');
  cLim = cs([1 end])+[-.5 .5]*(diff(cs([1 end]))/length(cs)+1e-6);
end%length(C)
theCols = colormap; % use the current colormap

for c = cs
  xyzs=[];
  ccol = theCols( ceil(size(theCols,1)*(c-cLim(1))/diff(cLim)) ,:);
  % find up-down neighbours on either side of contour
  [ki,kj] = find( sign(C(i,j)-c)~=sign(C(i,j+1)-c) );
  k = sub2ind([ni,nj],ki,kj);
  % find those with left-right on either side of c
  l = k( sign(C(k)-c)~=sign(C(k+1)-c) );
  xyzs=[xyzs; csegment(l,l+1,l+ni)];
  % find those with diag on either side of c
  l = k( sign(C(k+ni)-c)~=sign(C(k+1)-c) );
  xyzs=[xyzs; csegment(l+ni,l,l+1)];
  % find left-right and diag neighbours on either side of contour
  [ki,kj] = find( (sign(C(i+1,j)-c)~=sign(C(i,j)-c)) ...
                & (sign(C(i+1,j)-c)~=sign(C(i,j+1)-c)) );
  l = sub2ind([ni,nj],ki,kj);
  xyzs=[xyzs; csegment(l+1,l,l+ni)];
  % find left-right and up-down neighbours on either side of contour
  [ki,kj] = find( (sign(C(i+1,j+1)-c)~=sign(C(i+1,j)-c)) ...
                & (sign(C(i+1,j+1)-c)~=sign(C(i,j+1)-c)) );
  l = sub2ind([ni,nj],ki+1,kj+1);
  xyzs=[xyzs; csegment(l,l-1,l-ni)];
  % find up-down and diag neighbours on either side of contour
  [ki,kj] = find( (sign(C(i+1,j)-c)~=sign(C(i+1,j+1)-c)) ...
                & (sign(C(i+1,j)-c)~=sign(C(i,j+1)-c)) );
  l = sub2ind([ni,nj],ki+1,kj+1);
  xyzs=[xyzs; csegment(l-ni,l,l-1)];
  % find left-right and diag neighbours on either side of contour
  [ki,kj] = find( (sign(C(i,j+1)-c)~=sign(C(i+1,j+1)-c)) ...
                & (sign(C(i,j+1)-c)~=sign(C(i+1,j)-c)) );
  l = sub2ind([ni,nj],ki+1,kj+1);
  xyzs=[xyzs; csegment(l-1,l,l-ni)];
  % plot all line segments found for this c
  if ~isempty(xyzs)
    nans=nan+xyzs(:,1);
    xs = [xyzs(:,1:2) nans]';
    ys = [xyzs(:,3:4) nans]';
    zs = [xyzs(:,5:6) nans]';
    plot3(xs(:),ys(:),zs(:),'color',ccol)
    hold on
  end%if
end%for c
hold off
clim(cLim) % set color range used after the plots
hand = gca().Children;

function xyz=csegment(i,j,k)
% Sub-function to compute contour line segment between grid
% links (i,j) and (i,k), where i,j,k are 1D indices of the
% 2D grid of surface points.  AJR, 4 Oct 2023
if isempty(i), xyz=[]; else
dj = (c-C(i))./(C(j)-C(i));
dk = (c-C(i))./(C(k)-C(i));
xyz = cat(2 ,X(i)+(X(j)-X(i)).*dj ,X(i)+(X(k)-X(i)).*dk ...
            ,Y(i)+(Y(j)-Y(i)).*dj ,Y(i)+(Y(k)-Y(i)).*dk ...
            ,Z(i)+(Z(j)-Z(i)).*dj ,Z(i)+(Z(k)-Z(i)).*dk );
end%if empty
end%function csegment

end%function contour4

