% legacy interface patchSmooth2() auto-invokes new patchSys2()
function dudt=patchSmooth2(t,u,patches)
global smOOthCount
if isempty(smOOthCount), smOOthCount=1; 
   else smOOthCount=smOOthCount+1; end
l2=log2(smOOthCount);
if abs(l2-round(l2))<1e-9
   warning('Use new patchSys2 instead of old patchSmooth2')
end
if nargin<3, global patches, end
dudt=patchSys2(t,u,patches);
