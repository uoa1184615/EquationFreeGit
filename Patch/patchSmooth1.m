% legacy interface patchSmooth1() auto-invokes new patchSys1()
function dudt=patchSmooth1(t,u,patches)
global smOOthCount
if isempty(smOOthCount), smOOthCount=1; 
   else smOOthCount=smOOthCount+1; end
l2=log2(smOOthCount);
if abs(l2-round(l2))<1e-9
   warning('Use new patchSys1 instead of old patchSmooth1')
end
if nargin<3, global patches, end
dudt=patchSys1(t,u,patches);
