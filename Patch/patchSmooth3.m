% legacy interface patchSmooth3() auto-invokes new patchSys3()
function dudt=patchSmooth3(t,u,patches)
global smOOthCount
if isempty(smOOthCount), smOOthCount=1; 
   else smOOthCount=smOOthCount+1; end
l2=log2(smOOthCount);
if abs(l2-round(l2))<1e-9
   warning('Use new patchSys3 instead of old patchSmooth3')
end
if nargin<3, global patches, end
dudt=patchSys3(t,u,patches);
