% test chebyshev 'boundary layers' of full-domain separated
% by a central region of chebyshev distribution.  Arises
% because a chebyshev distribution concentrates points near
% the boundary and so otherwise the patches may well
% overlap.   AJR, 5 Jan 2023
nSubP=5
dx=0.04
for nP=2:9
    Xlim=cumsum(round(20*rand(1,2))/10)
    p=configPatches1(@sin,[0 1],'chebyshev' ...
        ,nP,0,dx,nSubP,'EdgyInt',false);
    xs=squeeze(p.x)
end
