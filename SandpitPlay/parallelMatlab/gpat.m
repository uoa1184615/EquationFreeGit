function y=gpat(x)
% also setfield | isfield | rmfield
global pat
y=isfield(pat,x);
if isfield(pat,x), y=getfield(pat,x); end
end
