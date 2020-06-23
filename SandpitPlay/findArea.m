function a = findArea(width,varargin)
% findArea(width)
% findArea(width,height)
% findArea(... 'shape',shape)
p = inputParser; 
%addRequired(p,'width');
addRequired(p,'width',@isnumeric);
%defaultHeight = 1; 
addOptional(p,'height',1,@isnumeric);
%defaultShape = 'rectangle';
checkString = @(s) any(strcmp(s,{'square','rectangle'})); addParamValue(p,'shape','rectangle',checkString);
parse(p,width,varargin{:}); 
width = p.Results.width
height = p.Results.height
shape = p.Results.shape
