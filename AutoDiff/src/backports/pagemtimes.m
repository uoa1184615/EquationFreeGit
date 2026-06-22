function [result] = pagemtimes(x,y)
warning('using AutoDiff/src/backports/pagemtimes: are you sure?')%AJR
    result = multiprod(x,y);
end

