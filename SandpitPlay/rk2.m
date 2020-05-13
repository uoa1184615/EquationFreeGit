% Second order Runge-Kutta. No adaptive stepping, no error control. Use when
% testing other schemes.
% Inputs: 
% differential equation dxdt = dxdt(t,x)
% vector of times tSpan
% initial condition xIC
function [tOut,xOut] = rk2(dxdt,tSpan,xIC) 
deltas = diff(tSpan);
tOut = tSpan;
xOut = zeros(length(tSpan),length(xIC));
xOut(1,:) = xIC;
xIC = xIC(:);

for n = 1:length(deltas)
    del = deltas(n);
    t = tSpan(n);
    k1 = dxdt(t,xIC);
    k2 = dxdt(t+del,xIC+del*k1);
    xIC = xIC+del*(k1+k2)/2;
    xOut(n+1,:)=xIC;
end