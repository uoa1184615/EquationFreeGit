%{
Attempt to use PIG in a telescoping fashion on problems with three or more
separated time scales.
Compare to the application of
PIG without telescoping.
%}
clear; close all

%Decide whether to compare to normal PIG
compareToNormalPIG = 1;


addpath('../ProjInt')

e1 = 1e-3; %first fast time scale
e2 = 1e-6; %second fast time scale

ny1 = 1; %number of fast variables on first scale
ny2 = 1; 

ee1 = e1 * ( 1+1e-1*randn(ny1,1) ); %spread out fast time scales slightly
ee2 = e2 * ( 1+1e-1*randn(ny2 ,1) );

dxdt = @(t,x) [-sum(x(2+ny1:end)) + x(1).*sum(x(2:1+ny1)); %slow variable
               (cos(x(1)) - x(2:1+ny1))./ee1; %first fast time scale
               (sin(x(1)) - x(2+ny1:end))./ee2; %second fast time scale
               ];

disp 'Output eigenvalues in the dynamical system:'
eigs = [sort(-1./ee2(:)') sort(-1./ee1(:)') 0] %roughly correct

           


%set up microsolver for slowest of the fast time scales
b1 = 2*max(ee1)*log(1/max(ee1));
del2 = min(ee2)/2; %time step
b12 = linspace(0, b1,1+ceil(b1/del2)); %bridge b1 with time steps <= del2
fb =  @(tb0, xb0) feval('rk2',dxdt,tb0+b12,xb0); %'fast burst'

tSpan = [0 5];
x0 = randn(1+ny1+ny2,1);

if compareToNormalPIG
    disp 'Simulate PIG normally (and time it):'
    tic
    [ts1,xs1,tms1,xms1] = PIG('ode45',fb,tSpan,x0); %boring ordinary PIG
    timeForNormalPIG = toc
end

disp 'Try a telescoping PIG (and time it):'

%set up microsolver for fastest time scales
b2 = 2*max(ee2)*log(1/max(ee2));
b22 = linspace(0, b2,1+ceil(b2/del2)); %bridge b2 with time steps <= del2
rfb = @(tb0, xb0) feval('rk2',dxdt,tb0+b22,xb0(:)); %'really fast burst'

del1 = min(ee1)/2;
b11 = linspace(0, b1,1+ceil(b1/del1)); %bridge b1 with time steps <= del1

%%%telescope with PIG
fPb =  @(tb0, xb0) feval('PIG','rk2',rfb,tb0 + b11,xb0,[],[],'cancelChristmas'); %'fast PI burst'

%%%telescope with PIRK2
% rfb = @(tb0, xb0,bObsolete) rfb(tb0,xb0);
% fPb =  @(tb0, xb0) feval('PIG',rfb,tb0 + b11,xb0,b2); %'fast PI burst'

tic
[ts2,xs2,tms2,xms2] = PIG('ode45',fPb,tSpan,x0);
timeForCoolPIG = round(toc,1)

ax1 = min(min(tms2))-0.1;
ax2 = tSpan(end)+0.1;


figure;
subplot(3,1,1) %slow variable
plot(ts2,xs2(:,1),'o')
xlim([ax1 ax2])
xlabel('Time $t$')
ylabel('Solution $x$')
title('Simulation by TelePIG')


subplot(3,1,2) %first fast variable
plot(tms2,xms2(:,2:1+ny1),'.')
xlim([ax1 ax2])
xlabel('Time $t$')
ylabel('Solution $y_1$')

subplot(3,1,3) %second fast variable
plot(tms2,xms2(:,2+ny1:end),'.')
xlim([ax1 ax2])
xlabel('Time $t$')
ylabel('Solution $y_2$')


if compareToNormalPIG
    figure;
    subplot(3,1,1) %slow variable
    plot(ts1,xs1(:,1),'o')
    xlim([ax1 ax2])
    xlabel('Time $t$')
    ylabel('Solution $x$')
    title('Simulation by PIG')
    
    subplot(3,1,2) %first fast variable
    plot(tms1,xms1(:,2:1+ny1),'.')
    xlim([ax1 ax2])
    xlabel('Time $t$')
    ylabel('Solution $y_1$')
    
    subplot(3,1,3) %second fast variable
    plot(tms1,xms1(:,2+ny1:end),'.')
    xlim([ax1 ax2])
    xlabel('Time $t$')
    ylabel('Solution $y_2$')
    
    figure; hold on
    title('Simulated macroscale solutions')
    oP = plot(ts1,xs1,'rx');
    tP = plot(ts2,xs2,'bo');
    legend([oP(1),tP(1)],'Ordinary PIG','TelePIG')
end





               
               