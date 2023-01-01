
w=linspace(0,5);
c=(730 + 110*w.^2 - 20*w.^4)./( - 100 - 51*w.^2 + 39*w.^4 - 9*w.^6 + w.^8);
c(abs(c)>100)=nan;
figure(5)
quasiLogAxes( plot(w,c) ,100,0.1)
