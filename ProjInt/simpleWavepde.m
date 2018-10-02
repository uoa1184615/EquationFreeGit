function Ut=simpleWavepde(t,U,x)
% global patches
dx=x(2)-x(1);
Ut=nan(size(U));
ht=Ut;
%{
\end{matlab}
Compute the PDE derivatives at points internal to the patches.
\begin{matlab}
%}
i=2:size(U,1)-1;
%{
\end{matlab}
Here `wastefully' compute time derivatives for both \pde{}s at all grid points---for `simplicity'---and then merges the staggered results.
Since \(\dot h_{ij} \approx-(u_{i+1,j}-u_{i-1,j})/(2\cdot dx) =-(U_{i+1,j}-U_{i-1,j})/(2\cdot dx)\) as adding\slash subtracting one from the index of a \(h\)-value is the location of the neighbouring \(u\)-value on the staggered micro-grid.
\begin{matlab}
%}
ht(i,:)=-(U(i+1,:)-U(i-1,:))/(2*dx);
%{
\end{matlab}
Since \(\dot u_{ij} \approx-(h_{i+1,j}-h_{i-1,j})/(2\cdot dx) =-(U_{i+1,j}-U_{i-1,j})/(2\cdot dx)\) as adding\slash subtracting one from the index of a \(u\)-value is the location of the neighbouring \(h\)-value on the staggered micro-grid.
\begin{matlab}
%}
Ut(i,:)=-(U(i+1,:)-U(i-1,:))/(2*dx);
%{
\end{matlab}
Then overwrite the unwanted~\(\dot u_{ij}\) with the corresponding wanted~\(\dot h_{ij}\).
\begin{matlab}
%}
% Ut(patches.hPts)=ht(patches.hPts);
Ut(1:2:end) = ht(1:2:end);
Ut([1 end]) = 0;
end