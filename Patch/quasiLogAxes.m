% quasiLogAxes() transforms selected axes of the given plot
% to a quasi-log axes (via asinh), including possibly
% transforming color axis.  AJR, 25 Sep 2021 -- 18 Apr 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{quasiLogAxes()}: transforms some axes of a
plot to quasi-log}
\label{sec:quasiLogAxes}

This function rescales some coordinates and labels the axes
of the given 2D or 3D~plot.  The original aim was to
effectively show the complex spectrum of multiscale systems
such as the patch scheme. The eigenvalues are over a wide
range of magnitudes, but are signed.  So we use a nonlinear
asinh transformation of the axes, and then label the axes
with reasonable ticks. The nonlinear rescaling is useful in
other scenarios also.

\begin{matlab}
%}
function quasiLogAxes(handle,xScale,yScale,zScale,cScale)
%{
\end{matlab}
\paragraph{Input} 
\begin{itemize}
\item \verb|handle|: handle to your plot to transform, for
example, obtained by \verb|handle=plot(...)|

\item \verb|xScale| (optional, default~inf): if inf, then no
transformation is done in the `x'-coordinate.  Otherwise, when
\verb|xScale| is not inf, transforms the plot \(x\)-coordinates 
with the \(\text{asinh}()\) function so that
\begin{itemize}
\item for \(|x|\lesssim x_{\text{scale}}\) the x-axis scaling 
is approximately linear, whereas
\item for \(|x|\gtrsim x_{\text{scale}}\) the x-axis scaling 
is approximately signed-logarithmic.
\end{itemize}

\item \verb|yScale| (optional, default~inf): corresponds to
\verb|xScale| for the second axis scaling.

\item \verb|zScale| (optional, default~inf): corresponds to
\verb|xScale| for a third axis scaling if it exists.

\item \verb|cScale| (optional, default~inf): corresponds to
\verb|xScale| but for a colormap, and colorbar scaling if 
one exists.
\end{itemize}

\paragraph{Output} None, just the transformed plot.


\paragraph{Example}
If invoked with no arguments, then execute an example.
\begin{matlab}
%}
if nargin==0  
   % generate some data
   n=99;  fast=(rand(n,1)<0.8);
   z = -rand(n,1).*(1+1e3*fast)+1i*randn(n,1).*(5+1e2*fast);
   % plot data and transform axes
   handle = plot(real(z),imag(z),'.');
   xlabel('real-part'), ylabel('imag-part')
   quasiLogAxes(handle,1,10);
   return
end% example
%{
\end{matlab}


Default values for scaling, \verb|inf| denotes no
transformation of that axis.
\begin{matlab}
%}
if nargin<5, cScale=inf; end
if nargin<4, zScale=inf; end
if nargin<3, yScale=inf; end
if nargin<2, xScale=inf; end
%{
\end{matlab}


\begin{devMan}
Get current limits of the plot to use if the user has set
them already.  And also get the pointer to the axes and to
the figure of the plot. 
\begin{matlab}
%}
xlim0=xlim; ylim0=ylim; zlim0=zlim; clim0=caxis;
theAxes = get(handle(1),'parent');
theFig = get(theAxes,'parent');
%{
\end{matlab}

Find overall factors so the data is nonlinearly mapped to
order oneish---so that then pgfplots et al.\ do not think
there is an overall scaling factor on the axes.
\begin{matlab}
%}
xFac=1e-99; yFac=xFac; zFac=xFac; cFac=xFac;
for k=1:length(handle)
    if ~isinf(xScale)
    temp = asinh(handle(k).XData/xScale);
    xFac = max(xFac, max(abs(temp(:)),[],'omitnan') );
    end
    if ~isinf(yScale)
    temp = asinh(handle(k).YData/yScale);
    yFac = max(yFac, max(abs(temp(:)),[],'omitnan') );
    end
    if ~isinf(zScale)
    temp = asinh(handle(k).ZData/zScale);
    zFac = max(zFac, max(abs(temp(:)),[],'omitnan') );
    end
    if ~isinf(cScale)
    temp = asinh(handle(k).CData/cScale);
    cFac = max(cFac, max(abs(temp(:)),[],'omitnan') );
    end
end%for
xFac=9/xFac; yFac=9/yFac; zFac=9/zFac; cFac=9/cFac;
%{
\end{matlab}

Scale the plot data in the plot \verb|handle|. Give an
error if it appears that the plot-data has already been
transformed.   Color data has to be transformed first
because usually there is automatic flow from z-data to c-data.
\begin{matlab}
%}
for k=1:length(handle)
    assert(~strcmp(handle(k).UserData,'quasiLogAxes'), ...
       'Replot graph---it appears plot data is already transformed')
    if ~isinf(cScale)
    handle(k).CData = cFac*asinh(handle(k).CData/cScale); 
    clim1=[min(handle(k).CData(:)) max(handle(k).CData(:))];
    end
    if ~isinf(xScale)
    handle(k).XData = xFac*asinh(handle(k).XData/xScale); 
    xlim1=[min(handle(k).XData(:)) max(handle(k).XData(:))];
    end
    if ~isinf(yScale)
    handle(k).YData = yFac*asinh(handle(k).YData/yScale); 
    ylim1=[min(handle(k).YData(:)) max(handle(k).YData(:))];
    end
    if ~isinf(zScale)
    handle(k).ZData = zFac*asinh(handle(k).ZData/zScale); 
    zlim1=[min(handle(k).ZData(:)) max(handle(k).ZData(:))];
    end
    handle(k).UserData = 'quasiLogAxes';
end%for
%{
\end{matlab}
Set 4\%~padding around all margins of transformed
data---crude but serviceable.  Unless the axis had already
been manually set, in which case use the transformed set
limits.
\begin{matlab}
%}
if ~isinf(xScale), 
    if xlim('mode')=="manual"
         xlim1=xFac*asinh(xlim0/xScale);
    else xlim1=xlim1+0.04*diff(xlim1)*[-1 1]; 
    end, end
if ~isinf(yScale), 
    if ylim('mode')=="manual"
         ylim1=yFac*asinh(ylim0/yScale);
    else ylim1=ylim1+0.04*diff(ylim1)*[-1 1];  
    end, end
if ~isinf(zScale), 
    if zlim('mode')=="manual"
         zlim1=zFac*asinh(zlim0/zScale);
    else zlim1=zlim1+0.04*diff(zlim1)*[-1 1];  
    end, end
if ~isinf(cScale), 
    if theAxes.CLimMode=="manual"
         clim1=cFac*asinh(clim0/cScale);
    else clim1=clim1+   0*diff(clim1)*[-1 1];  
    end, end
%{
\end{matlab}


\paragraph{Scale axes, and tick marks on axes}
\begin{matlab}
%}
if ~isinf(xScale)
    xlim(xlim1);
    tickingQuasiLogAxes(theAxes,'X',xlim1,xScale,xFac)
end%if
if ~isinf(yScale)
    ylim(ylim1);
    tickingQuasiLogAxes(theAxes,'Y',ylim1,yScale,yFac)
end%if
if ~isinf(zScale)
    zlim(zlim1);
    tickingQuasiLogAxes(theAxes,'Z',zlim1,zScale,zFac)
end%if
%{
\end{matlab}
But for color, only tick when we find a colorbar. 
\begin{matlab}
%}
if ~isinf(cScale)
  caxis(clim1);
  for p=1:numel(theFig.Children)
    ca = theFig.Children(p);
    if class(ca) == "matlab.graphics.illustration.ColorBar"
      tickingQuasiLogAxes(ca,'C',clim1,cScale,cFac)
      break  
    end
  end
end%if
%{
\end{matlab}
Turn the grid on by default.
\begin{matlab}
%}
grid on
end%function
%{
\end{matlab}


\subsection{\texttt{tickingQuasiLogAxes()}: typeset ticks
and labels on an axis}
\begin{matlab}
%}
function tickingQuasiLogAxes(ca,Q,qlim1,qScale,qFac)
%{
\end{matlab}
\paragraph{Input} 
\begin{itemize}
\item \verb|ca|: pointer to axes/colorbar dataset.
\item \verb|Q|: character, either \texttt{X,Y,Z,C}.
\item \verb|qlim1|: the scaled limits of the axis.
\item \verb|qScale|: the scaling parameter for the axis.
\item \verb|qFac|: the scaling factor for the axis.
\end{itemize}

\paragraph{Output} None, just the ticked and labelled axes.

Get the order of magnitude of the horizontal data.
\begin{matlab}
%}
    qmax=max(abs(qlim1));
    qmag=floor(log10(qScale*sinh(qmax/qFac)));
%{
\end{matlab}
Form a range of ticks, geometrically spaced, trim off the
small values that would be too dense near zero (omit those
within 6\% of \verb|qmax|).
\begin{matlab}
%}
    ticks=10.^(qmag+(-7:0));
    j=find(ticks>qScale*sinh(0.06*qmax/qFac));
    nj=length(j);
    if nj<3,     ticks=[1;2;5]*ticks(j);
    elseif nj<5, ticks=[1;3]*ticks(j);
    else         ticks=ticks(j);
    end
    ticks=sort([0;ticks(:);-ticks(:)]);
%{
\end{matlab}
Set the ticks in place according to the transformation.
\begin{matlab}
%}
    if Q=='C', p='s'; Q=''; else p=''; end
    set(ca,[Q 'Tick' p],qFac*asinh(ticks/qScale) ...
          ,[Q 'TickLabel' p],cellstr(num2str(ticks,4)))
    if Q=='X', set(ca,[Q 'TickLabelRotation'],40), end
end%function qScaling
%{
\end{matlab}
\end{devMan}
%}
