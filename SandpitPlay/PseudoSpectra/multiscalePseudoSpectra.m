% multiscalePseudoSpectra()  hierarchially contours the
% 2-norm pseudospectra of a given matrix A.  For matrices A
% obtained from multiscale systems there are often big
% boring gaps in the spectrum which it is useless to bother
% analysing.  So this hierarchial version heuristically
% concentrates on 'interesting' regions.  It is based upon
% psa.m which is typically about N/4 times faster than the
% obvious SVD method.  Comes with no guarantees! L.N.
% Trefethen, March 1999.  Ditto, AJR, 18 May 2023

function multiscalePseudoSpectra(A,sc,maxLevel,nInt) 

% Inputs:
% A - the matrix to contour the pseudospectra
if nargin<1, A=randn(20); disp("Using matrix randn(20)"), end
%
% sc - (optional, default Inf) parameter(s) of nonlinear
% scaling in the plot.  Let s:=|eigenvalue|, then for s<sc
% the plot is approximately linear, and for s>sc the plot is
% approximately logarithmic.
if nargin<2, sc=Inf; end
%
% maxLevel - (optional, default 2) depth of recursive refinement
if nargin<3, maxLevel=2; end
%
% nInt - (optional, default 22) the grid resolution of each level
% maxLevel=3; nInt=10; % similar resolution and operation count
if nargin<4, nInt=22; end

% Output: in the current figure, plots spectrum and plots
% the pseudospectra contours at levels 10^n for n in log10Vals

% Parameters:  
% fracZoom is the fraction of mesh points to refine about.
% With fracZoom=1/nInt, operation count is roughly propto
% nInt^(maxLevel+1)
fracZoom = 1/nInt; 
log10Vals = -8:2; 


% establish independent scaling in each direction
if numel(sc)==1, sc=[sc sc]; end
if sc(1)==Inf, xSpace = @(a,b,n) linspace(a,b,n+1);
else xSpace = @(a,b,n) sc(1)*sinh( ...
              linspace(asinh(a/sc(1)),asinh(b/sc(1)),n+1) );
end
if sc(2)==Inf, ySpace = @(a,b,n) linspace(a,b,n+1);
else ySpace = @(a,b,n) sc(2)*sinh( ...
              linspace(asinh(a/sc(2)),asinh(b/sc(2)),n+1) );
end

% Compute Schur form, from psa.m
  [U,T] = schur(A);
  if isreal(A), [U,T] = rsf2csf(U,T); end
  T = triu(T); eigA = diag(T);

% Initial plot of eigenvalues as dots
  s = 0.08*norm(A,1);
  xmin = min(real(eigA))-s;  ymin = min(imag(eigA))-s;
  xmax = max(real(eigA))+s;  ymax = max(imag(eigA))+s; 
  clf, plot(real(eigA),imag(eigA),'.'), hold on
  axis([xmin xmax ymin ymax]), grid on, drawnow
  title('spectrum (dots) and pseudo-spectra contours, via psa()') % AJR

% From psa.m, reorder Schur decomposition and compress to
% interesting subspace:
%  select = find(real(eigA)>-250) % <- ALTER SUBSPACE SELECTION
  select = 1:length(eigA); % Hmmm, lets do all subspaces??
  n = length(select);
  for ii = 1:n
    for kk = select(ii)-1:-1:ii
      G([2 1],[2 1]) = planerot([T(kk,kk+1) T(kk,kk)-T(kk+1,kk+1)]')';
      J = kk:kk+1; T(:,J) = T(:,J)*G; T(J,:) = G'*T(J,:);
    end
  end
  T = triu(T(1:n,1:n)); I = eye(n);
 
  % Set up mesh grid M{1} for coarse contour plot:
  M{1}.x = xSpace(xmin,xmax,nInt);
  M{1}.y = ySpace(ymin,ymax,nInt)';
  
  drawLevel(1)
  hold off
  xlabel('$\Re\lambda$'), ylabel('$\Im\lambda$')
  colormap(0.8*parula(numel(log10Vals)))
  % designed to use quasiLogAxes(), but need not 
  if exist('quasiLogAxes')==2, quasiLogAxes(sc(1),sc(2)); end
  return% end of this function execution



% resursive function to compute and draw the contours
function drawLevel(l)

  % iteration code from psa.m  Note i is y-thing, j is x-thing
  Z = M{l}.x+1i*M{l}.y;
  sigmin = nan(size(Z));
  for i = 1:length(M{l}.y)
    if isreal(A) & (abs(sum(M{l}.y([1 end])))<1e-9) & (i>(length(M{l}.y)+1)/2)
      sigmin(i,:) = sigmin(length(M{l}.y)+1-i,:);
    else
      for j = 1:length(M{l}.x)
        z = Z(i,j); T1 = z*I-T; T2 = T1';
        if real(z)<1e9 % was 100      % <- ALTER GRID POINT SELECTION
          sigold = 0; qold = zeros(n,1); beta = 0; H = [];
          q = randn(n,1) + sqrt(-1)*randn(n,1); q = q/norm(q);
          for k = 1:99
            v = T1\(T2\q) - beta*qold;
            alpha = real(q'*v); v = v - alpha*q;
            beta = norm(v); qold = q; q = v/beta;
            H(k+1,k) = beta; H(k,k+1) = beta; H(k,k) = alpha;
            if (alpha>1e100), sig = alpha; 
            else sig = max(eig(H(1:k,1:k))); 
            end%if
            if (abs(sigold/sig-1)<.001) | (sig<3 & k>2), break, end
            sigold = sig;
          end%for k
          sigmin(i,j) = 1/sqrt(sig);
        end%if real(z)
      end
    end%if
  end%for i
  sig0 = 10^(min(log10Vals)-2); % avoid log(zero)
  M{l}.h = log10(sigmin+sig0);
  
  % conditionally recursively refine
  if l<maxLevel
    % proceed to heuristically find grid indices i,j to refine
    i=2:nInt;
    curv=M{l}.h(i+1,i)+M{l}.h(i-1,i)+M{l}.h(i,i+1)+M{l}.h(i,i-1)-4*M{l}.h(i,i);
    [curvSort,ij] = sort(curv(:),'Descend');
    K=ceil(fracZoom*numel(curv));
    [i,j]=ind2sub([1 1]*(nInt-1),sort(ij(1:K)));
    i=i+1; j=j+1;
    % for each of these points 
    for k=1:K
      % mark point as nan at this level so no contours close
      M{l}.h(i(k),j(k))=nan;
      % Set the refined grid to use at next level.  When
      % neighbours have been refined, then halve the domain
      % of this refinement. Due to ordering of sorted ij,
      % just test N and W neighbours 
      if M{l}.h(i(k)-1,j(k))==nan
           yMin=M{l}.y(i(k)); ny=nInt/2;
      else yMin=M{l}.y(i(k)-1); ny=nInt;
      end
      if M{l}.h(i(k),j(k)-1)==nan
           xMin=M{l}.x(j(k)); nx=nInt/2;
      else xMin=M{l}.x(j(k)-1); nx=nInt;
      end
      M{l+1}.y = ySpace(yMin,M{l}.y(i(k)+1),ny)';
      M{l+1}.x = xSpace(xMin,M{l}.x(j(k)+1),nx);
      % recurse to next level
      drawLevel(l+1)   
    end% for k
  end% if l
  
  % now plot the contours for this level (so coarsest last)
  contour(M{l}.x,M{l}.y,M{l}.h,log10Vals,'ShowText','on');
end% function drawLevel
  
end% function multiscalePseudoSpectra
