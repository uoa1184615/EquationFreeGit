% psa.m - Simple code for 2-norm pseudospectra of given matrix A.
%         Typically about N/4 times faster than the obvious SVD method.
%         Comes with no guarantees!   - L. N. Trefethen, March 1999.
function psa(A) % made into function by AJR, Oct 2020

% Compute Schur form and plot eigenvalues:
  [U,T] = schur(A);
  if isreal(A), [U,T] = rsf2csf(U,T); end, T = triu(T); eigA = diag(T);

% Set up grid for contour plot:
  npts = 80;                                 % <- ALTER GRID RESOLUTION
  s = .8*norm(A,1);
%  s = 0.1*norm(A,1)+max(abs( [real(eigA);imag(eigA)] )); % AJR alteration
%  s=5 %%%%%%%%???? restricts window of analysis
  sc=1e-1; %%%%%%%????
  xmin = -s; xmax = s; ymin = -s; ymax = s;  % <- ALTER AXES
  xmin=-s/10+min(real(eigA)); xmax=+s/10+max(real(eigA)); % AJR
  ymin=-s/10+min(imag(eigA)); ymax=+s/10+max(imag(eigA)); % AJR
  x = linspace(xmin,xmax,npts);
  y = linspace(ymin,ymax,npts);
  x=sinh(linspace(asinh(xmin/sc),asinh(xmax/sc),npts))*sc; % AJR focus
  y=sinh(linspace(asinh(ymin/sc),asinh(ymax/sc),npts))*sc; % AJR focus
  [xx,yy] = meshgrid(x,y); zz = xx + sqrt(-1)*yy;
 
  hold off, plot(real(eigA),imag(eigA),'.'), hold on % AJR omit ,'markersize',15
  axis([xmin xmax ymin ymax]), axis square, grid on, drawnow
  title('spectrum (.) and pseudo-spectrum contours, via psa()') % AJR
 
% Reorder Schur decomposition and compress to interesting subspace:
  select = find(real(eigA)>-1e3);%-250            % <- ALTER SUBSPACE SELECTION
  n = length(select);
  for i = 1:n
    for k = select(i)-1:-1:i
      G([2 1],[2 1]) = planerot([T(k,k+1) T(k,k)-T(k+1,k+1)]')';
      J = k:k+1; T(:,J) = T(:,J)*G; T(J,:) = G'*T(J,:);
    end
  end
  T = triu(T(1:n,1:n)); I = eye(n);
 
% Compute resolvent norms by inverse Lanczos iteration and plot contours:
  sigmin = Inf*ones(length(y),length(x));
  for i = 1:length(y)
    if isreal(A) & (ymax==-ymin) & (i>length(y)/2)
      sigmin(i,:) = sigmin(length(y)+1-i,:);
    else
      for j = 1:length(x)
        z = zz(i,j); T1 = z*I-T; T2 = T1';
        if real(z)<100                       % <- ALTER GRID POINT SELECTION
          sigold = 0; qold = zeros(n,1); beta = 0; H = [];
          q = randn(n,1) + sqrt(-1)*randn(n,1); q = q/norm(q);
          for k = 1:99
            v = T1\(T2\q) - beta*qold;
            alpha = real(q'*v); v = v - alpha*q;
            beta = norm(v); qold = q; q = v/beta;
            H(k+1,k) = beta; H(k,k+1) = beta; H(k,k) = alpha;
            if (alpha>1e100), sig = alpha; else sig = max(eig(H(1:k,1:k))); end
            if (abs(sigold/sig-1)<.001) | (sig<3 & k>2) break, end
            sigold = sig;
            end
         %text(x(j),y(i),num2str(k))         % <- SHOW ITERATION COUNTS
          sigmin(i,j) = 1/sqrt(sig);
        end
      end
    end%AJR , disp(['finished line ' int2str(i) ' out of ' int2str(length(y))])
  end
  contour(x,y,log10(sigmin+1e-20),-8:1,'ShowText','on');    % <- ALTER LEVEL LINES
  colormap(0.8*parula)
  hold off;
  
end% function