function [x,y,w,z] = pathfollow(A,b,c)
% syntax: [x,y,w,z] = pdip(A,b,c)
%
% path-following primal-dual interior-point method for problem 
%
% PRIMAL: max c'x s.t. Ax+w=b, x>=0, 
% DUAL:   min b'y  s.t. A'y-z=c,  z>=0.
%
% input: A is an m x n SPARSE constraint matrix.
%        b is an m x 1 right-hand side vector
%        c is an n x 1 cost vector.
%
% output: x is the  n x 1 solution of the primal problem
%         y is the m x 1 dual solution
%         w is the n x 1 vector of "dual slacks"
%         f is the optimal objective value
%
% internal parameters: 
%         kmax is the maximum number of iterations allowed
%         tol is the convergence tolerance
%         M is the factor used to define the starting point
%         maxDiag is the element for the X^{-1}S matrix
%         zetaMin is the minimum value of the steplength scale parameter eta

% check validity of input arguments

[m,n] = size(A);
if m <= 0 || n <= 0
  error('input matrix A must be nontrivial');
end
if n ~= length(c)
  error('size of vector p must match number of columns in A');
end
if m ~= length(b)
  error('size of vector b must match number of rows in A');
end

% set the internal parameters
kmax = 1000; 
tol = 10.e-12;
M = 10.e+10;
maxDiag = 5.e+15;
zetaMin = .8;

% start the clock
t0=cputime;

% set initial point, based on largest element in (A,b,p)
bigM = max(max(abs(A)));
bigM = max([norm(b,inf), norm(c,inf), bigM]);
x = M*bigM*ones(n,1); w = x; y = zeros(m,1);

bc = 1+max([norm(b), norm(c)]);

for iter=1:kmax
  
% compute residuals
  sigma =c- A'*y+z;
  zeta = b- A*x - w;
  Rc = x.*w;
  gamma = mean(Rc);

% check relative decrease in residual, for purposes of convergence test
  relResidual = norm([sigma;zeta;Rc])/bc;
  fprintf(1,'iter %2i: mu = %9.2e, resid = %9.2e\n', iter, full(gamma), ...
	  full(relResidual));

% test for convergence
  if(relResidual <= tol & gamma <= tol) break; end;

% make a heuristic choice of the centering parameter, and adjust the 
% right-hand side
  sigma = min(0.1,100*gamma);
  Rc = Rc - sigma*gamma;
  
  % set up the scaling matrix and form the coef matrix for normal equations
  c = min(maxDiag, x./w);
  B = A*sparse(1:n,1:n,c)*A';
  % use the form of the Cholesky routine "cholinc" that's best
  % suited to interior-point methods
  R = cholinc(B(ordering,ordering),'inf');
  
  % set up the right-hand side
  t1 = x.*sigma-Rc;
  t2 = -(zeta+A*(t1./w));
  
  % solve the normal equations system for dy and recover the other 
  %  step components dx and ds
  dy = zeros(m,1);
  dy(ordering) = R\(R'\t2(ordering));
  dx = (x.*(A'*dy)+t1)./w;
  ds = -(w.*dx+Rc)./x;
 
  % set the parameter eta defining fraction of max step to boundary
  eta = max(zetaMin,1-gamma);
  [alpha, alphax, alphas] = steplength(x, w, dx, ds, eta);

  % take the step
  x = x + alphax * dx;
  w = w + alphas * ds;
  y = y + alphas * dy;
end

% calculate final primal objective value
z = c'*x;

% convert x,y,s to full data structures
x=full(x); w=full(w); y=full(y);

fprintf('Done!\t[m n] = [%g %g]\tCPU = %g\n', m, n, cputime-t0);
return;  