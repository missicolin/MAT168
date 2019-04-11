A = [1 3 0 1; 2 1 0 0; 0 1 4 1];
iter = 0;
[m,n] = size (A);
d = 0.1;
r = 0.99;
tol = 10^-12;
kmax = 1000;
c = [1 4 1 3]';
b = [4 3 3]';
w = ones(m,1);
x = ones(n,1);
y = ones(m,1);
z = ones(n,1);

zeta = b -A*x - w;
sigma = c- A'*y+z;
gamma = (x'*z +y'*w)/(n+m);

normzeta = norm(zeta, inf);
normsigma = norm(sigma, inf);
normgamma = norm(gamma, inf);

for iter=1:kmax
Z = diag(z);
X = diag(x);
Y = diag(y);
W = diag(w);  
delta = d*gamma;
Xz = X*z;
Yw = Y*w;

D=[A eye(m); zeros(n,n+m); diag(z) zeros(n, m);zeros(m,n) Y];

E=[zeros(m,m+n); A', -eye(n); zeros(n,m) X; W, zeros(m,n)];

F = [D E];

G = [zeta; sigma;delta - Xz; delta - Yw];

dv = linsolve (F,G);
dx = dv(1:n);
dw = dv(n+1:m+n);
dy = dv(n+m+1:2*m+n);
dz = dv(2*m+n+1:2*m+2*n);

theta = min( 1, r*(max(dv./G)));

x = x + theta * dx
w = w + theta * dw;
y = y + theta * dy;
z = z + theta * dz;

zeta = b -A*x - w;
sigma =c- A'*y+z;
gamma = (x'*z +y'*w)/(n+m);

convzeta = norm(zeta, inf);
convsigma = norm(sigma, inf);
convgamma = norm(gamma, inf);

if(convzeta <= tol*normzeta && convsigma <= tol*normsigma && convgamma <= tol*normgamma)
break;
end
iter = iter+1
end
