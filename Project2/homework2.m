A = [1 3 0 1; 2 1 0 0; 0 1 4 1];
c = [1 4 1 3]';
b = [4 3 3]';
iter = 0;
[m,n] = size (A);
d = 0.1; %d=0.1 and d=0.5 
r = 0.99; %r=0.8 and r=0.99
tol = 10^-12;
M = 10^8;
kmax = 1000;
w = ones(m,1);
x = ones(n,1);
y = ones(m,1);
z = ones(n,1);
zeta = b - A*x - w; %primal
sigma = c- A'*y + z; %dual
gamma = (x'*z +y'*w)/(n+m); %change
normzeta = norm(zeta, inf); %primal norm
normsigma = norm(sigma, inf); %dual norm
normgamma = norm(gamma, inf); %change norm

for iter=1:5
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

G = [zeta; sigma;delta - Xz; delta - Yw]
X*Z;
delta*diag(4);
X*z;
x'*z + y'*w;

delta - X*z;
delta - Y*w;
dv = linsolve (F,G);

dx = dv(1:n)
dw = dv(n+1:m+n)
dy = dv(n+m+1:2*m+n)
dz = dv(2*m+n+1:2*m+2*n)

xi = max(max(-dx./x));
yi = max(max(-dy./y));
wi = max(max(-dw./w));
zi = max(max(-dz./z));
di = [xi yi wi zi];
di = max(di);
theta = min( 1,r*(1./di));

x = x + theta * dx;
w = w + theta * dw;
y = y + theta * dy;
z = z + theta * dz;

zeta = b -A*x - w;
sigma =c- A'*y+z;
gamma = (x'*z +y'*w)./(n+m);

convzeta = norm(zeta, inf);
convsigma = norm(sigma, inf);
convgamma = norm(gamma, inf);
%xs = [3.200000000000000e+00 2.000000000000000e-01]

%fnorm = norm(x - xs, inf)./norm(x,inf)
if(convzeta <= tol*normzeta && convsigma <= tol*normsigma && convgamma <= tol*normgamma)
break;
end
if(norm(x,inf) > M || norm(y,inf) > M)
    fprintf('!!UNBOUNDED!!\n')
    break;
end
iter = iter+1
end
