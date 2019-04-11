A = [1 -1; -2 1];
b = [1 2]';
c = [1 1];
ini = 0;
iter = 0;
BA = sym('w',[1 n]); % Variables Basis
NA = sym('x',[1 n]); % Variables Constrain
NB = [NA BA]; 
A=-A;

while max(c) > 0,%choose largest coefficient
[cj, col] = max(c);%choose largest coefficient
%||| [cj, col] =  find(c > 0,1,'first') |||
% If bland's pivoting rule
Acol = A(:,col);
[i, row] = max(-Acol./b); %select leaving variable
if i < 0;
    opt = -1; %unbounded
    'unbounded'
    break;
end
Arow = A(row,:);%A matrix
a = A(row,col);
A = A - Acol*Arow/a;
A(row,:) = -Arow/a;
A(:,col) = Acol/a;
A(row,col) = 1./a; %A after pivoting 
brow = b(row);% b matrix
b = b - Acol*(brow)./a;
b(row) = -brow./a; %b after pivoting
b(row)*c(col);
P = b(row)*c(col) + ini; % P constraint add constant value. 
ini = P;
ccol = c(col);%c matrix 
c = c - ccol*Arow./a;
c(col) = ccol./a;
B = b; %Create D (dictionary) matrix
B(n+1) = P;
D = vertcat(A,c);

D = horzcat(B,D) %Print D matrix
BA(row) = NA(col)  %Print Basis
N = setdiff(NB, BA)  %Print Constrain

iter = iter+1
end