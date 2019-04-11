n=4;
m=3;
ini=0;
c = [1,4,1,3];
A = [1 3 0 1; 2 1 0 0; 0 1 4 1];
b = [4,3,3]';
iter = 0;
BA = sym('w',[1 3]); % Variables Basis
NA = sym('x',[1 4]); % Variables Constrain
NB = [NA BA]; 
A=-A;
while max(c) > 0,
[cj, col] =  find(c > 0,1,'first')%Bland's pivoting rule
Acol = A(:,col);
[i, row] = max(-Acol./b); %select leaving variable
if i < 0;
    opt = -1; %unbounded
    'unbounded'
    break;
end
Arow = A(row,:);
a = A(row,col);
A = A - Acol*Arow/a;
A(row,:) = -Arow/a;
A(:,col) = Acol/a;
A(row,col) = 1./a;

brow = b(row);
b = b - Acol*(brow)./a;
b(row) = -brow./a;

ccol = c(col);
c = c - ccol*Arow./a
c(col) = ccol./a;

iter = iter+1

end




