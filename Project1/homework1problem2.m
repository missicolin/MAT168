n=3;
A=zeros(n);
for x = 1:n
    for y = 1:n
        if x-y > 0;
            A(x,y) = 2*(10.^ [x-1]);
        else if x == y;
            A(x,y) =   1;  
        end
        end
    end
end
A
b= 100.^[0:n-1]';
c= 10.^[0:n-1];
c=fliplr(c);
iter = 0;
A=-A;
while max(c) > 0,%choose largest coefficient
[cj, col] = max(c);
Acol = A(:,col);
[i, row] = max(-Acol./b) %select leaving variable
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
b(row) = -brow./a

ccol = c(col);
c = c - ccol*Arow./a;
c(col) = ccol./a

iter = iter+1

end

