function [x,iteration] = Gauss_Seidel(A,b)

[n,~]=size(A);

% defining W
W=tril(A);
for j=1:n
    W(j,j)=A(j,j);
end

N=n+1;
h=1/N;
x{1}=zeros(n,1);
j=1;
iteration=0;

while ~(max(b-A*x{j})<(.1*h))
    x{j+1}=(eye(n)-(eye(n)/W)*A)*x{j}+(eye(n)/W)*b';
    j=j+1;
    iteration=iteration+1;
end

x=x{j};


end

