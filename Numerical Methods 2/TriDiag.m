function x = TriDiag(A,d)

% Getting the size of all the matrices and initializing vars, matrices
[~,n]=size(A);
m=[];
c=[];
m(1)=A(1,1);
U=zeros(n);
L=zeros(n);

% Coefficients in the factorization
for j=1:n-1
    % initialize the b's and c's of the A matrix
    % for clarity's sake
    b(j)=A(j,j+1);
    c(j)=A(j+1,j);
    l(j)=c(j)/m(j);
    m(j+1)=A(j+1,j+1)-l(j)*b(j);
    
    % Define Upper and Lower Matrices
    U(j,j+1)=b(j);
    U(j,j)=m(j);
    L(j+1,j)=l(j);
    L(j,j)=1;
end
U(n,n)=m(n);
L(n,n)=1;

% Forward sub on Ly=d (from previous HW)
y=forwardsub(L,d);

% Backward sub on Ux=y (from previous HW)
x=backwardsub(U,y);


end

