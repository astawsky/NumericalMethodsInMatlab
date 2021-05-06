function x = forwardsub(A,b)
    [n,~]=size(A);
    x=zeros(n,1);
    for j=1:n
        x(j)=(b(j)-A(j,:)*x)/A(j,j);
    end
end