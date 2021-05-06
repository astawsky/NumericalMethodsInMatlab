function x = backwardsub(U,b)
    [n,~]=size(U);
    x=zeros(n,1);
    for j=n:-1:1
        x(j)=(-U(j,:)*x+b(j))/U(j,j);
    end
end