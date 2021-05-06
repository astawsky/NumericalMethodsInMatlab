function [A,F,TruVals,X,h] = CentralDiff(n)

% setting the second derivative and the exact solution
f=@(x) 9*sin(3*x)-6*x;
u=@(x) sin(3*x)+x^3;

% n is the amount of grid points
h=pi/(n-1);

for i=1:n-1
    
    % making the tridiagonal matrix A, (n-1)x(n-1) matrix
    A(i,i)=2*(1/h^2);
    if i~=1
        A(i,i-1)=-1*(1/h^2);
    end
    if i~=n-1
        A(i,i+1)=-1*(1/h^2);
    end
    
    % setting the grid points
    X(i)=i*h;
    
    % making the F vector
    if i==n-1
        F(i)=f(X(i))+pi^3/h^2;
    else
        F(i)=f(X(i));
    end
    
    % making the vector of true values
    TruVals(i)=u(X(i));
    
end

end

