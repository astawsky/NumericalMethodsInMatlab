function [U] = Thomas(A,F)

[n,~]=size(A);

for i=1:n-1
    
    % The new coeffiecients
    if i==1
        C(i)=A(i,i+1)/A(i,i);
        D(i)=F(i)/A(i,i);
    else
        C(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*C(i-1));
        D(i)=(F(i)-A(i,i-1)*D(i-1))/(A(i,i)-A(i,i-1)*C(i-1));
    end
    
end

D(n)=(F(n)-A(n,n-1)*D(n-1))/(A(n,n)-A(n,n-1)*C(n-1));

% Back sub to solve for U
U(n)=D(n);

for i=n-1:-1:1
    U(i)=D(i)-C(i)*U(i+1);
end

end

