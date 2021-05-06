function BestApprox = Romberg(tol,a,b,f)
h=(b-a);
R=[];
R(1,1)=0; % to initialize the loop
R(2,2)=(.5)*(b-a)*(f(a)+f(b));
l=1; % l is the tol variable
while abs(R(l+1,l+1)-R(l,l))>tol
    h=h/2;
    Summa=0;
    for j=1:2^(l-1)
        Summa = Summa + f(a+(2*j-1)*h);
    end
    
    R(l+2,2)=(1/2)*(R(l+1,2))+h*Summa;
    m=1;
    while m<=l
        R(l+2,m+2)=R(l+2,m+1)+(1/(4^(m+2)-1))*(R(l+2,m+1)-R(l+1,m+1));
        m=m+1;
    end
    BestApprox=R(l+2,l+2);
    l=l+1;
end



