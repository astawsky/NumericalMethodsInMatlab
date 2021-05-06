function x = NewtonsMethod1(x,g,dg,tol)

error=100;
while error>tol
    previous=x;
    x=previous+g(previous)/dg(previous);
    error=abs(x-previous);
end

end
