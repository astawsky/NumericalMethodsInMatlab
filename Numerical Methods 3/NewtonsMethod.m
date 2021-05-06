function x = NewtonsMethod(x,t,g,dg,tol)

error=100;
while error>tol
    previous=x;
    x=previous+g(previous,t)/dg(previous,t);
    error=abs(x-previous);
end

end

