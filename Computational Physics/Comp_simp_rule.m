function total = Comp_simp_rule(a,b,N,f)

% N must be even!

h = (b-a)/N;

total = 0;
for j=0:2:(N-2)
    total = total + (h/3)*(f(a+j*h)+4*f(a+(j+1)*h)+f(a+(j+2)*h));
end

end

