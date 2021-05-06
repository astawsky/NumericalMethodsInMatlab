function total = Comp_trap_rule(a,b,N,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h = (b-a)/N;

innersum = 0;
for j=1:N-1
    innersum = innersum + f(a+j*h);
end

total = h*((.5)*(f(a)+f(b))+innersum);

end

