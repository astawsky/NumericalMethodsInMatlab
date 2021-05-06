function P = BaryFormula(A,B,C,x)
%UNTITLED7 Summary of this function goes here
%   A is the set of nodes, B is the baryweights,
%   C is the set of function values at the nodes,
%   x is the evaluated point
P=0;
T=0;
Q=0;
for j=1:length(B)
    T = T + C(j).*(B(j)./(x - A(j)));
end

for j=1:length(B)
    Q = Q + (B(j)./(x - A(j)));
end

P = T/Q;



end

