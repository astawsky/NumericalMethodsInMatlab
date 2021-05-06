function [C,p] = NewtonDivDifCoeff(Nodes,f,x)
% Nodes is the set of nodes, f is the function we are
  % approximating and x is the value we are approximating at
    % C is the array of coefficients and p is the polynomial

C=[];
for i=1:length(Nodes) % all the zero divided differences
    C=[C, f(Nodes(i))];
end
for k=2:length(Nodes) % Computing the set of coefficients of the int. poly.
    for j=length(Nodes):-1:k
        C(j)=(C(j)-C(j-1))/(Nodes(j)-Nodes(j-k+1));
    end
end
p = C(length(Nodes)); % computing the integration polynomial p
for j=(length(Nodes)-1):-1:1
    p = C(j) + (x - Nodes(j))*p;
end

end

