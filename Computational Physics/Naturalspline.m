function [s,A,B,C,D] = Naturalspline(X,Y)
% X is the set of x values and Y is the corresponding set of y values
%  from the collection of given points that need interpolating
%   s is the final cubic spline and ABCD are the arrays containing
%    the coefficients of the respective cubic splines

h = []; % defining all the distances between consecutive x values
for j=1:(length(X)-1)
    h = [h, X(j+1)-X(j)];
end

d = []; % defining the RHS of tridiagonal matrix
for j=1:length(X)-2
    d = [d, (-6./h(j)).*(Y(j+1)-Y(j))+(6./h(j+1)).*(Y(j+2)-Y(j+1))];
end

a = [];
for j=2:(length(X)-1)
    a = [a, 2.*(h(j-1)+h(j))];
end
b = [];
for j=2:(length(X)-1)
    b = [b, h(j)];
end
c = [];
for j=2:(length(X)-1)
    c = [c, h(j)];
end

l = [];
m = [a(1)];
for j=1:length(X)-3
    l = [l, c(j)./m(j)];
    m = [m, a(j+1) - l(j).*b(j)];
end

y=[d(1)];
for j=2:length(X)-2
    y = [y, d(j) - l(j-1).*y(j-1)];
end

z = zeros(1,length(X));
z(length(z)-1)=(y(length(X)-2)./m(length(X)-2));
z(1)=0;
z(length(X)) = 0;
for j=length(X)-3:-1:1
    z(j+1) = (y(j)-b(j).*z(j+2))./m(j);
end

A = [];
for j=1:length(X)-1
    A = [A, (z(j+1)-z(j))./(6.*h(j))];
end
B = [];
for j=1:length(X)-1
    B = [B, (.5).*z(j)];
end
C = [];
for j=1:length(X)-1
    C = [C, (Y(j+1)-Y(j))./h(j) - ((z(j+1)+2.*z(j)).*h(j))./(6)];
end
D = [];
for j=1:length(X)-1
    D = [D, Y(j)];
end

syms s(u)
s(u) = piecewise(X(1)<=u<X(2), A(1).*(u-X(1)).^3 + B(1).*(u-X(1)).^2 + C(1).*(u-X(1)) + D(1));
for j=2:length(X)-1
    if j == length(X)-1
        s(u) = piecewise(X(j)<=u<=X(j+1), A(j).*(u-X(j)).^3 + B(j).*(u-X(j)).^2 + C(j).*(u-X(j)) + D(j), s);
    end
    if j ~= length(X)-1
        s(u) = piecewise(X(j)<=u<X(j+1), A(j).*(u-X(j)).^3 + B(j).*(u-X(j)).^2 + C(j).*(u-X(j)) + D(j), s);
    end
end

end

