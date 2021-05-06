function F = complexderivcoeffs(N)

% X is the array of the equidistant nodes
X=[];
for j=0:N-1
    X=[X,(j.*2.*pi)./N];
end

syms f(x) g(x)
f(x) = exp(sin(x));
g(x) = diff(exp(sin(x)));

% C contains the complex coefficients
C=[];
for j=1:N
    C = [C,fft(exp(sin(X(j))))];
end

% C1 contains the derivatives of the complex coefficients
C1 = [];

% split the sum into two sums because of indices

for k=(-N/2):0
    C1 = [C1, (1i*k)*C(k+N)];
end
for k=1:N/2
    C1 = [C1, (1i*k)*C(k)];
end

% F contains the approximated value points of the derivative of f(x)
F = ifft(C1);

end

