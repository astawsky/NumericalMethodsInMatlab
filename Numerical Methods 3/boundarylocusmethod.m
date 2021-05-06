function boundarylocusmethod(alpha,beta)

% Boundary locus method
    m = length(alpha);
    theta = linspace(0, 2*pi, 200);
    n = length(theta);
    eVal = exp(1i*theta);

    num = zeros(1,n);
    den = zeros(1,n);

    % main for loop
    for j=1:m
        num = num + alpha(j)*eVal.^j;
        den = den + beta(j)*eVal.^j;
    end

    rho = num./den;
    x = real(rho);
    y = imag(rho);
    figure();
    plot(x,y);
    

end

