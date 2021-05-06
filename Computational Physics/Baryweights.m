function B = Baryweights(A,e,c)
%UNTITLED5 Summary of this function goes here
%   A is the set of nodes
%   e=1 means the nodes are equidistributed
%   c=1 means the nodes are Chebyshev

B=[]; % The set of Bary weights

if e==1 % If the nodes are equispaced
    for j=1:length(A)
        Lambda = (1)^(j)*nchoosek(length(A)-1,j-1);
        B = [B, Lambda];
    end
end

if c==1 % If the nodes are Chebyshev
    for j=1:length(A)
        if j== 1 || length(A)
            Lambda = .5*(-1)^(j);
        end
        if j~=1 && j~=length(A)
            Lambda = (-1)^(j);
        end
        B = [B, Lambda];
    end   
end

if e==0 && c==0 % If the nodes are arbitrary
    for j=1:1:length(A)
    P=1;
    for k=1:length(A)
        if k ~= j
            P = P*(A(j) - A(k));
        end
    end
    Lambda = (1/P);
    B = [B, Lambda];
    end
end
end

