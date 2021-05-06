function [L,P] = LowerUpper(A)

if det(A)==0
    disp('A is non-singular!');
else
    P=A; % P is a dummy A
n=size(A);
L=eye(n); % initialize the LTM
for i=1:n-1 % number of reductions, L_i matrices, columns
    K=eye(n); % K is re-initialized
    for j=i+1:n % number of m_{i,j}s, rows
        K(j,i)=(-P(j,i)/P(i,i)); % K is the L_i LTM
        L(j,i)=(P(j,i)/P(i,i)); % Adding the same elements to overall L
    end
    P=K*P; % Reducing K by that LTM
end
end

end