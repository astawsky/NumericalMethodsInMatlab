% input is A and we want to switch the ith and jth row
function P = PermMat(A,i,j)
    I=eye(size(A));
    P(:,j)=I(:,i);
    P(:,i)=I(:,j);
end