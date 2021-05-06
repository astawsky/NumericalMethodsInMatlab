% Alejandro Stawsky, Math 104C HW#2, April 28 2018
format long

% 1. k=.1,.05,.025, f(u,t)=-3tu-u^2, u(1)=0.5, u(2)=?

% % a. Forward Euler Method

N=[10,20,40];%partitions
K=[.1,.05,.025];%steps
U=[];%solutions
u=[];%collection of all solutions

for i=1:3%partition index
    t=1;%starting time
    for j=1:N(i)%step index
        t=t+K(i);%t_n
        U(1)=.5;%IC
        U(j+1)=U(j)-K(i)*(3*U(j)*t+(U(j))^2);%rest of solutions
        if j==N(i)
            u(i)=[U(j+1)];%collection of last solutions, \approx u(2)
        end
    end%step
end%partition

u%allsolutions

r1=(u(1)-u(2))/(u(2)-u(3))%ratio

% % b. Backward Euler Method
% % (We use the data calculated in a. and run Newton's Method)

j=0;
Q=1;
for i=1:3%different steps
    L(1)=u(i);
    while Q~=0 % computer tolerance, find the zero, run Newton
        j=j+1;
        L(j+1)=L(j)-(L(j)-u(i)-K(i)*(-2*(2-K(i))*L(j)-L(j)^2))/(1+K(i)*(2*(2-K(i))+L(j)^2));
        Q=L(j+1)-L(j);
    end
    P(i)=L(j+1);%the solutions
    j=0;
end
P%all solutions 

r2=(P(1)-P(2))/(P(2)-P(3))%ratio

% % c. Trapezoidal Method

f = @(x) -6*x-x^2;
g = @(u,x,k) u-x-.5*k*(f(x)+f(u));

Y

% % d. Leapfrog Method



% 2.

% % a.


@(t) f=-e^-t*u(t);



% 4.

% % a. 
%% We know that \[ x'(t)=-3tx-x^2 \] \[ x''(t)=-3(1+tx'(t))-2xx' \]




