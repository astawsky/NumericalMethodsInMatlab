% Alejandro Stawsky, Math 104C, Spring 2018, HW 3

format long

% 4. (Problem 3 and 4? Did he mean 2 and 3?)

% % a. Problem 2 Method

k=.01; % length of step
T=1; % Distance
N=T/k; % Number of steps
U(1)=1; % initial condition
U(2)=exp(k); % initial condition

for i=2:N % method from problem 2
    U(i+1)=(3/2)*U(i)-(1/2)*U(i-1)+k*((5/4)*U(i)-(3/4)*U(i-1));
end
x=0:k:1;
figure();
plot(x,U,'b--o');
hold on
plot(x,exp(x),'r');
legend('Approximation','True Solution');
hold off

% % b. Problem 3 Method

for i=2:N % method from problem 3
    U(i+1)=(3)*U(i)-(2)*U(i-1)+k*((1/2)*U(i)-(3/2)*U(i-1));
end
x=0:k:1;
figure();
plot(x,U,'b--o');
hold on
plot(x,exp(x),'r');
legend('Approximation','True Solution');
hold off

% % c. 
% We can see that this method diverges and is worse than the previous one
% (zero-stable?)



