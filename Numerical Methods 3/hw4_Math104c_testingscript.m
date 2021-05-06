h=1/17; % timestep
N=17; % number of steps taking, i.e. t is in [0,1]
K1=3;
K2=1;
u=[3; 4; 2]; % initial values
for i=1:N-1
    A=[-K1*u(2,i) -K1*u(1,i) K2; -K1*u(2,i) -K1*u(1,i) K2; K1*u(2,i) K1*u(1,i) -K2]; % Jacobian Matrix
    u(:,i+1)=u(:,i)+h*(A*u(:,i));
end

figure();
x=linspace(0,1,N);
plot(x,u(1,1:N));
hold on
plot(x,u(2,1:N));
hold on
plot(x,u(3,1:N));
title('Numerical solutions with timestep h=1/17');
legend('component 1', 'component 2', 'component 3');
hold off