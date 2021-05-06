function [Traj1,Traj2,Vel1,Vel2] = gravitational_motion_RK4(T,n,x_init,y_init,vx_init,vy_init)
format long
tau=T/n; % time step, as specified
Traj1=zeros(1,n);
Traj2=zeros(1,n);% x and y plot values of position
Vel1=zeros(1,n);
Vel2=zeros(1,n); % x and y plot values of velocity

% the acceleration of x and y
Ax=@(x,y) -x./((x^2+y^2).^(3/2)); 
Ay=@(x,y) -y./((x^2+y^2).^(3/2));

% choose values to make it circular elliptic or parabolic
Traj1(1)=x_init; % x initial position
Traj2(1)=y_init; % y initial position
Vel1(1)=vx_init; % x initial velocity
Vel2(1)=vy_init; % y initial velocity, circle

% Runge-Kutta 4th Order Method
for j=2:n
    k1x=Ax(Traj1(j-1),Traj2(j-1));
    k1y=Ay(Traj1(j-1),Traj2(j-1));
    l1x=Vel1(j-1);
    l1y=Vel2(j-1);
    k2x=Ax(Traj1(j-1)+tau*l1x/2,Traj2(j-1)+tau*l1y/2);
    k2y=Ay(Traj1(j-1)+tau*l1x/2,Traj2(j-1)+tau*l1y/2);
    l2x=(Vel1(j-1)+tau*k1x/2);
    l2y=(Vel2(j-1)+tau*k1y/2);
    k3x=Ax(Traj1(j-1)+tau*l2x/2,Traj2(j-1)+tau*l2y/2);
    k3y=Ay(Traj1(j-1)+tau*l2x/2,Traj2(j-1)+tau*l2y/2);
    l3x=(Vel1(j-1)+tau*k2x/2);
    l3y=(Vel2(j-1)+tau*k2y/2);
    k4x=Ax(Traj1(j-1)+tau*l3x,Traj2(j-1)+tau*l3y);
    k4y=Ay(Traj1(j-1)+tau*l3x,Traj2(j-1)+tau*l3y);
    l4x=(Vel1(j-1)+tau*k3x);
    l4y=(Vel2(j-1)+tau*k3y);
    
    % the approximation step
    Vel1(j)=Vel1(j-1)+tau.*((1/6)*k1x+(1/3)*k2x+(1/3)*k3x+(1/6)*k4x);
    % the approximation step
    Vel2(j)=Vel2(j-1)+tau.*((1/6)*k1y+(1/3)*k2y+(1/3)*k3y+(1/6)*k4y);
    % the approximation step
    Traj1(j)=Traj1(j-1)+tau.*((1/6)*l1x+(1/3)*l2x+(1/3)*l3x+(1/6)*l4x);
    % the approximation step
    Traj2(j)=Traj2(j-1)+tau.*((1/6)*l1y+(1/3)*l2y+(1/3)*l3y+(1/6)*l4y);
end
end

