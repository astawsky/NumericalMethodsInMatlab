function [x_new,y_new,vx_new,vy_new] = gravitational_motion_single_step_RK4(T,n,x_init,y_init,vx_init,vy_init)
format long
tau=T/n; % time step, as specified

% the acceleration of x and y
Ax=@(x,y) -x./((x^2+y^2).^(3/2)); 
Ay=@(x,y) -y./((x^2+y^2).^(3/2));

% one step of the Runge-Kutta 4th Order Method
k1x=Ax(x_init,y_init);
k1y=Ay(x_init,y_init);
l1x=vx_init;
l1y=vy_init;
k2x=Ax(x_init+tau*l1x/2,y_init+tau*l1y/2);
k2y=Ay(x_init+tau*l1x/2,y_init+tau*l1y/2);
l2x=(vx_init+tau*k1x/2);
l2y=(vy_init+tau*k1y/2);
k3x=Ax(x_init+tau*l2x/2,y_init+tau*l2y/2);
k3y=Ay(x_init+tau*l2x/2,y_init+tau*l2y/2);
l3x=(vx_init+tau*k2x/2);
l3y=(vy_init+tau*k2y/2);
k4x=Ax(x_init+tau*l3x,y_init+tau*l3y);
k4y=Ay(x_init+tau*l3x,y_init+tau*l3y);
l4x=(vx_init+tau*k3x);
l4y=(vy_init+tau*k3y);

% the approximation step
vx_new=vx_init+tau.*((1/6)*k1x+(1/3)*k2x+(1/3)*k3x+(1/6)*k4x);
% the approximation step
vy_new=vy_init+tau.*((1/6)*k1y+(1/3)*k2y+(1/3)*k3y+(1/6)*k4y);
% the approximation step
x_new=x_init+tau.*((1/6)*l1x+(1/3)*l2x+(1/3)*l3x+(1/6)*l4x);
% the approximation step
y_new=y_init+tau.*((1/6)*l1y+(1/3)*l2y+(1/3)*l3y+(1/6)*l4y);
end

