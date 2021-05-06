function [rho_matrix,v_matrix] = FTCS_FixedWallsBC...
    (t_end,rho_init,v_init,g,number_of_timesteps,number_of_spacesteps,h)
format long
% this function is for the periodic BCs
% t_end,rho_init,v_init,g,number_of_spacesteps,h
% in order are: the last timestep, initial rho values,
% initial v values, the gravitational
% constant, the number of timesteps, the number of spacesteps, 
% the space stepsize, and the accuracy for controlling fluid mass

% Initial Conditions
% rows are moving up in time and columns is space
rho_matrix=zeros(number_of_timesteps,number_of_spacesteps);
rho_matrix(1,:)=rho_init;

v_matrix=zeros(number_of_timesteps,number_of_spacesteps);
v_matrix(1,:)=v_init;

% Periodic Boundary Conditions
rho_matrix(:,1)=rho_init(1);
v_matrix(:,1)=v_init(1);
rho_matrix(:,number_of_spacesteps)=rho_init(number_of_spacesteps);
v_matrix(:,number_of_spacesteps)=v_init(number_of_spacesteps);

tau=t_end/number_of_timesteps;

for time=2:number_of_timesteps % excluding boundaries
    for space=1:number_of_spacesteps
        if space==1
            rho_matrix(time,space)=rho_matrix(time-1,space)-(tau/(h))*...
                ( rho_matrix(time-1,space+1)*v_matrix(time-1,space+1) );
            v_matrix(time,space)=0;
        end
        if space==number_of_spacesteps
            rho_matrix(time,space)=rho_matrix(time-1,space)+(tau/(h))*...
                ( rho_matrix(time-1,space-1)*v_matrix(time-1,space-1) );
            v_matrix(time,space)=0;
        end
        if space~=number_of_spacesteps && space~=1
            rho_matrix(time,space)=rho_matrix(time-1,space)-(tau/(2*h))*...
                ( rho_matrix(time-1,space+1)*v_matrix(time-1,space+1)-...
                rho_matrix(time-1,space-1)*v_matrix(time-1,space-1) );
        end
        if rho_matrix(time,space)~=0 && space~=number_of_spacesteps && space~=1
            v_matrix(time,space)=.5*(v_matrix(time-1,space+1)+...
                v_matrix(time-1,space-1))*(1-tau/(2*h))*...
                (v_matrix(time-1,space+1)-v_matrix(time-1,space-1))...
                -tau/(rho_matrix(time-1,space)*2*h)*(...
                rho_matrix(time-1,space+1)-rho_matrix(time-1,space-1))...
                +tau*g;
%             v_matrix(time,space)=(.5*(v_matrix(time-1,space+1)+...
%                 v_matrix(time-1,space-1))-(tau/(2*h))*...
%                 (1/rho_matrix(time,space))*(rho_matrix(time-1,space+1)...
%                 -rho_matrix(time-1,space-1)) +g*tau)/...
%                 ( (1+(tau/(2*h))*( v_matrix(time-1,space+1)-...
%                 v_matrix(time-1,space-1) ) ) );
%             v_matrix(time,space)=v_matrix(time-1,space)-tau/(2*h)*(...
%                 v_matrix(time-1,space)*(v_matrix(time-1,space+1)-...
%                 v_matrix(time-1,space-1))+(1/rho_matrix(time,space)...
%                 )*(rho_matrix(time-1,space+1)-rho_matrix(time-1,...
%                 space-1)))+g*tau;
        end
    end
    
    if isnan(rho_matrix(time,space))
        keyboard
    elseif isnan(v_matrix(time,space))
        keyboard
    end

end





end






