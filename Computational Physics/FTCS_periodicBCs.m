function [rho_matrix,v_matrix,plot_arrays] = FTCS_periodicBCs...
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

% Do we need this really?
% k=1;
% P_matrix=k*rho_matrix;

tau=t_end/number_of_timesteps;

photoshots=[ceil(number_of_timesteps/5),ceil(2*number_of_timesteps/5)...
    ceil(3*number_of_timesteps/5) ceil(4*number_of_timesteps/5)...
    ceil(5*number_of_timesteps/5)];

plot_arrays=cell(2,length(photoshots));
index=1;

total_fluid_mass=sum(length(find(rho_init)))/number_of_spacesteps;
FTCS_fluid_mass=total_fluid_mass;

% x is discretized from x_1 to x_N
for time=2:number_of_timesteps % excluding boundaries
    if find(photoshots==time)
        plot_arrays(1,index)={rho_matrix(time,:)};
        plot_arrays(2,index)={v_matrix(time,:)};
        index=index+1;
    end
    for space=2:number_of_spacesteps % excluding initial
        % solve for the rho variable first
        
        % the periodic BCs
        if space==number_of_spacesteps
            
%             if rho_matrix(time-1,space-1)~=0 || rho_matrix(time-1,2)~=0 ...
%                     || rho_matrix(time-1,space)~=0
% %                 rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*...
% %                     (rho_matrix(time-1,2)*v_matrix(time-1,2)...
% %                     +rho_matrix(time-1,space-1)*v_matrix(time-1,space-1));
% %                 rho_matrix(time,1)=rho_matrix(time,space);
%                 rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*(...
%                     rho_matrix(time-1,space)*(v_matrix(time-1,2)+...
%                     v_matrix(time-1,space-1))+v_matrix(time-1,space)*(...
%                     rho_matrix(time-1,2)+rho_matrix(time-1,space-1)));
%             end

                rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*...
                    (rho_matrix(time-1,2)*v_matrix(time-1,2)...
                    +rho_matrix(time-1,space-1)*v_matrix(time-1,space-1));
                rho_matrix(time,1)=rho_matrix(time,space);
%             rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*(...
%                 rho_matrix(time-1,space)*(v_matrix(time-1,2)+...
%                 v_matrix(time-1,space-1))+v_matrix(time-1,space)*(...
%                 rho_matrix(time-1,2)+rho_matrix(time-1,space-1)));
            
        else
            rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*...
                    (rho_matrix(time-1,space+1)*v_matrix(time-1,space+1)...
                    +rho_matrix(time-1,space-1)*v_matrix(time-1,space-1));
%             if rho_matrix(time-1,space-1)~=0 || rho_matrix(time-1,space+1)~=0 ...
%                     || rho_matrix(time-1,space)~=0
%                 rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*...
%                     (rho_matrix(time-1,space+1)*v_matrix(time-1,space+1)...
%                     +rho_matrix(time-1,space-1)*v_matrix(time-1,space-1));
%             end


%             rho_matrix(time,space)=rho_matrix(time-1,space)-tau/(2*h)*(...
%                 rho_matrix(time-1,space)*(v_matrix(time-1,space+1)+...
%                 v_matrix(time-1,space-1))+v_matrix(time-1,space)*(...
%                 rho_matrix(time-1,space+1)+rho_matrix(time-1,space-1)));
        end
        
        % solve for the v variable second
        % Periodic BCs
        if rho_matrix(time,space)~=0
            if space==number_of_spacesteps
%                 && ...
%                     (rho_matrix(time-1,space-1)~=0 || rho_matrix(time-1,2)~=0)
                v_matrix(time,space)=v_matrix(time-1,space)-tau/(2*h)*(...
                    v_matrix(time-1,space)*(v_matrix(time-1,2)+...
                    v_matrix(time-1,space-1))+(1/rho_matrix(time,space)...
                    )*(rho_matrix(time-1,2)+rho_matrix(time-1,space-1)))+g*tau;
                v_matrix(time,1)=v_matrix(time,space);
            else
%                 if (rho_matrix(time-1,space-1)~=0 || rho_matrix(time-1,space+1)~=0)
                v_matrix(time,space)=v_matrix(time-1,space)-tau/(2*h)*(...
                    v_matrix(time-1,space)*(v_matrix(time-1,space+1)+...
                    v_matrix(time-1,space-1))+(1/rho_matrix(time,space)...
                    )*(rho_matrix(time-1,space+1)+rho_matrix(time-1,space-1)))+g*tau;
            end
            
        end
        
        if isnan(rho_matrix(time,space))
            keyboard
        elseif isnan(v_matrix(time,space))
            keyboard
        end
    end
end

end

