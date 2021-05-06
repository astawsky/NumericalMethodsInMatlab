function [x_halfstep_final,y_halfstep_final,vx_halfstep_final,vy_halfstep_final,n_new] = gravitational_motion_single_step_RK4_adaptive_2(ideal_accuracy,T,n,x_init,y_init,vx_init,vy_init)
format long

S1=.9; S2=.4; % safety parameters

max_repetitions=10000;

n_new=n;

for j=1:max_repetitions
    
    n=n_new;
    
    tau=T/n;
    
    % take the wholestep
    [x_wholestep,y_wholestep,vx_wholestep,vy_wholestep]=gravitational_motion_single_step_RK4(T,n,x_init,y_init,vx_init,vy_init);
%     r_wholestep=sqrt(x_wholestep^2+y_wholestep^2);
%     vel_wholestep=sqrt(vx_wholestep^2+vy_wholestep^2);
    
    % take the halfstep
    n_halftimestep=2*n; % so that tau halves
    [x_halfstep,y_halfstep,vx_halfstep,vy_halfstep]=gravitational_motion_single_step_RK4(T,n_halftimestep,x_init,y_init,vx_init,vy_init);
    [x_halfstep_final,y_halfstep_final,vx_halfstep_final,vy_halfstep_final]=gravitational_motion_single_step_RK4(T,n_halftimestep,x_halfstep,y_halfstep,vx_halfstep,vy_halfstep);
%     r_halfstep_final=sqrt(x_halfstep_final^2+y_halfstep_final^2);
%     vel_halfstep_final=sqrt(vx_halfstep_final^2+vy_halfstep_final^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% another try %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Estimated Truncation Error for x position
    desired_LTE_accuracy_x=ideal_accuracy*(abs(x_wholestep)+abs(x_halfstep_final))/2;
    Current_LTE_estimate_x=x_wholestep-x_halfstep_final;
    LTE_proportional_change_x=max(abs(Current_LTE_estimate_x)./(desired_LTE_accuracy_x+eps));
    % tau_estimate_x=LTE_proportional_change_x^(-1/5)*tau;
    
    % Estimated Truncation Error for y position
    desired_LTE_accuracy_y=ideal_accuracy*(abs(y_wholestep)+abs(y_halfstep_final))/2;
    Current_LTE_estimate_y=y_wholestep-y_halfstep_final;
    LTE_proportional_change_y=max(abs(Current_LTE_estimate_y)./(desired_LTE_accuracy_y+eps));
    % tau_estimate_y=LTE_proportional_change_y^(-1/5)*tau;
    
    % Estimated Truncation Error for x velocity
    desired_LTE_accuracy_vx=ideal_accuracy*(abs(vx_wholestep)+abs(vx_halfstep_final))/2;
    Current_LTE_estimate_vx=vx_wholestep-vx_halfstep_final;
    LTE_proportional_change_vx=max(abs(Current_LTE_estimate_vx)./(desired_LTE_accuracy_vx+eps));
    % tau_estimate_vx=LTE_proportional_change_vx^(-1/5)*tau;
    
    % Estimated Truncation Error for y velocity
    desired_LTE_accuracy_vy=ideal_accuracy*(abs(vy_wholestep)+abs(vy_halfstep_final+eps))/2;
    Current_LTE_estimate_vy=vy_wholestep-vy_halfstep_final;
    LTE_proportional_change_vy=max(abs(Current_LTE_estimate_vy)./(desired_LTE_accuracy_vy));
    % tau_estimate_vy=LTE_proportional_change_vy^(-1/5)*tau;
    
    tau_new_x=S1*tau*LTE_proportional_change_x^(-1/5);
    tau_new_x=max(tau_new_x,tau/S2);
    tau_new_x=min(tau_new_x,tau*S2);
    
    tau_new_y=S1*tau*LTE_proportional_change_y^(-1/5);
    tau_new_y=max(tau_new_y,tau/S2);
    tau_new_y=min(tau_new_y,tau*S2);
    
    tau_new_vx=S1*tau*LTE_proportional_change_vx^(-1/5);
    tau_new_vx=max(tau_new_vx,tau/S2);
    tau_new_vx=min(tau_new_vx,tau*S2);
    
    tau_new_vy=S1*tau*LTE_proportional_change_vy^(-1/5);
    tau_new_vy=max(tau_new_vy,tau/S2);
    tau_new_vy=min(tau_new_vy,tau*S2);
    
    tau_new=min([tau_new_x,tau_new_vx,tau_new_y,tau_new_vy]);
    
    n_new=T/tau_new;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     % Estimated Truncation Error for velocity
%     desired_LTE_accuracy_vel=ideal_accuracy*(abs(vel_wholestep)+abs(vel_halfstep_final))/2;
%     Current_LTE_estimate_vel=vel_wholestep-vel_halfstep_final;
%     LTE_proportional_change_vel=max(abs(Current_LTE_estimate_vel)./abs(desired_LTE_accuracy_vel));
%     tau_estimate_vel=LTE_proportional_change_vel^(-1/5)*tau;
% 
%     % Estimated Truncation Error for position
%     desired_LTE_accuracy=ideal_accuracy*(abs(r_wholestep)+abs(r_halfstep_final))/2;
%     Current_LTE_estimate=r_wholestep-r_halfstep_final;
%     LTE_proportional_change=max(abs(Current_LTE_estimate)./abs(desired_LTE_accuracy));
%     tau_estimate=LTE_proportional_change^(-1/5)*tau;
    
    
    
    % get the new tau
%     
%     if S1*tau_estimate>S2*tau
%         tau_new=S2*tau;
%     elseif S1*tau_estimate<tau/S2
%         tau_new=tau/S2;
%     else
%         tau_new=S1*tau_estimate;
%     end
    
    
%     tau_new=S1*tau*LTE_proportional_change^(-1/5);
%     tau_new=max(tau_new,tau/S2);
%     tau_new=min(tau_new,tau*S2);
%     
%     tau_new_vel=S1*tau*LTE_proportional_change_vel^(-1/5);
%     tau_new_vel=max(tau_new_vel,tau/S2);
%     tau_new_vel=min(tau_new_vel,tau*S2);
%     
%     tau_new=min(tau_new,tau_new_vel);
%     
%     n_new=T/tau_new;
    
    % check if the new timestep meets the conditions
    if LTE_proportional_change_x<1 && LTE_proportional_change_vx<1 && LTE_proportional_change_y<1 && LTE_proportional_change_vy<1
        return % i.e. EUREKA, WE FOUND IT!
    end
end

% couldn't find it :_(
error('Could not find suitable timestep in 10000 tries');


end

