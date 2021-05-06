function [x_halfstep_final,y_halfstep_final,vx_halfstep_final,vy_halfstep_final,n_new] = gravitational_motion_single_step_RK4_adaptive(ideal_accuracy,T,n,x_init,y_init,vx_init,vy_init)
format long

S1=.9; S2=.4; % safety parameters

max_repetitions=10000;

n_new=n;

for j=1:max_repetitions
    
    n=n_new;
    
    tau=T/n;
    
    % take the wholestep
    [x_wholestep,y_wholestep,vx_wholestep,vy_wholestep]=gravitational_motion_single_step_RK4(T,n,x_init,y_init,vx_init,vy_init);
    r_wholestep=sqrt(x_wholestep^2+y_wholestep^2);
    % vel_wholestep=sqrt(vx_wholestep^2+vy_wholestep^2);
    
    % take the halfstep
    n_halftimestep=2*n; % so that tau halves
    [x_halfstep,y_halfstep,vx_halfstep,vy_halfstep]=gravitational_motion_single_step_RK4(T,n_halftimestep,x_init,y_init,vx_init,vy_init);
    [x_halfstep_final,y_halfstep_final,vx_halfstep_final,vy_halfstep_final]=gravitational_motion_single_step_RK4(T,n_halftimestep,x_halfstep,y_halfstep,vx_halfstep,vy_halfstep);
    r_halfstep_final=sqrt(x_halfstep_final^2+y_halfstep_final^2);
    % vel_halfstep_final=sqrt(vx_halfstep_final^2+vy_halfstep_final^2);
    
%     % Estimated Truncation Error
%     farness=ideal_accuracy*(abs(r_wholestep)+abs(r_halfstep_final))/2;
%     error=r_wholestep-r_halfstep_final;
%     LTE=max(abs(error)./(eps+farness)); % we use eps for double precision accuracy
    
%     % Estimated Truncation Error for velocity
%     desired_LTE_accuracy_vel=ideal_accuracy*(abs(vel_wholestep)+abs(vel_halfstep_final))/2;
%     Current_LTE_estimate_vel=vel_wholestep-vel_halfstep_final;
%     LTE_proportional_change_vel=max(abs(Current_LTE_estimate_vel)./abs(desired_LTE_accuracy_vel));
%     tau_estimate_vel=LTE_proportional_change_vel^(-1/5)*tau;

    % Estimated Truncation Error for position
    desired_LTE_accuracy=ideal_accuracy*(abs(r_wholestep)+abs(r_halfstep_final))/2;
    Current_LTE_estimate=r_wholestep-r_halfstep_final;
    LTE_proportional_change=max(abs(Current_LTE_estimate)./abs(desired_LTE_accuracy));
    tau_estimate=LTE_proportional_change^(-1/5)*tau;
    
    
    
    % get the new tau
%     
%     if S1*tau_estimate>S2*tau
%         tau_new=S2*tau;
%     elseif S1*tau_estimate<tau/S2
%         tau_new=tau/S2;
%     else
%         tau_new=S1*tau_estimate;
%     end
    
    
    tau_new=S1*tau*LTE_proportional_change^(-1/5);
    tau_new=max(tau_new,tau/S2);
    tau_new=min(tau_new,tau*S2);
    n_new=T/tau_new;
    
    % check if the new timestep meets the conditions
    if LTE_proportional_change<1
        return % i.e. EUREKA, WE FOUND IT!
    end
end

% couldn't find it :_(
error('Could not find suitable timestep in 5000 tries');


end

