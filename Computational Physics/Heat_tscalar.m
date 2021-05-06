function [Time_matrix,charac_time] = Heat_tscalar(t, h, tau, number_of_timesteps, K, L)
% t is a scalar
% Nt is number of time steps
% K is one
% h is the adaptive scaling factor
% tau
% L is length equal to one

number_of_spacesteps=ceil(L/h + 1);
Time_matrix=zeros(number_of_spacesteps, number_of_timesteps); % Time
Time_matrix(1,:)=t;
desired_accuracy=10^(-4);

for n = 1:(number_of_timesteps-1)
    
    % the numerical steps for space param, x
    for i = 2:(number_of_spacesteps-1)
        Time_matrix(i,n+1)=Time_matrix(i,n)+(K*tau*(Time_matrix(i-1,n)-2*Time_matrix(i,n)+Time_matrix(i+1,n)))/h^2;
    end
    
    % figure out the characteristic time
    if abs(Time_matrix(ceil(number_of_spacesteps/2),n+1)- (Time_matrix(1,1) - L*h*(number_of_spacesteps-1)/(K*2))) < desired_accuracy
        charac_time = tau*(n+1);
        Time_matrix(:,n+2:end) = [];
        break
    end
end