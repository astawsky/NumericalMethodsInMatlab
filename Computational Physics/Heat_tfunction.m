function [Time_array,charac_time] = Heat_tfunction(t, h, tau, number_of_steps, K, L)

% Nt is number of time steps
% K is one
% h is the adaptive scaling factor
% tau
% L is length equal to one

Nx = ceil(L/h + 1);
Time_array = zeros(Nx, number_of_steps); % Time

% don't really need isa
if isa(t, 'double') % Check whether T_0 is a double (-precision array)
    if (length(t) == number_of_steps || length(t) == 1)
        Time_array(1,:) = t;
    end
elseif isa(t, 'function_handle')
    Time_array(1,:) = t(linspace(0,tau*number_of_steps,size(Time_array,2)));
end

for n = 1:(number_of_steps-1)
    for i = 2:(Nx-1)
        Time_array(i,n+1) = Time_array(i,n) + K *tau/h^2*(Time_array(i-1,n)-2*Time_array(i,n)+Time_array(i+1,n));
        if i > (n+1)
            break
        end
    end
    %Accuracy verification
    if abs(Time_array(round(Nx*1/2),n+1)- (Time_array(1,1) - L/K*h*1/2*(Nx-1))) < 10^-4
        charac_time = tau*(n+1);
        Time_array(:,n+2:end) = [];
        break
    end
end

return