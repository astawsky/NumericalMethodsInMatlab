function [A_mean,B_mean,C_mean,A_saves,B_saves,C_saves,...
    time_saves,time_mean]=stochastic_simulations

k1=100; k2=1000; k3=0.01; k4=9900; k6 = 100;
k5=1000;% sec-1, for t<10 min
k5_halved=500;% sec-1, for t>=10 min

A_init=0;
B_init=0;
C_init=0;

numberofrealisations=10;
run_time = 60*40;
change_k5_time = 60*10;

A_saves = zeros(10^7,numberofrealisations);
B_saves = zeros(10^7,numberofrealisations);
C_saves = zeros(10^7,numberofrealisations);
time_saves = zeros(10^7,numberofrealisations);

for n=1:numberofrealisations
    
    k5=1000;
    A_init=0;
    B_init=0;
    C_init=0;
    IC=[A_init, B_init, C_init];
    time=0;
    k=1;
    A_saves(1,n)=IC(1);
    B_saves(1,n)=IC(2);
    C_saves(1,n)=IC(3);
    time_saves(1,n)=0;
    
    while (time<run_time)
        if time>=change_k5_time
            k5=k5_halved;
        end
        rr=rand(2,1);
        propensities=[k1; % C is created 
                      IC(3)*k2; % C is degraded and turned into B
                      IC(2)*k3; % B is degraded
                      IC(1)*IC(3)*k4; % C is degraded and A stays the same
                      k5; % A is created
                      IC(1)*k6];% A is degraded
        Propensities = sum(propensities);
        cumProp = cumsum(propensities/Propensities);
         
        reaction_index = find(cumProp>=rr(2),1,'first');
         
        if numel(reaction_index)~=1
            break
        end
        switch reaction_index
            case 1
                IC(3)=IC(3)+1;
            case 2
                IC(2)=IC(2)+1;
                IC(3)=IC(3)-1;
            case 3
                IC(2)=IC(2)-1;
            case 4
                IC(3)=IC(3)-1;
            case 5
                IC(1)=IC(1)+1;
            case 6
                IC(1)=IC(1)-1;
        end
              
        time=time+(1/Propensities)*log(1/rr(1));
        k=k+1;
        A_saves(k,n)=IC(1);
        B_saves(k,n)=IC(2);
        C_saves(k,n)=IC(3);
        time_saves(k,n)=time;
    end
 
end


for n=1:numberofrealisations    
   ind = find(time_saves(:,n)>0,1,'last');
   time_saves(ind+1:end,n) = time_saves(ind,n);
   A_saves(ind+1:end,n) = A_saves(ind,n);
   B_saves(ind+1:end,n) = B_saves(ind,n);
   C_saves(ind+1:end,n) = C_saves(ind,n);
end

% to plot the means
A_mean = mean(A_saves,2);
B_mean = mean(B_saves,2);
C_mean = mean(C_saves,2);
time_mean = mean(time_saves,2);

end