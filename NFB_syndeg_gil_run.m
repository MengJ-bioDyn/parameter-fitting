close all %close plots
clear %clear data


% initial parameters
T = 50.0;   % time of simulation
dt = 0.05;   % time step
C0 = 5.0;  % negative feedback scale
n0 = 3;     % negative feedback Hill coefficient
alpha = 62.0;  % maximal production rate
g = 2;        % degradation rate

% initial array that passes parameters to simulation
parms(1) = C0;
parms(2) = n0;
parms(3) = alpha;
parms(4) = g;

NN=50; %number of trajectories to generate for each parameter combo

XSAVE = nan(NN,T/dt);
tic
for ii=1:NN
    [time X] = NFB_syndeg_gil(T, dt, parms);
    XSAVE(ii,:) = X(1,:);
end
toc

	% compute statistics
    MEAN = mean(XSAVE,1);
    STD = std(XSAVE,1);
   
    %plots of results
    hold on

    % plot one trajectory
    plot(time,XSAVE(10,:), 'o')        
    xlim([0 max(time)]);
    box on
    xlabel('time')
    ylabel('x')
    
    plot(time, MEAN,time, STD, 'linewidth',1);
    
    MEAN(end)
    STD(end)
    