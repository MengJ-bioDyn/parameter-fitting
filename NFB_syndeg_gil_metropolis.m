% SDCSB workshop 2017
% Meng Jin

% purpose: fits statistic of steady state to Gillespie negative feedback model 

close all %close plots
clc
clear % clear memory

% load target data (statistics that we try to fit)
targetMean=7.6;
targetStd=1.6;

NN=50;
%number of trajectories to generate for each parameter combo

% initial parameters
T = 50.0;   % time of simulation
dt = 0.05;   % time step

C0 = 5;  % negative feedback scale
n0 = 3;     % negative feedback Hill coefficient

alpha = 5.0;  % maximal production rate
g = 7;        % degradation rate

% initial array that passes parameters to simulation
parms(1) = C0;
parms(2) = n0;
parms(3) = alpha;
parms(4) = g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metropolis parameters

metroIter = 200; % number of iterations

% array that chooses which parameters to perturb (0=do not perturb)
vPerturb = zeros(1,length(parms));
vPerturb(1) = 0;
vPerturb(2) = 0;
vPerturb(3) = 1;
vPerturb(4) = 1;

% strength of perturbation (<1)
perturbAmp = 0.3; 
% scale of constant parameter noise
perturbAddAmp = 0.01;

% array to store iterations of parameters
parmsSave = zeros(length(parms),metroIter);
energySave = zeros(1,metroIter);

% initial energy variable (start with very high energy)
energy = 1e10;

% thermal energy
kT = 0.01;

%%%%
% metropolis algorithm

for iMETRO=1:metroIter
        
    % perturb the parameters, but keep old parameters and energy
    energyOld = energy;
    parmsOld = parms;

    r1 =(rand(1,length(parms))-0.5);
    r2 =(rand(1,length(parms))-0.5);
    parms = abs(parms .* (1+vPerturb.*perturbAmp.*r1) + perturbAddAmp*vPerturb.*r2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run many simulations and return the trajectories
    XSAVE = nan(NN,T/dt);

    % tic
    for ii=1:NN 
        [time X] = NFB_syndeg_gil(T, dt, parms);
        XSAVE(ii,:) = X(1,:);
    end
    % toc
    
	% compute statistics    
    MEAN = mean(XSAVE,1); 
    STD = std(XSAVE,1);
        
    % compute energy (currently mean square difference of mean from target mean)
    energy = mean((MEAN-targetMean).^2) + 1*mean((STD-targetStd).^2);
    
    % accept or reject step
    if (energy>energyOld)        
        if (exp((energyOld-energy)/kT) < rand)
            parms = parmsOld;
            energy = energyOld;
        end
    end
 
	% store metropolis iteration
    energySave(iMETRO) = energy;
    parmsSave(:,iMETRO) = parms;
    
    %save run with lowest error
    if  iMETRO==1
        best_x=XSAVE(ii,:);
        best_mean=MEAN;
        best_std=STD;
    elseif energy<=nanmin(energySave(1:iMETRO))
        best_x=XSAVE(ii,:);
        best_mean=MEAN;
        best_std=STD;
    end
    
    % plot results
    clf; % clear figure
    
	% plot trajectory for trial and target systems
    subplot(2,2,[1 2]);    
    plot(time,XSAVE,'o')        
    xlim([0 max(time)]);
    box on
    xlabel('time')
    ylabel('x');
    plot(time, best_mean,'b-', 'LineWidth', 1.5)
       
    
    % plot energy from fitting iterations
    subplot(2,2,[3 4]);
    plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
    xlabel('iteration')
    ylabel('log energy')
    
    drawnow;
    
end

%plotting at the end
figure()

% plot best and target trajectory
subplot(2,2,[1 2]);    
plot(time,best_x,'o',time,best_mean,'b-', 'LineWidth', 1)        
xlim([0 max(time)]);
box on
xlabel('time')
ylabel('x')
title(['alpha=' num2str(parms(3)) '; g=' num2str(parms(4)) '; mean= ',num2str(best_mean(end)) '; std= ',num2str(best_std(end)) ])
 
    
% plot energy from fitting iterations
subplot(2,2,[3 4]);
plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
xlabel('iteration')
ylabel('log energy')
    