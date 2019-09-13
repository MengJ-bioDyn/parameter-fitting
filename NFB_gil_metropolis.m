% SDCSB workshop 2016
% Bart Borek

% purpose: fits Gillespie negative feedback circuit period stats to data

close all %close plots
clc
clear % clear memory

% load target data (statistics that we try to fit)
targetMean=3;
targetStd=0.4;

NN=1; %number of trajectories to generate for each parameter combo

% initial parameters
T = 500.0;   % time of simulation
dt = 0.01;   % time step
C0 = 10.0;  % negative feedback scale
n0 = 4;     % negative feedback Hill coefficient
f = 0.0;   % positive feedback strength
C1 = 30.0;  % positive feedback scale
n1 = 0;     % positive feedback Hill coefficient
alpha = 1000.0;  % maximal production rate
g = 7;        % degradation rate
tauDelay = 1; % time delay

% initial array that passes parameters to simulation
parms(1) = C0;
parms(2) = n0;
parms(3) = f;
parms(4) = C1;
parms(5) = n1;
parms(6) = alpha;
parms(7) = g;
parms(8) = tauDelay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metropolis parameters

metroIter = 100; % number of iterations

% array that chooses which parameters to perturb (0=do not perturb)
vPerturb = zeros(1,length(parms));
vPerturb(6) = 1;
vPerturb(7) = 1;

% strength of perturbation (<1)
perturbAmp = 0.4; 
% scale of constant parameter noise
perturbAddAmp = 0.0;

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
    tic
    for ii=1:NN
        [time X] = NFB_gil(T, dt, parms);
        XSAVE(ii,:) = X(2,:);
    end
    toc
    
	% compute statistics
    [peakAmp,peakTime]=findpeaks(smooth(XSAVE),'minpeakdistance',100,'minpeakheight',10);
    if numel(peakTime)>1
        period=diff(time(peakTime));
        MEAN=mean(period);
        STD=std(period);    
    else
        MEAN=0;
        STD=0;  
    end
        
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
        best_peakAmp=peakAmp;
        best_peakTime=peakTime;
        best_fd=fitdist(period','normal');
    elseif energy<=nanmin(energySave(1:iMETRO))
        best_x=XSAVE(ii,:);
        best_peakAmp=peakAmp;
        best_peakTime=peakTime;
        best_fd=fitdist(period','normal');
    end
    
    % plot results
    clf; % clear figure
    
	% plot trajectory for trial and target systems
    subplot(2,2,1);    
    plot(time,XSAVE,time(peakTime),peakAmp,'ro', 'LineWidth', 1.5)        
    xlim([0 max(time)]);
    box on
    xlabel('time')
    ylabel('x')
    
    % plot period distribution
    subplot(2,2,2);    
    hold on
    fd=fitdist(period','normal');
	X=0:0.05:5;    
    plot(X,normpdf(X,fd.mu,fd.sigma));
    plot(X,normpdf(X,targetMean,targetStd),'r');    
    box on
    xlabel('period')    
    
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
subplot(2,2,1);    
plot(time,best_x,time(best_peakTime),best_peakAmp,'ro', 'LineWidth', 1)        
xlim([0 max(time)]);
box on
xlabel('time')
ylabel('x')
title(['a=' num2str(parms(6)) '; b=' num2str(parms(7)) ])
    
% plot normal distributions with target and best means and stds 
subplot(2,2,2);    
hold on
X=0:0.05:5;    
plot(X,normpdf(X,best_fd.mu,best_fd.sigma), 'LineWidth', 2);
plot(X,normpdf(X,targetMean,targetStd),'rx', 'LineWidth', 2);    
box on
xlabel('period')    
    
% plot energy from fitting iterations
subplot(2,2,[3 4]);
plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
xlabel('iteration')
ylabel('log energy')
    