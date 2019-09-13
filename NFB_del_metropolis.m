% SDCSB workshop 2017
% Meng Jin

% fits a negative feedback oscillator trajectory to data



clear % clear memory
clc
close all %close plots

% load target data (statistics that we try to fit)
load NFB_del_target %[50,10,2,0.5]
targetY=Y0; 


% number of trajectories to generate per iteration
NN = 1;

% initial parameters
a = 50;
b = 10;
n = 2;
delay = 0.5;

% initial array that passes parameters to simulation
parms(1) = a;
parms(2) = b;
parms(3) = n;
parms(4) = delay;

% initialize time arrary for calculating error function
T = T0;  % time is set the same as "experimental time"

% initialize Y
Y = nan(size(T));

% Metropolis parameters

% number of iterations
metroIter = 500;

% array that chooses which parameters to perturb (0=do not perturb)
vPerturb(1) = 1;
vPerturb(2) = 1;
vPerturb(3) = 1;
vPerturb(4) = 1;

% strength of perturbation (less than 1)
perturbAmp = 0.2;
% scale of constant parameter noise
perturbAddAmp = 0.1;

% array to store iterations of parameters
parmsSave = nan(length(parms),metroIter);
energySave = nan(metroIter);

% initial energy variable (start with very high energy)
energy = 1e10;

for iMETRO=1:metroIter
        
    % perturb the parameters, but keep old parameters and energy
    energyOld = energy;
    parmsOld = parms;
    Y_Old = Y;

    r1 =(rand(1,length(parms))-0.5);
    r2 =(rand(1,length(parms))-0.5);
    parms = abs(parms .* (1+vPerturb.*perturbAmp.*r1) + perturbAddAmp*vPerturb.*r2);
  
    % run many simulations and return the trajectories
    x_init = 0; %initial condition
    
    % tic
    del = parms(4);    
    sol = dde23(@NFB_del,del,x_init,T,[],parms);    
    Y=deval(sol,T);
    % toc

    % thermal energy
    kT = 0.01;
    
    % compute energy 
    energy = mean((Y-targetY).^2);   
    
    
    % accept or reject step
    if energy>energyOld
        if (exp((energyOld-energy)/kT) < rand)
            parms = parmsOld;
            energy = energyOld;
            Y = Y_Old;
        end
    end    
    
        %keep best trajectory
    if  iMETRO==1
        bestY=Y;
        best_parms = parms;
    elseif energy<=nanmin(energySave(1:iMETRO))
        bestY=Y;
        best_parms = parms;
    end

	% store metropolis iteration
    energySave(iMETRO) = energy;
    parmsSave(:,iMETRO) = parms;
    
    % plot results
    figure(1)
    subplot(2,1,1)
	% plot trajectory for trial and target systems
    plot(T,Y,'g--', 'LineWidth', 2);    hold on
    plot(T,targetY,'ro', 'LineWidth', 2);
    plot(T,bestY,'b-', 'LineWidth', 2); hold off
    xlim([0 max(T)]);
    box on
    xlabel('time')
    ylabel('x')
    
    % plot energy from fitting iterations
    subplot(2,1,2);
    plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
    xlabel('iteration')
    ylabel('log energy')
    
    
    drawnow;  
    
end


% plot best and target trajectory
subplot(2,1,1);
hold on
plot(T,targetY,'ro', 'LineWidth', 2);
plot(T,bestY,'b-', 'LineWidth', 2);
hold off
xlim([0 max(T)]);
box on
set(gca,'fontsize',16)
xlabel('time')
ylabel('y')
 title(['a=' num2str(parms(1)) '; b=' num2str(parms(2)) '; n=' num2str(parms(3)) '; delay=' num2str(parms(4)) ])

% plot energy from fitting iterations
subplot(2,1,2);
plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
set(gca,'fontsize',16)
xlabel('iteration')
ylabel('log energy')



