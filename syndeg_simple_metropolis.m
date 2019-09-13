% SDCSB Lecture Series 2017
% Meng Jin

% purpose: fits a synthesis/degradation model to data


clear % clear memory
clc

close all % close all plots

% load target data 
load syndeg_simple_data 
targetY=Y0; 


% initial parameters
a = 10;
b = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial array that passes parameters to simulation
parms(1) = a;
parms(2) = b;

% initiate time array, the same as in the loaded data

T = T0; %time

% initiate Y

Y = nan(size(T));

%%% Metropolis parameters
% number of iterations
metroIter = 500;

% array that chooses which parameters to perturb (0=do not perturb)
vPerturb(1) = 1;
vPerturb(2) = 1;

%strength of perturbation (less than 1)
perturbAmp = 0.2;
% scale of constant parameter noise
perturbAddAmp = 0.01;
%
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
  

    % run simulations and return the trajectories
    x_init = 0; %initial condition
    
    % time array T has already been defined, the same as in loaded data

    tic
    sol = ode23(@syndeg_simple,T,x_init,[],parms);    
    Y=deval(sol,T);
    toc

    % thermal energy
    kT = 0.01;

    % compute energy 
    energy = mean((Y-targetY).^2);     

    % accept or reject step
    if energy>energyOld
        if (exp((energyOld-energy)/kT) < rand)
            parms = parmsOld;
            energy = energyOld;
            Y= Y_Old;
        end
    end    

	% store metropolis iteration
    energySave(iMETRO) = energy;
    parmsSave(:,iMETRO) = parms;

    %keep best trajectory
    if  iMETRO==1
        bestY=Y;
        best_parms = parms;
    elseif energy<=nanmin(energySave(1:iMETRO))
        bestY=Y;
        best_parms = parms;
    end
    
    % plot results
    figure(1)
    subplot(2,1,1)
	% plot mean trajectory for trial and target systems
    plot(T,Y,'g-', 'LineWidth', 2);    hold on
    plot(T,targetY,'ro', 'LineWidth', 2);
    plot(T,bestY,'b-', 'LineWidth', 2); hold off
    xlim([0 max(T)]);
%     hold off
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


%show recorded samples and best fit

figure;

% plot trajectory for trial and target systems
subplot(2,1,1);
hold on
plot(T,targetY,'ro', 'LineWidth', 2);
plot(T,bestY,'b-', 'LineWidth', 2);
hold off
xlim([0 max(T)]);
box on
xlabel('time')
ylabel('x')
 title(['a=' num2str(best_parms(1)) '; b=' num2str(best_parms(2)) ])

% plot energy from fitting iterations
subplot(2,1,2);
plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
xlabel('iteration')
ylabel('log energy')
