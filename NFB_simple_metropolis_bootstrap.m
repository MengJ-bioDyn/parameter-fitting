% SDCSB Lecture Series 2017
% Meng Jin

% purpose: use bootstrap method to estimate confidence interval of
% parameters (fitted by metroplis method)


clear all
clc

close all % close all plots

% load target data (statistics that we try to fit)
load NFB_simple_data

expY=[T0', Y0']; 

% Number of parameters we want to know confidence interval
para_size =3;

% number of bootstrap sampling
samplingN = 100;

% allocate space for parameter distribution
para_distribution = zeros(para_size, samplingN);

% bootstrap loop. use parallel computing

parfor j =1:samplingN

% generate bootstrap sampling

YBB = genBootstrap(expY);


% initial parameters
a = 10;
b = 1;
k = 1;
parms = zeros(1,3);
parms(1) = a;
parms(2) = b;
parms(3) = k;

% initialize time arrary for calculating error function
T = T0;               % the same as "experimental time"

% initialize Y
Y = nan(size(T));


%%% Metropolis parameters
% number of iterations
metroIter = 300;

% array that chooses which parameters to perturb (0=do not perturb)
vPerturb = zeros(1,3);
vPerturb(1) = 1;
vPerturb(2) = 0;
vPerturb(3) = 1;

%strength of perturbation (less than 1)
perturbAmp = 0.2;
% scale of constant parameter noise
perturbAddAmp = 0.01;

% array to store iterations of parameters
parmsSave = nan(length(parms),metroIter);
energySave = nan(metroIter);


% initial energy variable (start with very high energy)
energy = 1e10;

parmsSave = nan(length(parms),metroIter);
energySave = nan(metroIter);
    
    
    targetT = YBB(:,1);
    targetY = YBB(:,2);
    [notused  time_index] = ismember(targetT, T);

    
    for iMETRO=1:metroIter

        % perturb the parameters, but keep old parameters and energy
        energyOld = energy;
        parmsOld = parms;
        Y_Old = Y; 

        r1 =(rand(1,length(parms))-0.5);
        r2 =(rand(1,length(parms))-0.5);
        parms = abs(parms .* (1+vPerturb.*perturbAmp.*r1) + perturbAddAmp*vPerturb.*r2);


        % run simulations of your model and return the trajectories
        
        %initial condition
        x_init = 0; 
        sol = ode23(@NFB_simple,T,x_init,[],parms);    
        Y=deval(sol,T);

        % thermal energy
        kT = 0.01;

        % compute energy 
        energy = mean((Y(time_index)'-targetY).^2);     

        % accept or reject step
        if energy>energyOld
            if (exp((energyOld-energy)/kT) < rand)
                parms = parmsOld;
                energy = energyOld;
                Y = Y_Old;
            end
        end    

        % store metropolis iteration
        energySave(iMETRO) = energy;
        parmsSave(:,iMETRO) = parms';

        %keep best trajectory
        if  iMETRO==1
            bestY=Y;
            best_parms = parms;
        elseif energy<=nanmin(energySave(1:iMETRO))
            bestY=Y;
            best_parms = parms;
        end

    
    end
   
    
    para_distribution(:, j) = parms'; 
    j
end




    % get the confidence interval
    disp('bootstrap confidence intervals of parameters a, k ')
    CI=prctile(para_distribution,[5 95],2)
    mean_parms = mean(para_distribution,2);
    
    % plot the distribution
    figure;
    subplot(2,1,1);
    hist(para_distribution(1, :),20); set(gca,'fontsize',14); 
    title(['distribution of parameter a, mean= ',num2str(mean_parms(1)),' 90% CI=[',num2str(CI(1,1)),',',num2str(CI(1,2)),']'])
    subplot(2,1,2);
    hist(para_distribution(3, :),20); set(gca,'fontsize',14); 
    title(['distribution of parameter k, mean= ',num2str(mean_parms(3)),' 90% CI=[',num2str(CI(3,1)),',',num2str(CI(3,2)),']'])
    

    
