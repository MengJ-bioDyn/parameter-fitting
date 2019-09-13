% SDCSB Lecture Series 2017
% Meng Jin

% purpose: use crossvalidation to test whether model is adequate


clear all % clear memory
clc
close all % close all plots

% load target data 

load NFB_simple_data
% load syndeg_simple_data_cv


expY=[T0', Y0']; 


% initiate training/validation set
Kfold_ini =6; 
Kfold = Kfold_ini-1;
% here because of "rolling origin", real fold is Kfold-1
tsn = length(expY);
BN  = fix(tsn/Kfold_ini);
train_set = cell(1,Kfold);
valid_set = cell(1, Kfold);


for j=1:Kfold
    if j<Kfold
        valid_id = [1+(j-1)*BN: j*BN]';
        train_id = setdiff([1:tsn]', valid_id);
              
        valid_set{j} = expY(valid_id,:);
        train_set{j} = expY(train_id,:);
        
    elseif j==Kfold
        
        valid_id = [1+(j-1)*BN: tsn]';
        train_id = setdiff([1:tsn]', valid_id);
        
        valid_set{j} = expY(valid_id,:);
        train_set{j} = expY(train_id,:);
        
    end
end



cross_val_error = zeros(1, Kfold-1);


for j =1:Kfold
% number of trajectories to generate per iteration
YBB = train_set{j};

% initial parameters
a = 10;
b = 1;

parms = zeros(1,2);
parms(1) = a;
parms(2) = b;

% initialize time arrary for calculating error function
T = T0;       % the same as "experimental time"

% initialize Y
Y = nan(size(T));

%%% Metropolis parameters
% number of iterations
metroIter = 400;

% array that chooses which parameters to perturb (0=do not perturb)
vPerturb = zeros(1,2);
vPerturb(1) = 1;
vPerturb(2) = 1;

%strength of perturbation (less than 1)
perturbAmp = 0.2;
% scale of constant parameter noise
perturbAddAmp = 0.01;


x_init = 0; %initial condition

        
        
    % initial energy variable (start with very high energy)
    energy = 1e10;
    
    parmsSave = nan(length(parms),metroIter);
    energySave = nan(metroIter);
  
    % use trainging set
    targetT = YBB(:,1);
    [notused  time_index] = ismember(targetT, T);
    targetY = YBB(:,2);
    
    for iMETRO=1:metroIter

        % perturb the parameters, but keep old parameters and energy
        energyOld = energy;
        parmsOld = parms;
        Y_Old = Y; 
        
        r1 =(rand(1,length(parms))-0.5);
        r2 =(rand(1,length(parms))-0.5);
        parms = abs(parms .* (1+vPerturb.*perturbAmp.*r1) + perturbAddAmp*vPerturb.*r2);


        % test your model  and return the trajectories
        
        sol = ode23(@syndeg_simple,T,x_init,[],parms);    
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

%         % plot results
%         figure(1)
%         subplot(2,1,1)
%         % plot mean trajectory for trial and target systems
%         plot(T,Y,'g-', 'LineWidth', 2);    hold on
%         plot(targetT,targetY,'ro', 'LineWidth', 2);
%         plot(T,bestY,'b-', 'LineWidth', 2); hold off
%         xlim([0 max(T)]);
%     %     hold off
%         box on
%         xlabel('time')
%         ylabel('x')
% 
%         % plot energy from fitting iterations
%         subplot(2,1,2);
%         plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
%         xlabel('iteration')
%         ylabel('log energy')
% 
%         drawnow; 
    
    end
    
    figure;

    % plot trajectory for trial and target systems
    subplot(2,1,1);
    hold on
    plot(T(time_index),targetY,'ro', 'LineWidth', 2);
    plot(T,bestY,'b-', 'LineWidth', 2);
    hold off
    xlim([0 max(T)]);
    box on
    xlabel('time')
    ylabel('x')
     title(['a=' num2str(parms(1)) '; b=' num2str(parms(2)) ])

    % plot energy from fitting iterations
    subplot(2,1,2);
    plot(log10(energySave(1:iMETRO)),'r-', 'LineWidth', 1.5);
    xlabel('iteration')
    ylabel('log energy')
    
    para_distribution(:, j) = best_parms'; 
    
    
    %%%%%% validation step %%%%%%%%
    
    % validate your model use the parameters you just got
    sol2 = ode23(@syndeg_simple,T,x_init,[],best_parms);    
    Y2=deval(sol2,T);
    
    YVV = valid_set{j};
        
    % calculate the error difference between your model "prediction" and
    % real validation data set
    [notused  vtime_index] = ismember(YVV(:,1), T);
    YVV_y = YVV(:,2);
    errork = sum( (Y2(vtime_index)' - YVV_y).^2)./length(YVV_y);
    cross_val_error(j) = errork;
    j
end


% plot cross-validation error
figure;
plot(cross_val_error,'bo-')
mean(cross_val_error)
