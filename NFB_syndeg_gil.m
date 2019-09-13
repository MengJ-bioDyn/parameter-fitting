% SDCSB Lecture Series 2013
% Will Mather
% purpose: simulates a one component feedback system with delay in the feedback
%
% returns:
%		Time: array of time values
%		XTraj: array of system variables at each of these times
% accepts:
%		Ttot: total time to run the simulation
%		tSamp: the sampling time for output (smaller leads to more frequent output)
%		parms: an array of parameters to be used
%

function [Time XTraj] = NFB_syndeg_gil(Ttot,tSamp,parms)

% time variable for simulation
t = 0;

% output iteration variable
nIter = 0;

% variables
iR = 1;
% iR1 = 2;
% iR2 = 3;
% iR3 = 4;
% iR4 = 5;
% iR5 = 6;
x(iR) = 0;
% x(iR1) = 0;
% x(iR2) = 0;
% x(iR3) = 0;
% x(iR4) = 0;
% x(iR5) = 0;

% trajectory
Time = zeros(1,ceil(Ttot/tSamp));
XTraj = zeros(length(x),ceil(Ttot/tSamp));

% parameters
C0 = parms(1);			% scales when repression occurs
n0 = parms(2);			% Hill coefficient for repression
% f = parms(3)+1;			% strength of positive feedback (greater than zero)
% C1 = parms(4);			% scales when activation occurs
% n1 = parms(5);			% Hill coefficient for activation
alpha = parms(3);		% maximal production rate
g = parms(4);			% degradation rate of proteins
% tauDelay = parms(8);	% delay for feedback
% finv = 1/f;

% rate for delay chain
mDelay = 5;   % number of intermediate steps (keep fixed)
% delayRate = mDelay/tauDelay;

while (t<Ttot)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute rates for reactions
    
    % reaction 1: birth
     Kv(1) = alpha / (1+(x(iR)/C0)^n0);
%      Kv(1) = alpha;

    % reaction 2: death
    Kv(2) = g*x(iR);
    
%     % reactions 3 through 3+mDelay-1 are the delay chain
%     Kv(3) = x(iR1)*delayRate;
%     Kv(4) = x(iR2)*delayRate;
%     Kv(5) = x(iR3)*delayRate;
%     Kv(6) = x(iR4)*delayRate;
%     Kv(7) = x(iR5)*delayRate;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % number of reactions
    NDim = length(Kv);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute cumulative vector
    Kvcum = cumsum(Kv);
    KvTot = Kvcum(NDim);
    Kvcum = Kvcum ./ KvTot;

    % find time of next reaction
    tau = -(1/KvTot)*log(rand);
    t = t+tau;

    % output trajectory when time exceeds sampling time
    while ((t>nIter*tSamp) && (nIter*tSamp < Ttot))
        nIter = nIter+1;
        Time(nIter) = nIter*tSamp;
        XTraj(:,nIter) = x;
    end

    % find out which reaction to choose
    RR = rand;
    nrxn = 1;
    while (Kvcum(nrxn)<RR)
        nrxn = nrxn+1;
    end

    switch (nrxn)
        
        case (1)
            % reaction 1: birth (at front of delay chain)
%             x(iR1) = x(iR1)+1;
            x(iR) = x(iR)+1;
            
        case (2)
            % reaction 2: death
            x(iR) = x(iR)-1;
            
            
%         % delay chain
%         case (3)
%             x(iR1) = x(iR1)-1;
%             x(iR2) = x(iR2)+1;
%         case (4)
%             x(iR2) = x(iR2)-1;
%             x(iR3) = x(iR3)+1;
%         case (5)
%             x(iR3) = x(iR3)-1;
%             x(iR4) = x(iR4)+1;
%         case (6)
%             x(iR4) = x(iR4)-1;
%             x(iR5) = x(iR5)+1;
%         case (7)
%             x(iR5) = x(iR5)-1;
%             x(iR) = x(iR)+1;
    end
    
    
end

end

