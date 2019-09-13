%SDCSB workshop 2016
%Bart Borek

clear, clc

%generate target dataset %[5 1 0.5]

%% Specify Parameters
a1 = 5;   % unbound production
b1 = 1;   % degradation
k1=0.5; % repression strength

N=1; %# of species
T = 0:0.2:5; %time of simulation vector
Y=zeros(N,length(T));
Y2=zeros(N,length(T));

    parameters = [a1,b1,k1]; % create parameter vector    
    x_init= 0; %initial value of x
    sol = ode23(@NFB_simple,T,x_init,[],parameters);    
    solution=deval(sol,T);
    noise=0.2*(rand(1,length(T))-0.1);
    Y(:)=solution(1,:) +noise;    

plot(T,Y,'ro') %,'LineWidth',2
xlabel('time','FontSize',12)
ylabel('X','FontSize',12)
set(gca,'FontSize',12)

hold on


%%%%%
%% parameter guess simulation
a1 = 10;   % unbound production
b1 = 0.6;   % degradation
k1=1.5;   % repression strength


    parameters = [a1,b1,k1]; % create parameter vector    
    x_init= 0; %initial value of x
    sol = ode23(@NFB_simple,T,x_init,[],parameters);    
    solution=deval(sol,T);
    Y2(:)=solution(1,:)';    
    
    
    plot(T,Y2,'o-') %,'LineWidth',2
%legend('n=1','n=2','Location','se')
xlabel('time','FontSize',12)
ylabel('X','FontSize',12)
set(gca,'FontSize',12)

%calculate and show error value for this instance
E =sum((Y-Y2).^2)




