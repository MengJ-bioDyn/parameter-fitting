%SDCSB workshop 2016
%Bart Borek

% clear

%generate target dataset %[5 1 0.5]

%% Specify Parameters
a1 = 4;   % unbound production
b1 = 2;   % degradation

N=1; %# of species
T = 0:0.2:5; %time of simulation vector
Y0=zeros(N,length(T));
Y2=zeros(N,length(T));

    parameters = [a1,b1]; % create parameter vector    
    x_init= 0; %initial value of x
    sol = ode23(@syndeg_simple,T,x_init,[],parameters);    
    solution=deval(sol,T);
    noise=0.2*(rand(1,length(T))-0.5);
%     noise=[0.0,0.121566234062458,-0.0538864902329159,0.0777389450887783,-0.164406656094219,0.103023044009804,-0.234083576811290,-0.111538507519555,-0.226914304684423,-0.201434109382076,0.161728914163646,0.0974143114879085,-0.0914502599695697,0.225111024419177,-0.232776959748546,-0.0306278201718009,-0.0592207714534958,0.132758394074501,0.147599950568532,-0.156563697722811,-0.00511780210588447,-0.0272068996445503,0.0731565050556323,0.104682415429036,0.127343340991180,-0.111987461500711];
    Y0(:)=solution(1,:) +noise;    

plot(T,Y0,'r+') %,'LineWidth',2
xlabel('time','FontSize',14)
ylabel('X','FontSize',14)
set(gca,'FontSize',14)

hold on
% plot(T,solution(1,:))

%%%%%
%% parameter guess simulation
a1 = 4;   % unbound production
b1 = 1;   % degradation


    parameters = [a1,b1]; % create parameter vector    
    x_init= 0; %initial value of x
    sol = ode23(@syndeg_simple,T,x_init,[],parameters);    
    solution=deval(sol,T);
    Y2(:)=solution(1,:)';    
    
    
    plot(T,Y2,'o-') %,'LineWidth',2
%legend('n=1','n=2','Location','se')
xlabel('time','FontSize',12)
ylabel('X','FontSize',12)
set(gca,'FontSize',12)

%calculate and show error value for this instance
E =sum((Y0-Y2).^2)




