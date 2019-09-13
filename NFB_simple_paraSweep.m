%SDCSB workshop 2016
%Bart Borek

clear %clear data
close all %close plots


load NFB_simple_data  

%set up parameter space to step through
nump1=20;
nump2=20;
a1_vec=linspace(0.1,10,nump2);
k1_vec=linspace(0.01,1,nump1);

N=1;%# of species

T= T0;
Y2=zeros(N,length(T));



for i=1:nump1
    for j=1:nump2
        
%assign parameter value for step
a1 = a1_vec(i);   % unbound production
b1 = 1;   % degradation
k1=k1_vec(j); % repression strength


    parameters = [a1,b1,k1]; % create parameter vector    
    x_init= 0; %initial value of x
    sol = ode23(@NFB_simple,T,x_init,[],parameters);    
    solution=deval(sol,T);
    Y2(:)=solution(1,:);    
    

E(i,j)=sum((Y0-Y2).^2);



    end
end

surf(a1_vec,k1_vec, E);
zlabel('E','FontSize',12)
xlabel('a','FontSize',12)
ylabel('k','FontSize',12)
set(gca,'FontSize',12)

%display lowest energy 
[r,c]=find(E==min(min(E)));
LowestE=min(min(E))
LowestaforE=a1_vec(r)
LowestkforE=k1_vec(c)

%plot best parameter trajectory
    best_parameters = [a1_vec(r),b1,k1_vec(c)]; % create parameter vector    
    x_init= 0; %initial value of x
    sol = ode23(@NFB_simple,T,x_init,[],best_parameters);    
    solution=deval(sol,T);
    Y3=solution(1,:);  

figure
hold on 

plot(T0,Y3)
plot(T0,Y0,'ro') %,'LineWidth',2
xlabel('time','FontSize',12)
ylabel('X','FontSize',12)
set(gca,'FontSize',12)    