function dy = syndeg_simple(t,y,P)
%Parameters
a = P(1); %mRNA production
b = P(2); %mRNA degradation


%variables
x=y(1);

%ODE
dy(1) = a - b*x;