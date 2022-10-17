% Stability analysis lecture

clc;clear;close all;


r = 0.2;
K=1;
yls=[-r/2,2*r];

A=0.25;

x = 0:0.01:1;


exponential = x.*r;
logistic = x.*(1-x./K).*r;
allee = x.*(1-x./K).*(x./A - 1).*r;


%% exponential
figure(1); hold on;
plot(x,exponential,'-k','LineWidth',3);
plot(x,x.*0,':k','LineWidth',3);

% add fixed points:
plot(0,0,'.k','MarkerSize',60);
ylim(yls);

%% logistic
figure(2); hold on;
plot(x,logistic,'-k','LineWidth',3);
plot(x,x.*0,':k','LineWidth',3);

% add fixed points:
plot(0,0,'.k','MarkerSize',60);
plot(K,0,'.k','MarkerSize',60);
ylim(yls);

%% logistic w/ allee effect
figure(3); hold on;
plot(x,allee,'-k','LineWidth',3);
plot(x,x.*0,':k','LineWidth',3);

% add fixed points:
plot(0,0,'.k','MarkerSize',60);
plot(A,0,'o','MarkerSize',18,'MarkerFaceColor',[1,1,0.9725],'MarkerEdgeColor','k','LineWidth',3);
plot(K,0,'.k','MarkerSize',60);
ylim(yls);
tspan = [0 10];

%% add some initial conditions
for x0 = 0:0.1:1
    
    x00 = [x0,x0,x0]; %same initial condition for all 3 eqn's
    [t, y]=ode45(@(t,x)all_odes(t,x,r,K,A), tspan, x00);
    
    figure(4); hold on;
    plot(t,y(:,1),'-k','LineWidth',3); hold on;

    figure(5); hold on;
    plot(t,y(:,2),'-k','LineWidth',3); hold on;


    figure(6); hold on;
    plot(t,y(:,3),'-k','LineWidth',3); hold on;


end

function xdot = all_odes(t,x,r,K,A)

    xdot = zeros(3,1);

    % exponential
    X=x(1);
    xdot(1,1) = X*r;
    
    % logistic
    X=x(2);
    xdot(2,1) = X*r*(1-X / K);
    
    % RESISTANT
    X=x(3);
    xdot(3,1) = X*r*(1-X / K)*(X / A - 1);
    
end



