clc;clear;close all;


a=3; b=0;
c=5; d=1;

A1=[a,b;c,d];
tspan = [0 10];
x0 = [0.5,0.5];

% A1:
[t, x]=ode45(@(t,n)replicator(t,n,A1), tspan, x0);
plot(t,x(:,1),'-k','LineWidth',5); hold on;
plot(t,x(:,2),'-r','LineWidth',5);


% A2:
A2=[0,0;-(a-c),-(b-d)]
[t, x]=ode45(@(t,n)replicator(t,n,A2), tspan, x0);
plot(t,x(:,1),'--m','LineWidth',3); hold on;
plot(t,x(:,2),'--g','LineWidth',3);

clean();
ylim([0 1])

function xdot = replicator(t,x,A)
    x1=x(1); x2=x(2);
    f1 = A(1,1)*x1 + A(1,2)*x2;
    f2 = A(2,1)*x1 + A(2,2)*x2;
    phi = x1*f1 + x2*f2;

    xdot(1,1) = x1*(f1- phi); 
    xdot(2,1) = x2*(f2- phi);    
end



% function xdot = replicator(t,x,A)
%     x1=x(1); x2=x(2);
%     f1 = A(1,1)*x1 + A(1,2)*x2;
%     f2 = A(2,1)*x1 + A(2,2)*x2;
%     phi = x1*f1 + x2*f2;
% 
%     xdot(1,1) = x1*x2*(f1- f2); 
%     xdot(2,1) = x2*x1*(f2- f1);    
% end