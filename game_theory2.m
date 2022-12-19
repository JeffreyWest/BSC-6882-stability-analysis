clc;clear;close all;

A1=[0,0;1,1];
A2=[0,0;1,-1];
A3=[0,0;-1,-1];
A4=[0,0;-1,1];
% pairwise_fixedpoint(A); hold on;


A=A4;
tspan=[0,10];
for x00 = 0.05:0.1:0.95
    [t, x]=ode45(@(t,n)replicator(t,n,A), tspan, [x00,1-x00]);
    plot(t,x(:,1),'-r','LineWidth',3); hold on;
    plot(t,x(:,2),'--b','LineWidth',3);
end
clean();
printPNG(figure(1),'t4.png');


function xdot = replicator(t,x,A)
    x1=x(1); x2=x(2);
    f1 = A(1,1)*x1 + A(1,2)*x2;
    f2 = A(2,1)*x1 + A(2,2)*x2;
    phi = x1*f1 + x2*f2;

    xdot(1,1) = x1*(f1- phi); 
    xdot(2,1) = x2*(f2- phi);    
end