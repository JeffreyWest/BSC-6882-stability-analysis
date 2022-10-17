%% Stability analysis lecture
% step 1: Find all fixed points ($\dot{x} = 0$)
% step 2: Find the generalized Jacobian matrix of $F$, $\matr{J}(\matr{x})$
% step 3: Find Jacobian matrix eigenvalues evaluated at stationary points, $\matr{J}(\matr{x}^*)$
% step 4: Classify stability type for each zero according to the sign of eigenvalues

clc;clear;close all;

%% treated, with mutual competition:
c = 2; r1 = 0.1; r2 = 0.1;
K1 = 1; K2 = 1;
d1 = 0.05; d2 = 0.05;
m12 = 1.2; m21 = 1.2;


%% untreated, with mutual competition
c = 0; r1 = 0.1; r2 = 0.1;
K1 = 1; K2 = 1;
d1 = 0.05; d2 = 0.05;
m12 = 2; m21 = 3;



fixed_points = getFixedPoints(c,m12,m21,d1,d2,r1,r2,K1,K2);

figure(1); hold on;
Quiver(c,m12,m21,d1,d2,r1,r2,K1,K2)

% iterate over all fixed points
for fi = 1:1:size(fixed_points,1)
    fixed_point = fixed_points(fi,:)
    Jstar = jacobian(fixed_point,c,m12,m21,d1,d2,r1,r2,K1,K2);

    [lambda1,lambda2] = real_part_eigenvalues(Jstar);
    
    [lambda1,lambda2]
    check_fixedpoint_type_and_plot(lambda1,lambda2,fixed_point);
    
end



%% step 5: plot
clean();

function [fixedpoints] = getFixedPoints(c,m12,m21,d1,d2,r1,r2,K1,K2)

    % we know there are 4 fixed points:
    fixedpoints = zeros(4,2);

    R1 = r1*(1-c); % for conveniece of not having to re-type this
    
    %% origin
    x1 = 0;
    x2 = 0;
    fixedpoints(1,:) = [x1,x2];

    %% on y-axis: (0,(r1-d1)/r1*k))
    x1=(R1-d1)/R1*K1;
    x2=0;
    fixedpoints(2,:) = [x1,x2];
    
    %% on x-axis: (0,(r2-d2)/r2*K))
    x1=0;
    x2=(r2-d2)/r2*K2;
    fixedpoints(3,:) = [x1,x2];
    
    %% internal fixed point
    del1 = (R1-d1)*K1/R1;
    del2 = (r2-d2)*K2/r2;
    x1 = (del1-m21*del2)/(1-m12*m21);
    x2 = (del2-m12*del1)/(1-m12*m21);
    fixedpoints(4,:) = [x1,x2];
    
end

function [j] = jacobian(fixed_point,c,m12,m21,d1,d2,r1,r2,K1,K2)

x1 = fixed_point(1);
x2 = fixed_point(2);

j = [r1*(1-c)-d1 - r1*(1-c)*m21/K1*x2 - 2*r1*(1-c)*x1/K1, -r1*(1-c)*m21*x1/K1;
     -r2*m12*x2/K2, r2-d2-r2*m12*x1/K2 - 2*r2*x2/K2];

end

function [e1,e2] = real_part_eigenvalues(J)
    [eigenvectors,eigenvalues] = eig(J);

    e1=real(eigenvalues(1,1));
    e2=real(eigenvalues(2,2));

end

function [] = check_fixedpoint_type_and_plot(lambda1,lambda2,fixed_point)
    if (lambda1 == 0) || (lambda2 == 0)
        disp('This fixed point has a zero eigenvalue');
        plot(fixed_point(1),fixed_point(2),'o','MarkerSize',18,'MarkerFaceColor',[.5,.5,.5],'MarkerEdgeColor','k','LineWidth',3);

    elseif max([lambda1,lambda2]>0)
        disp('This fixed point is unstable.');
        plot(fixed_point(1),fixed_point(2),'o','MarkerSize',18,'MarkerFaceColor',[1,1,0.9725],'MarkerEdgeColor','k','LineWidth',3);
    else
        disp('This fixed point is a stable attractor');
        plot(fixed_point(1),fixed_point(2),'.k','MarkerSize',60);
    end
end

function [] = Quiver(c,m12,m21,d1,d2,r1,r2,K1,K2)

    color = [0,0,0];

    fp = getFixedPoints(c,m12,m21,d1,d2,r1,r2,K1,K2);
    Keffective1 = max(fp(:,1));
    Keffective2 = max(fp(:,2));
    Kmax = max([Keffective1,Keffective2]);
    Kmax = max(Kmax,1);

    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    
    %% set up meshgrid
    step = Kmax/20;
    step0 = step/2;
    x_grid = step0:step:(Kmax-step0);

    step = Kmax/20;
    step0 = step/2;
    y_grid = step0:step:(Kmax-step0);

    [X1,X2] = meshgrid(x_grid,y_grid); % Generate domain.
    
    %% SENSITIVE    
    x1dot = r1.*X1.*( 1 - (X1 + m21.*X2) ./ K1 ).*(1 - c) - d1.*X1;
    x2dot = r2.*X2.*( 1 - (X2 + m12.*X1) ./ K2 ) - d2.*X2;

    
    
    %% arrow length:

    len = 0.0375;
    HL = 10;
    HW = 8;

    for i = 1:length(X1)
        for j = 1:length(X2)                
                
            
            %% calculate direction & magnitude
            V = x2dot(i,j);
            U = x1dot(i,j); 
            mag = sqrt(V^2 + U^2);
            
            if (mag > 0)
            %% set direction & magnitude
                ah = annotation('arrow','headStyle','cback1','HeadLength',HL,'HeadWidth',HW);

                set(ah,'parent',gca);
                set(ah,'position',[X1(i,j), X2(i,j), len*U/mag, len*V/mag]);
                set(ah,'Color',color);
            end
        end
    end
    
    
    xlim([0 Kmax]);
    ylim([0 Kmax]);
end

% make the plots pretty
function clean()
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','Bold', 'LineWidth', 2);
    
    p = get(gcf,'Position');
    set(gcf,'Position',[p(1),p(2),600,600]);
    box on; set(gcf,'color','w');
    xlabel('x_1 (Sensitive)')
    ylabel('x_2 (Resistant)')
end





