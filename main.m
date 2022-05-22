%% Clear
clear;close all;clc;

%% Initializing Inputs
inputs = struct(...
    'chord', 1,...
    'Vinf', 100,...
    'max_thickness', 0.07,...
    'max_camber', 0.05,...
    'alpha_deg', 8,...
    'i_max', 61,...
    'j_max', 121,...
    'n', 30000,...
    'R', 5);

%% Initializing Airfoil Class
airfoil = Airfoil(inputs);

% Calculate Joukowski Airfoil
[~, outerCircle, joukowski] = airfoil.joukowskiAirfoil();

% Plot Joukowski Airfoil
figure;
plot(joukowski.x(1:end/2), joukowski.y(1:end/2),'g-','LineWidth',1); hold on
plot(joukowski.x(end/2: end), joukowski.y(end/2: end), 'b-','LineWidth',1);
xlabel('$X_{B}$ $for$ $Airfoil$','interpreter','latex','FontSize',14);
ylabel('$Y_{B}$ $for$ $Airfoil$','interpreter','latex','FontSize',14);
title('$Airfoil$ $shape$ $approximate$','interpreter','latex','FontSize',14);
axis equal
grid on

% Calculate & Plot O-Grid
figure;
hold on;
plot(outerCircle.x, outerCircle.y, '-r', 'LineWidth', 3);
plot(joukowski.x, joukowski.y, '-b', 'LineWidth', 2);
axis equal;
[xGrid, yGrid, ~] = airfoil.generatePhysicalGrid(outerCircle, joukowski);

% Calculate & Plot Computational Grid
figure;
[eta1Grid, eta2Grid] = airfoil.generateComputationalGrid([0, 1], [0, 1]);

% Calculate Metric derivatives x/eta1, y/eta1, x/eta2, y/eta2 
% And Calculate Values of C11, C22, C12, J
airfoil.transformationMetrics(xGrid, yGrid, eta1Grid, eta2Grid);

% Calculate & return Boundary Conditions
psi_ = airfoil.calculateDirichletBoundary();
n = 1;
errorLog10 = ones(1, airfoil.inputs.n);
error = ones(1, airfoil.inputs.n);
while((n <= airfoil.inputs.n) && (min(error) > 1e-8))
    
    psi_new = airfoil.iterate(psi_);
    psi_new = psi_ + 1.1 .* (psi_new - psi_);
    
    error(n) = max(max(abs(psi_new-psi_)));
    if error(n) > 0; errorLog10(n) = log10(error(n)); end
    
    psi_ = psi_new;
    psi_(:,1) = psi_(1,2);
    n = n + 1;
end

%% Plot Error History
figure;
plot(errorLog10,'k')
grid on
xlabel('Iteration number', 'fontsize',14)
ylabel('Log_1_0 (Error)', 'fontsize',14)
title('Convergence history using Point-SOR for alpha = 8^o (O-Grid)','fontsize',12)

%% Calculate velocity and Cp

dpsi_deta1 = airfoil.zerosImaxJmax();
dpsi_deta2 = airfoil.zerosImaxJmax();
for i=1:airfoil.inputs.i_max
     if i ==1
         dpsi_deta1(i,:) = (psi_(i+1,:)-psi_(i,:))./(eta1Grid(i+1,:)-eta1Grid(i,:));
     elseif i==airfoil.inputs.i_max
         dpsi_deta1(i,:) = (psi_(i,:)-psi_(i-1,:))./(eta1Grid(i,:)-eta1Grid(i-1,:));
     else 
         dpsi_deta1(i,:) = (psi_(i+1,:)-psi_(i-1,:))./(eta1Grid(i+1,:)-eta1Grid(i,:));
     end 
 end 
for j=1:airfoil.inputs.j_max
     if j ==1
         dpsi_deta2(:,j) = (psi_(:,j+1)-psi_(:,j))./(eta2Grid(:,j+1)-eta2Grid(:,j));
     elseif j==airfoil.inputs.j_max
          dpsi_deta2(:,j) = (psi_(:,j)-psi_(:,j-1))./(eta2Grid(:,j)-eta2Grid(:,j-1));
     else 
         dpsi_deta2(:,j) = (psi_(:,j+1)-psi_(:,j-1))./(eta2Grid(:,j+1)-eta2Grid(:,j-1));
     end 
end 
u = dpsi_deta1.*airfoil.deta1_dy + dpsi_deta2.*airfoil.deta2_dy;
v = -(dpsi_deta1.*airfoil.deta1_dx + dpsi_deta2.* airfoil.deta2_dx);
 
V = sqrt(u.^2+v.^2);
nonDimV = V/airfoil.inputs.Vinf;
C_p = 1-(nonDimV).^2;

%% Plotting Results
figure
quiver(xGrid,yGrid,u, v, 'k-', 'LineWidth', 0.5);hold on
plot(joukowski.x, joukowski.y,'r-','LineWidth',2) ; hold on
axis equal
xlim([-(airfoil.inputs.chord+0.5)/2 (airfoil.inputs.chord+0.5)/2])
xlabel('X-axis', 'fontsize',14)
ylabel('Y-axis', 'fontsize',14)
title('Velocity vector for angle of attack =8^o (Numerical Solution - O Grid)','fontsize',12)

figure;
contour(xGrid, yGrid, C_p, 5000);
title('Contours of Pressure coefficient (angle of attack =8^o) (Numerical Solution - O Grid)','fontsize',12)

figure
plot(xGrid(1:end/2,1),nonDimV(1:end/2,1), 'LineWidth', 1.5);
hold on;
plot(xGrid(end/2:end,1),nonDimV(end/2:end,1), 'LineWidth', 1.5);
xlabel('$X_{B}$','interpreter','latex','FontSize',14);
ylabel('$\frac{V}{V_{inf}}$','interpreter','latex','FontSize',14);
title('$Velocity$ $over$ $airfoil$ Using Numerical Solution (O-Grid)','interpreter','latex','FontSize',14);

figure;
plot(xGrid(1:end/2,1),C_p(1:end/2,1), 'LineWidth', 1.5);
hold on;
plot(xGrid(end/2:end,1),C_p(end/2:end,1), 'LineWidth', 1.5);
ylim([-4 1.3]);
xlabel('$X_{B}$','interpreter','latex','FontSize',14);
ylabel('$C_{p}$','interpreter','latex','FontSize',14);
title('$Pressure$ $over$ $airfoil$ Using Numerical Solution (O-Grid)','interpreter','latex','FontSize',14);
