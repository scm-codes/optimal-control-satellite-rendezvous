% Combined Plots
%%
set(groot,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0 0],'DefaultLineLineWidth',2);
set(groot,'DefaultAxesFontSize',20,'defaultAxesFontName', 'Times','defaultTextFontName', 'Times')
%%
cheb = load("chebyshev.mat");
flat = load("Flatness.mat");
tgrid = linspace(0,1590.21866726184,100);

% Unpack Case
Chebyshev = cheb.soln.interp;
Flatness = flat.soln.interp;

% Unpack LQR
X_Cheb = Chebyshev.state(tgrid);
U_Cheb = Chebyshev.control(tgrid);
t_Cheb = tgrid;

% Unpack MPC
X_Flat = Flatness.state(tgrid)';
U_Flat = Flatness.control(tgrid)';
t_Flat = tgrid;

% Relative Trajectory 
figure;
hold on
grid on
view(3);
start = plot3(X_Cheb(1,1),X_Cheb(2,1),X_Cheb(3,1),'mo','LineWidth',4);
p1 = plot3(X_Cheb(1,:),X_Cheb(2,:),X_Cheb(3,:),'r--','LineWidth',2);
p3 = plot3(X_Flat(1,:),X_Flat(2,:),X_Flat(3,:),'k:','LineWidth',2);
stop = plot3(0,0,0,'mx','LineWidth',4);
title('Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
zlabel('z (in Km)');
legend('Start','Chebyshev','Flatness','End')
hold off

% Inplane Relative Motion
figure;
hold on
grid on
start = plot(X_Cheb(1,1),X_Cheb(2,1),'mo','LineWidth',4);
p1 = plot(X_Cheb(1,:),X_Cheb(2,:),'r--','LineWidth',2);
p3 = plot3(X_Flat(1,:),X_Flat(2,:),X_Flat(3,:),'k:','LineWidth',2);
stop = plot3(0,0,0,'mx','LineWidth',4);
title('Inplane Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
legend('Start','Chebyshev','Flatness','End')
hold off

% Position 
figure;
suptitle('Relative Position');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_Cheb,X_Cheb(1,:),'r--','LineWidth',2);
p3 = plot(t_Flat,X_Flat(1,:),'k:','LineWidth',2);
legend('Chebyshev','Flatness')
xlabel('t (in s)');
ylabel('x (in Km)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_Cheb,X_Cheb(2,:),'r--','LineWidth',2);
p3 = plot(t_Flat,X_Flat(2,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('y (in Km)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_Cheb,X_Cheb(3,:),'r--','LineWidth',2);
p3 = plot(t_Flat,X_Flat(3,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('z (in Km)');
hold off

% Velocity
figure;
suptitle('Relative Velocity');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_Cheb,X_Cheb(4,:),'r--','LineWidth',2);
p3 = plot(t_Flat,X_Flat(4,:),'k:','LineWidth',2);
legend('Chebyshev','Flatness')
xlabel('t (in s)');
ylabel({"$\dot{x}$ (in Km/s)"},'Interpreter','latex');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_Cheb,X_Cheb(5,:),'r--','LineWidth',2);
p3 = plot(t_Flat,X_Flat(5,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel({"$\dot{y}$ (in Km/s)"},'Interpreter','latex');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_Cheb,X_Cheb(6,:),'r--','LineWidth',2);
p3 = plot(t_Flat,X_Flat(6,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel({"$\dot{z}$ (in Km/s)"},'Interpreter','latex');
hold off

% Control
figure;
suptitle('Control Input');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_Cheb,U_Cheb(1,:),'r--','LineWidth',2);
p3 = plot(t_Flat,U_Flat(1,:),'k:','LineWidth',2);
legend('Chebyshev','Flatness')
xlabel('t (in s)');
ylabel('Ux (in KN)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_Cheb,U_Cheb(2,:),'r--','LineWidth',2);
p3 = plot(t_Flat,U_Flat(2,:),'k:','LineWidth',2);xlabel('t (in s)');
xlabel('t (in s)');
ylabel('Uy (in KN)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_Cheb,U_Cheb(3,:),'r--','LineWidth',2);
p3 = plot(t_Flat,U_Flat(3,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('Uz (in KN)');
hold off