% Combined Plots 
set(groot,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0 0],'DefaultLineLineWidth',2);
set(groot,'DefaultAxesFontSize',20,'defaultAxesFontName', 'Times','defaultTextFontName', 'Times')

%%%% Case 1 
load("Case1.mat");

% Unpack Case
LQR = Case1.LQR;
LF = Case1.LF;
MPC = Case1.MPC;

% Unpack LQR
X_LQR = LQR.X;
U_LQR = LQR.U;
t_LQR = LQR.t;

% Unpack LF
X_LF = LF.X;
U_LF = LF.U;
t_LF = LF.t;

% Unpack MPC
X_MPC = MPC.X;
U_MPC = MPC.U;
t_MPC = MPC.t;

% Relative Trajectory 
figure;
hold on
grid on
view(3);
start = plot3(X_LQR(1,1),X_LQR(2,1),X_LQR(3,1),'mo','LineWidth',4);
p1 = plot3(X_LQR(1,:),X_LQR(2,:),X_LQR(3,:),'r--','LineWidth',2);
p2 = plot3(X_LF(1,:),X_LF(2,:),X_LF(3,:),'b-.','LineWidth',2);
p3 = plot3(X_MPC(1,:),X_MPC(2,:),X_MPC(3,:),'k:','LineWidth',2);
stop = plot3(0,0,0,'mx','LineWidth',4);
title('Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
zlabel('z (in Km)');
legend('Start','LQR','LF','MPC','End')
hold off

% Inplane Relative Motion
figure;
hold on
grid on
start = plot(X_LQR(1,1),X_LQR(2,1),'mo','LineWidth',4);
p1 = plot(X_LQR(1,:),X_LQR(2,:),'r--','LineWidth',2);
p2 = plot(X_LF(1,:),X_LF(2,:),'b-.','LineWidth',2);
p3 = plot3(X_MPC(1,:),X_MPC(2,:),X_MPC(3,:),'k:','LineWidth',2);
stop = plot3(0,0,0,'mx','LineWidth',4);
title('Inplane Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
legend('Start','LQR','LF','MPC','End')
hold off

% Position 
figure;
suptitle('Relative Position');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_LQR,X_LQR(1,:),'r--','LineWidth',2);
p2 = plot(t_LF,X_LF(1,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,X_MPC(1,:),'k:','LineWidth',2);
legend('LQR','LF','MPC')
xlabel('t (in s)');
ylabel('x (in Km)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_LQR,X_LQR(2,:),'r--','LineWidth',2);
p2 = plot(t_LF,X_LF(2,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,X_MPC(2,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('y (in Km)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_LQR,X_LQR(3,:),'r--','LineWidth',2);
p2 = plot(t_LF,X_LF(3,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,X_MPC(3,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('z (in Km)');
hold off

% Velocity
figure;
suptitle('Relative Velocity');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_LQR,X_LQR(4,:),'r--','LineWidth',2);
p2 = plot(t_LF,X_LF(4,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,X_MPC(4,:),'k:','LineWidth',2);
legend('LQR','LF','MPC')
xlabel('t (in s)');
ylabel({"$\dot{x}$ (in Km/s)"},'Interpreter','latex');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_LQR,X_LQR(5,:),'r--','LineWidth',2);
p2 = plot(t_LF,X_LF(5,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,X_MPC(5,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel({"$\dot{y}$ (in Km/s)"},'Interpreter','latex');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_LQR,X_LQR(6,:),'r--','LineWidth',2);
p2 = plot(t_LF,X_LF(6,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,X_MPC(6,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel({"$\dot{z}$ (in Km/s)"},'Interpreter','latex');
hold off

% Control
figure;
suptitle('Control Input');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_LQR,U_LQR(1,:),'r--','LineWidth',2);
p2 = plot(t_LF,U_LF(1,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,U_MPC(1,:),'k:','LineWidth',2);
legend('LQR','LF','MPC')
xlabel('t (in s)');
ylabel('Ux (in KN)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_LQR,U_LQR(2,:),'r--','LineWidth',2);
p2 = plot(t_LF,U_LF(2,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,U_MPC(2,:),'k:','LineWidth',2);xlabel('t (in s)');
xlabel('t (in s)');
ylabel('Uy (in KN)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_LQR,U_LQR(3,:),'r--','LineWidth',2);
p2 = plot(t_LF,U_LF(3,:),'b-.','LineWidth',2);
p3 = plot(t_MPC,U_MPC(3,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('Uz (in KN)');
hold off


%%%% Case 2
load("Case2.mat");

% Unpack Case
LQR = Case2.LQR;
MPC = Case2.MPC;

% Unpack LQR
X_LQR = LQR.X;
U_LQR = LQR.U;
t_LQR = LQR.t;

% Unpack MPC
X_MPC = MPC.X;
U_MPC = MPC.U;
t_MPC = MPC.t;

% Relative Trajectory 
figure;
hold on
grid on
view(3);
start = plot3(X_LQR(1,1),X_LQR(2,1),X_LQR(3,1),'mo','LineWidth',4);
p1 = plot3(X_LQR(1,:),X_LQR(2,:),X_LQR(3,:),'r--','LineWidth',2);
p3 = plot3(X_MPC(1,:),X_MPC(2,:),X_MPC(3,:),'k:','LineWidth',2);
stop = plot3(0,0,0,'mx','LineWidth',4);
title('Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
zlabel('z (in Km)');
legend('Start','LQR','MPC','End')
hold off

% Inplane Relative Motion
figure;
hold on
grid on
start = plot(X_LQR(1,1),X_LQR(2,1),'mo','LineWidth',4);
p1 = plot(X_LQR(1,:),X_LQR(2,:),'r--','LineWidth',2);
p3 = plot3(X_MPC(1,:),X_MPC(2,:),X_MPC(3,:),'k:','LineWidth',2);
stop = plot3(0,0,0,'mx','LineWidth',4);
title('Inplane Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
legend('Start','LQR','MPC','End')
hold off

% Position 
figure;
suptitle('Relative Position');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_LQR,X_LQR(1,:),'r--','LineWidth',2);
p3 = plot(t_MPC,X_MPC(1,:),'k:','LineWidth',2);
legend('LQR','MPC')
xlabel('t (in s)');
ylabel('x (in Km)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_LQR,X_LQR(2,:),'r--','LineWidth',2);
p3 = plot(t_MPC,X_MPC(2,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('y (in Km)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_LQR,X_LQR(3,:),'r--','LineWidth',2);
p3 = plot(t_MPC,X_MPC(3,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('z (in Km)');
hold off

% Velocity
figure;
suptitle('Relative Velocity');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_LQR,X_LQR(4,:),'r--','LineWidth',2);
p3 = plot(t_MPC,X_MPC(4,:),'k:','LineWidth',2);
legend('LQR','MPC')
xlabel('t (in s)');
ylabel({"$\dot{x}$ (in Km/s)"},'Interpreter','latex');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_LQR,X_LQR(5,:),'r--','LineWidth',2);
p3 = plot(t_MPC,X_MPC(5,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel({"$\dot{y}$ (in Km/s)"},'Interpreter','latex');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_LQR,X_LQR(6,:),'r--','LineWidth',2);
p3 = plot(t_MPC,X_MPC(6,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel({"$\dot{z}$ (in Km/s)"},'Interpreter','latex');
hold off

% Control
figure;
suptitle('Control Input');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_LQR,U_LQR(1,:),'r--','LineWidth',2);
p3 = plot(t_MPC,U_MPC(1,:),'k:','LineWidth',2);
legend('LQR','MPC')
xlabel('t (in s)');
ylabel('Ux (in KN)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_LQR,U_LQR(2,:),'r--','LineWidth',2);
p3 = plot(t_MPC,U_MPC(2,:),'k:','LineWidth',2);xlabel('t (in s)');
xlabel('t (in s)');
ylabel('Uy (in KN)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_LQR,U_LQR(3,:),'r--','LineWidth',2);
p3 = plot(t_MPC,U_MPC(3,:),'k:','LineWidth',2);
xlabel('t (in s)');
ylabel('Uz (in KN)');
hold off