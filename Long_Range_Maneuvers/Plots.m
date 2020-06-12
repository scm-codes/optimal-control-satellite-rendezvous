%%% Plots
%%
clear all
set(groot,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0 0],'DefaultLineLineWidth',2);
set(groot,'DefaultAxesFontSize',20,'defaultAxesFontName', 'Times','defaultTextFontName', 'Times')
%%


figure;clf;

tgrid = linspace(0,1590.21866726184,1000);
% 
% trap = load("Trapezoid.mat");
% state = trap.soln.interp.state(tgrid);
% control = trap.soln.interp.control(tgrid);
% suptitle("Result: Trapezoidal Collocation");
% 
% cheb = load("chebyshev.mat");
% state = cheb.soln.interp.state(tgrid);
% control = cheb.soln.interp.control(tgrid);
% suptitle("Result: Chebyshev Collocation");

% hesi = load("HermiteSimpson.mat");
% state = hesi.soln.interp.state(tgrid);
% control = hesi.soln.interp.control(tgrid);
% suptitle("Result: Hermite Simpson Collocation");

rk4 = load("rungekutta.mat");
state = rk4.soln.interp.state(tgrid);
control = rk4.soln.interp.control(tgrid);
suptitle("Result: RK-4 Collocation");
%%%% x

subplot(3,2,1)
hold on
plot(tgrid,state(1,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$x$ [in Km]"},'Interpreter','latex');
hold off

subplot(3,2,3)
hold on
plot(tgrid,state(4,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$\dot{x}$ [in Km/s]"},'Interpreter','latex');
hold off



subplot(3,2,5)
hold on
plot(tgrid,control(1,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$U_{x}$ [in KN]"},'Interpreter','latex');
hold off

%%%% y
subplot(3,2,2)
hold on
plot(tgrid,state(2,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$y$ [in Km]"},'Interpreter','latex');
hold off

subplot(3,2,4)
hold on 
plot(tgrid,state(5,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$\dot{y}$ [in Km/s]"},'Interpreter','latex');
hold off

subplot(3,2,6)
hold on 
plot(tgrid,control(2,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$U_{y}$ [in KN]"},'Interpreter','latex');
hold off


%%%% z
figure;clf;
subplot(2,2,1)
hold on
plot(tgrid,state(3,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"z [in Km]"},'Interpreter','latex');
hold off

subplot(2,2,2)
hold on
plot(tgrid,state(6,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$\dot{z}$ [in Km/s]"},'Interpreter','latex');
hold off

subplot(2,2,[3,4])
hold on 
plot(tgrid,control(3,:),'LineWidth',2);
xlabel({"t [in sec]"},'Interpreter','latex');
ylabel({"$U_{z}$ [in KN]"},'Interpreter','latex');
hold off


% Flatness Stats Plots
load('Flatness_Stats.mat')
J1_cost = Flatness_Stats.minVel.objVal(5:20);
J2_cost = Flatness_Stats.minAcc.objVal(5:20);
J1_nlpTime = Flatness_Stats.minVel.nlpTime(5:20);
J2_nlpTime = Flatness_Stats.minAcc.nlpTime(5:20);
harmonics = 5:20;

figure;
subplot(2,1,1)
hold on
grid on
p1 = plot(harmonics,J1_cost,'ro','LineWidth',2);
p2 = plot(harmonics,J2_cost,'kx','LineWidth',2);
p3 = plot(harmonics(9),J2_cost(9),'bs','LineWidth',2,'MarkerSize',20);
xlabel('No of Harmonics');
ylabel({'Area[$U^2$ vs t](in KN-s)'},'Interpreter','latex');
legend("J_1","J_2","Lowest Value")
hold off

subplot(2,1,2)
hold on
grid on
p1 = plot(harmonics,J1_nlpTime,'ro','LineWidth',2);
p2 = plot(harmonics,J2_nlpTime,'kx','LineWidth',2);
xlabel('No of Harmonics');
ylabel('NLP time (in s)');
hold off
