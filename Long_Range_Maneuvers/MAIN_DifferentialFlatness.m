% MAIN.m
%
% Uses differential flatness to solve for the non-linear rendezvous   

%
clc;clear;          
global mu;
mu = 3.986004415000000e+05;

%-------------------------------------------------------------------------%
%                      Initial orbital parameters                         %
%-------------------------------------------------------------------------%

% Target Satellite
params.target.mass = 100;
params.target.a = 7420;
params.target.e = 0.12;
params.target.i = 0;
params.target.raan = 0;
params.target.ap = 0;
params.target.theta0 = 0;

% Chaser Satellite
params.chaser.mass = 100;
params.chaser.a = 7440;
params.chaser.e = 0.1;
params.chaser.i = 0;
params.chaser.raan = 0;
params.chaser.ap = 0;
params.chaser.theta0 = 0;


%-------------------------------------------------------------------------%
%                                Inferred                                 %
%-------------------------------------------------------------------------%

[params.target.h,...
    params.target.T,...
    params.target.n,...
    params.target.tp] = InferParams(params.target);

[params.chaser.h,...
    params.chaser.T,...
    params.chaser.n,...
    params.chaser.tp] = InferParams(params.chaser);

X0 = rho0_info(params);

%-------------------------------------------------------------------------%
%                          Solver Inputs                                 %
%-------------------------------------------------------------------------%

%%%% Problem

% Fourier Method
% .fourier.Objective
%         = 1  (min Velocity)
%         = 2  (min Acceleration)
problem.fourier.harmonics = 13;
problem.fourier.frequency = params.target.n;
problem.fourier.Objective = 2;

% Time
duration = 0.25*params.target.T;
problem.time.limits = [0, duration];
problem.time.segments = 40; 
   
% Optmizer
problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'TolFun',1e-5,...
    'MaxFunEvals',1e6,...
    'MaxIter',1e5);


% State
problem.bounds.initialstate = X0;
problem.bounds.finalstate.low = zeros(6,1);
problem.bounds.finalstate.upp = [0.01;0.001;0.001;0.1;0;0];
problem.bounds.state.low = [-2*X0(1);-0.2*X0(1);-0.2*X0(1);-1000;-1000;-1000];
problem.bounds.state.upp = [2*X0(1);0.2*X0(1);0.2*X0(1);1000;1000;1000];

% Control
maxForce = 1;
problem.bounds.control.low = -maxForce;
problem.bounds.control.upp = maxForce;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% nlpTime = zeros(20,2);
% objVal = zeros(20,2);
% for j = 1:2
%     problem.fourier.Objective = j;
%     for i = 5:1:20
%         problem.fourier.harmonics = i;
%         soln = Fourier(problem,params);
%         nlpTime(i,j) = soln.info.nlpTime;
%         objVal(i,j) = soln.info.objVal;
%     end
% end

soln = Fourier(problem,params);
soln.interp.state= @(t) (Fourier_StateInterpolant(problem,t,soln.grid.fouriercoeffs));
soln.interp.control = @(t)(Fourier_ControlInterpolant(problem,params,t,soln.grid.fouriercoeffs));