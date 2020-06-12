% MAIN.m
%

clc; clear;
CurrentPath = pwd;
addpath('CurrentPath\OptimTraj-master')
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

% info
numpoints = 5000;
start = 0;
stop = 0.251*params.target.T;
params.target.tspan = linspace(start,stop,numpoints);

[params.target.theta,...
    params.target.dtheta,...
    params.target.ddtheta] = target_info(params.target);

[params.chaser.h,...
    params.chaser.T,...
    params.chaser.n,...
    params.chaser.tp] = InferParams(params.chaser);

X0 = rho0_info(params);


%-------------------------------------------------------------------------% 
%                             Setup function                              %
%-------------------------------------------------------------------------%

problem.func.dynamics = @(t,X,U) (RendezvousDynamics(t,X,U,params));
problem.func.pathObj = @(t,X,U) (pathObj_forcesquared(U));


%-------------------------------------------------------------------------% 
%                             Setup bounds                                %
%-------------------------------------------------------------------------%

maxForce = 1;
duration = 0.25*params.target.T;

problem.bounds.initialTime.low = start;
problem.bounds.initialTime.upp = start;
problem.bounds.finalTime.low = duration;
problem.bounds.finalTime.upp = duration;

problem.bounds.initialState.low = X0;
problem.bounds.initialState.upp = X0;
problem.bounds.finalState.low = zeros(6,1);
problem.bounds.finalState.upp = [0.01;0.001;0.001;0.1;0;0];

problem.bounds.state.low = [-2*X0(1);-0.2*X0(1);-0.2*X0(1);-inf;-inf;-inf];
problem.bounds.state.upp = [2*X0(1);0.2*X0(1);0.2*X0(1);inf;inf;inf];

problem.bounds.control.low = -maxForce*ones(3,1);
problem.bounds.control.upp = maxForce*ones(3,1);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,duration];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = zeros(3,2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'TolFun',1e-5,...
    'MaxFunEvals',1e6);

% problem.options(1).method = 'trapezoid';
% problem.options(1).trapezoid.nGrid = 40;
% 
% problem.options(1).method = 'hermiteSimpson';
% problem.options(1).hermiteSimpson.nSegment = 40;
% 
% problem.options(1).method = 'rungeKutta';
% problem.options(1).rungeKutta.nSegment = 40;

problem.options.method = 'chebyshev';
problem.options.chebyshev.nColPts = 40;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);

