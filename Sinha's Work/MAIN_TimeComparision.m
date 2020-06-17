% Time Comparision MAIN.m
clc; clear all;


%%%%% Using Sinha's Method
m = 13;
ti = 0;
tf = 2*pi;
casenum = "demo3";

switch casenum
    case "demo1"
        n = 1;
        aperiodic_mat = 0;
        periodic_mat = {@(t) cos(t)};
        
    case "demo2"
        % Demo System 2
        n = 2;
        aperiodic_mat = [0, 1;...
                                 0, 0]; 
        periodic_mat = {@(t)cos(3*t),@(t) sin(t);...
                                @(t) 0,@(t) 2*sin(t)+cos(t)};
    case "demo3"
        % Demo System 3
        n = 2;
        aperiodic_mat = [0, 0;...
                                 0, 0]; 
        periodic_mat = {@(t)cos(t),@(t) sin(t);...
                                @(t) -sin(t),@(t) cos(t)};
                            
    case "demo4"
        global alp beta w;
        alp = sqrt(2.5);
        beta = 2;
        w = 1;
        n = 2;
        aperiodic_mat = [0, 1;...
                         -alp^2, 0]; 
        periodic_mat = {@(t)0,@(t) 0;...
                        @(t) -beta*cos(w*t),@(t) 0};
end

G_firstkind = double(OpMat_Int_firstkind(m));
[PhiT_cheb_firstkind, time_cheb_firstkind] = getPhiT_firstkind(m,n,aperiodic_mat,periodic_mat,ti,tf,G_firstkind);
G_secondkind = double(OpMat_Int_secondkind(m));
[PhiT_cheb_secondkind, time_cheb_secondkind] = getPhiT_secondkind(m,n,aperiodic_mat,periodic_mat,ti,tf,G_secondkind);

%%%%% Using ode45

% Using n initial conditions to find PhiT
PhiT_ode45 = zeros(n,n);
tic;
X0 = eye(n,n);
for i = 1:n
    [~,xt] = ode45(@(t,X) system_for_ode(t,X,casenum), [ti, tf],X0(:,i));
    PhiT_ode45(:,i) = xt(end,:)';
end
time_ode45 = toc;

function dX = system_for_ode(t,X,casenum)

    switch casenum
        case "demo1"
             % demo 1 system 
            dX = cos(t)*X;
        case "demo2"
            % demo 2 system
            dX(1,1) = cos(3*t)*X(1) + (1+sin(t))*X(2);
            dX(2,1) = 0*X(1) + (2*sin(t)+cos(t))*X(2);
        case "demo3"
            % demo 3 system
            dX(1,1) = cos(t)*X(1) + sin(t)*X(2);
            dX(2,1) = -sin(t)*X(1) + cos(t)*X(2);
        case "demo4"
            global alp beta w;
            dX(1,1) = 0*X(1) + 1*X(2);
            dX(2,1) = -(alp^2+beta*cos(w*t))*X(1) + 0*X(2);
    end
end
