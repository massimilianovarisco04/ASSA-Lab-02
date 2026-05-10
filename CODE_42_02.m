clear;
close all;
clc;

%% Task 1.1 - Nonlineare Simulation
% Parametri per il sistema e condizioni iniziali
scenari = ['A', 'B', 'C', 'D', 'E'];
proximityP = proximity_parameters();
n = proximityP.n;
T = proximityP.T;
%mettiamo tutte le condizioni iniziali relative alla stessa coordinata di
%stato un un vettore
x_0 = [0.1, 0.1, 0, 0.1, 15];
y_0 = [0, 0, 0.05, 0, 0];
z_0 = [0, 0, 0.1, 0.05, 5];
x_dot_0 = [0, 0, 0, 0, 0];
y_dot_0 = [-2*n*x_0(1), -1.5*n*x_0(2), 0, -1.5*n*x_0(4), -2*n*x_0(5)];
z_dot_0 = [n*x_0(1), 0, 0, 0, 0];
%costruiamo la matrice delle condizioni iniziali
x0 = [x_0; x_dot_0; y_0; y_dot_0; x_0; x_dot_0];

t_0 = 0;
t_f = 2*T;

% Risoluzione sistema non lineare
for i = 1 : size(x0,2)
    ODE_obj = ode;  
    ODE_obj.ODEFcn = @(t,x) proximityP_f(t,x, proximityP);
    ODE_obj.InitialValue = x0(:,i).*10^3;
    ODE_obj.Solver = 'ode45';
    ODEResults_obj = solve(ODE_obj, t_0, t_f);
    tt = ODEResults_obj.Time'; 
    xx = (ODEResults_obj.Solution'); 

    % Solo per verifica
    if i == 1
        figure('Name',['Trajectory ', scenari(i)]);
        plot3(xx(:,1), xx(:,3), xx(:,5), 'LineWidth',1.5);
        grid on;
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        title(['Trajectory', scenari(i)]);
    end

    figure('Name',['Scenario ', scenari(i)]);
    subplot(2,1,1)
    plot(tt, xx(:,1), 'LineWidth',2);
    hold on;
    plot(tt, xx(:,3), 'LineWidth',2);
    plot(tt, xx(:,5), 'LineWidth',2);
    grid on;
    xlabel('Time[s]');
    ylabel('Position[m]');
    legend('x(t)', 'y(t)', 'z(t)');
    title(['Coordinates case ', scenari(i)]);

    subplot(2,1,2)
    plot(tt, xx(:,2), 'LineWidth',2);
    hold on;
    plot(tt, xx(:,4), 'LineWidth',2);
    plot(tt, xx(:,6), 'LineWidth',2);
    grid on;
    xlabel('Time[s]');
    ylabel('Speed[m/s]');
    legend('x_dot', 'y_dot', 'z_dot');
    title(['Coordinates of speed case ', scenari(i)]);
end

%% 

%% ----------------------- Definizione Funzioni ---------------------------
% Struttura contenente dati del problema
function proximityP = proximity_parameters()
    mu = 398600.44; %km^3/s^2
    Re = 6378.14; %km
    h0 = 500;  %km
    R0 = Re + h0;
    n = sqrt(mu/R0^3);
    T = 2*pi/n;
    
    proximityP.mu = mu;
    proximityP.Re = Re;
    proximityP.h0 = h0;
    proximityP.R0 = R0;
    proximityP.n = n;
    proximityP.T = T;
end

% Funzione per scrivere il sistema in spazio di stato
function xdot = proximityP_f (t,x, proximityP)
    mu = proximityP.mu;
    R0 = proximityP.R0;
    n = proximityP.n;

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    x6 = x(6);

    rc = sqrt((R0+x1)^2 + x3^2 + x6^2);

    xdot_1 = x2;
    xdot_2 = 2*n*x4 + n^2*(R0+x1) - (mu*(R0+x1))/rc^3;
    xdot_3 = x4;
    xdot_4 = -2*n*x2 + n^2*x3-(mu*x3)/(rc^3);
    xdot_5 = x6;
    xdot_6 = -(mu*x5)/rc^3;

    xdot = [xdot_1, xdot_2, xdot_3, xdot_4, xdot_5, xdot_6]';
end