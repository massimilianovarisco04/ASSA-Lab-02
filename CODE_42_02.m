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
%stato un un vettore. Sono tutte in [Km]!!
x_0 = [0.1, 0.1, 0, 0.1, 15];
y_0 = [0, 0, 0.05, 0, 0];
z_0 = [0, 0, 0.1, 0.05, 5];
x_dot_0 = [0, 0, 0, 0, 0];
y_dot_0 = [-2*n*x_0(1), -1.5*n*x_0(2), 0, -1.5*n*x_0(4), -2*n*x_0(5)];
z_dot_0 = [n*x_0(1), 0, 0, 0, 0];
%costruiamo la matrice delle condizioni iniziali, è 6 righe e 5 colonne
x0 = [x_0; x_dot_0; y_0; y_dot_0; z_0; z_dot_0];

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
        xlim([-100 100]);
        ylim([-200 200]);
        zlim([-150 150]);
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

%% Task 2 - Analisys of the Linearized Dynamical System
A = [ 0 1 0 0 0 0;
     3*n^2 0 0 2*n 0 0;
     0 0 0 1 0 0;
     0 -2*n 0 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 -n^2 0];
%n non va aggiustata come unità di misura perchè è in [1/s]

eigenpairs = eig(A);

%% Task 3 - Numerical Simulation and Model Comparison
B=zeros(size(A,1), 1);
C=[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 0];
%C è costruita in questo modo perchè non ci interessano le velocità ma solo
%le posizioni in questa simulazione
D=zeros(size(C,1), 1);

sys=ss(A,B,C,D);
t_linear=linspace(0,2*proximityP.T, 10000);

t_limite = NaN(1, size(x0,2));

for i=1:size(x0,2)

    ODE_obj = ode;  
    ODE_obj.ODEFcn = @(t,x) proximityP_f(t,x, proximityP);
    ODE_obj.InitialValue = x0(:,i).*10^3;
    ODE_obj.Solver = 'ode45';
    ODEResults_obj = solve(ODE_obj, t_0, t_f);
    tt = ODEResults_obj.Time'; 
    xx = (ODEResults_obj.Solution');

    [y]=initial(sys, x0(:,i).*10^3, t_linear);
    figure('Name',['Comparison linear vs non linear ', scenari(i)]);
    plot(tt, xx(:,1), 'LineWidth',2);
    hold on;
    plot(tt, xx(:,3), 'LineWidth',2);
    plot(tt, xx(:,5), 'LineWidth',2);
    plot(t_linear, y(:,1),'--','Color','black', LineWidth=2);
    plot(t_linear, y(:,3),'--', LineWidth=2);
    plot(t_linear, y(:,5),'--', LineWidth=2);
    grid on;
    xlabel('Time[s]');
    ylabel('Position[m]');
    legend('x(t)', 'y(t)', 'z(t)', 'x_l(t)', 'y_l(t)', 'z_l(t)');
    title(['Coordinates case ', scenari(i)]);

    y_err_1=interp1(t_linear, y(:,1), tt, 'linear');
    y_err_3=interp1(t_linear, y(:,3), tt, 'linear');
    y_err_5=interp1(t_linear, y(:,5), tt, 'linear');

%plottiamo gli errori
figure('Name','Relative error between linear and not linear')
    plot(tt, (abs(y_err_1(:)-xx(:,1)))./abs(xx(:,1)), LineWidth=2);
    hold on;
    plot(tt, abs(y_err_3(:)-xx(:,3))./abs(xx(:,3)), LineWidth=2);
    plot(tt, abs(y_err_5(:)-xx(:,5))./abs(xx(:,5)), LineWidth=2);
    grid on;
    xlabel('Time[s]');
    ylabel('Relative Error in Position');
    yline(0.01, 'r--', LineWidth=2.5);
    ylim([0 0.02]);
    xlim([0 1000]);
    legend('err. x(t)', 'err. y(t)', 'err. z(t)');
    title(['Coordinates case ', scenari(i)]);

end

%quale può essere un errore accettabile nel confronto tra lineare e non
%lineare?
%sapendo che il target è una nave spaziale, e che il chaser è un
%modulo che deve approcciare la navicella, il margine di
%errore sulla posizione è molto ristretto. nel caso di approccio fra i due
%oggetti, il margine di errore potrebbe essere dell'1% sulla posizione
%relativa, per ogni coordinata. 
% è quello che è stato plottato nel grafico,
%estraiamo a mano circa un tempo limite accettabile per l'approssimazione
%lineare, 



%% task 4 - Analytical Solution Validation and Notable Cases
%una volta risolta l'equazione analiticamente, matlab serve a risolvere il
%problema di cauchy, nel senso che risolvo un sistema lineare per trovare
%le varie condizioni al contorno.
S=[1 1 0 2/n 0 0;
    1i*n, -1i*n, 0 0 0 0;
    2*1i, -2*1i, 1 0 0 0;
    -2*n, -2*n, 0, -3, 0 0;
    0 0 0 0 1 1;
    0 0 0 0 1i*n, -1i*n]; %è in [m]
coeff=zeros(5,6);

coefficienti=zeros(5,6);
for k=1:size(x0,2)
    coeff_curr=S\(x0(:,k).*10^3); %i coefficienti sono messi in ordine da x_0 a z_1, in ordine, in modo da
    %essere direttamente sostituiti nella formula analitica
    coefficienti(k,:)=coeff_curr; %serve per controllo

x=@(t,coeff) coeff(1)*exp(1i*n*t)+coeff(2)*exp(-1i*n*t)+2*(coeff(4)/n);
y=@(t,coeff) coeff(3)-3*coeff(4).*t+2*1i*(coeff(1)*exp(1i*n*t)-coeff(2)*exp(-1i*n*t));
z=@(t,coeff) coeff(5)*exp(1i*n*t)+coeff(6)*exp(-1i*n*t);
t_vect=0:1:2*proximityP.T;
x_vect=real(x(t_vect, coeff_curr));
y_vect=real(y(t_vect, coeff_curr));
z_vect=real(z(t_vect,coeff_curr)); %real assicura che la parte immaginaria venga messa a zero, essendo molto molto piccola

[y]=initial(sys, x0(:,k).*10^3, t_vect);

figure('Name',['Comparison analityc vs non linear ', scenari(k)]);
    plot(t_vect, y(:,1), 'LineWidth',5);
    hold on;
    plot(t_vect, y(:,3), 'LineWidth',5);
    plot(t_vect, y(:,5), 'LineWidth',5);
    plot(t_vect, x_vect,'--','Color','black', LineWidth=5);
    plot(t_vect, y_vect,'--', LineWidth=5);
    plot(t_vect, z_vect,'--', LineWidth=5);
    grid on;
    xlabel('Time[s]');
    ylabel('Position[m]');
    legend('x(t)', 'y(t)', 'z(t)', 'x_a(t)', 'y_a(t)', 'z_a(t)');
    title(['Coordinates case ', scenari(k)]);

end

%% Task 5
% 5.1 - Linear Optimization and Nonlinear Verification
% Homing : avviciniamo il chaser all'obbiettivo ponendolo su un'orbita di
% avvicinamento-
% Partiamo dalle condizioni iniziali D, con il moto che inizia a t = 0
% In questa sezione cerchiamo l'impulso ottimale DeltaV1 per arrivare al target

% 1. Impostazioni del problema
scenario_idx = 4; % Scegliamo lo Scenario D 
scenario_target = xx(:,1);
x0_start = x0(:, scenario_idx) * 10^3; % Stato iniziale in metri [x; vx; y; vy; z; vz]
t_transfer = proximityP.T ; % Tempo di trasferimento (un'orbita)
target_pos = [scenario_target(1); scenario_target(3); scenario_target(5)]; % Obiettivo: origine (posizione relativa zero)

% 2. Definizione delle variabili decisionali
% Cerchiamo le 3 componenti dell'impulso iniziale: u = [dv1x, dv1y, dv1z]
u0 = [2; 2; 2]; % Punto di partenza per l'ottimizzatore (m/s)

% 3. Opzioni dell'algoritmo fmincon
% Usiamo 'sqp' perché è molto robusto per problemi di meccanica orbitale
options = optimoptions('fmincon', 'Display', 'iter-detailed', 'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-9);

% 4. Chiamata a fmincon
% La funzione obiettivo minimizza la somma dei moduli dei due impulsi (partenza + arrivo)
[u_opt, J_min] = fmincon(@(u) objective_dv(u, x0_start, S, n, t_transfer), ...
                         u0, [], [], [], [], [], [], ...
                         @(u) constraints_pos(u, x0_start, S, n, t_transfer, target_pos), ...
                         options);

% 5. Visualizzazione dei risultati
fprintf('\n--- RISULTATI OTTIMIZZAZIONE TASK 5 ---\n');
fprintf('DeltaV1 Ottimale (m/s): [%.4f, %.4f, %.4f]\n', u_opt);
fprintf('Costo Totale (DeltaV1 + DeltaV2): %.4f m/s\n', J_min);



% FUNZIONE OBIETTIVO: Minimizza il DeltaV totale (J = |dv1| + |dv2|)
function J = objective_dv(u, x0_start, S, n, t)
    % u è il DeltaV1 applicato all'istante iniziale
    dv1_mag = norm(u);
    
    % Calcoliamo dove arriveremmo applicando questo impulso
    x0_plus = x0_start;
    x0_plus([2, 4, 6]) = x0_plus([2, 4, 6]) + u(:); % Aggiungo dv1 alle velocità iniziali
    coeff = S \ x0_plus; % Trovo i coefficienti analitici con il nuovo impulso
    
    % Calcolo la velocità finale al tempo t (derivate delle tue formule Task 4)
    vx_f = real(1i*n*coeff(1)*exp(1i*n*t) - 1i*n*coeff(2)*exp(-1i*n*t));
    vy_f = real(-3*coeff(4) - 2*n*coeff(1)*exp(1i*n*t) - 2*n*coeff(2)*exp(-1i*n*t));
    vz_f = real(1i*n*coeff(5)*exp(1i*n*t) - 1i*n*coeff(6)*exp(-1i*n*t));
    
    % Per il rendezvous, all'arrivo dobbiamo azzerare la velocità relativa
    % Quindi dv2 è pari al modulo della velocità residua
    dv2_mag = norm([vx_f; vy_f; vz_f]);
    
    J = dv1_mag + dv2_mag;
end

% VINCOLI: Garantisce che al tempo t_transfer la posizione sia quella target
function [c, ceq] = constraints_pos(u, x0_start, S, n, t, target_pos)
    % Calcolo traiettoria analitica con l'impulso scelto dall'ottimizzatore
    x0_plus = x0_start;
    x0_plus([2, 4, 6]) = x0_plus([2, 4, 6]) + u(:);
    coeff = S \ x0_plus;
    
    % Posizione calcolata con le tue formule analitiche (Task 4)
    xf = real(coeff(1)* exp(1i*n*t) + coeff(2)*exp(-1i*n*t) + 2*(coeff(4)/n));
    yf = real(coeff(3) - 3*coeff(4)*t + 2*1i*(coeff(1)*exp(1i*n*t) - coeff(2)*exp(-1i*n*t)));
    zf = real(coeff(5)*exp(1i*n*t) + coeff(6)*exp(-1i*n*t));
    
    % Vincolo di uguaglianza: Posizione finale - Target = 0
    ceq = [xf; yf; zf] - target_pos; 
    c = []; % Nessun vincolo di disuguaglianza
end

% 6. Calcolo traiettoria ottimale per il plot
% Applico l'impulso ottimale alle condizioni iniziali
x0_plus = x0_start;
x0_plus([2, 4, 6]) = x0_plus([2, 4, 6]) + u_opt(:);
coeff_opt = S \ x0_plus;

% Vettore temporale per la propagazione
t_vec = linspace(0, t_transfer, 1000);

% Allocazione
x_traj = zeros(size(t_vec));
y_traj = zeros(size(t_vec));
z_traj = zeros(size(t_vec));

for k = 1:length(t_vec)
    t_k = t_vec(k);
    x_traj(k) = real(coeff_opt(1)*exp(1i*n*t_k) + coeff_opt(2)*exp(-1i*n*t_k) + 2*(coeff_opt(4)/n));
    y_traj(k) = real(coeff_opt(3) - 3*coeff_opt(4)*t_k + 2*1i*(coeff_opt(1)*exp(1i*n*t_k) - coeff_opt(2)*exp(-1i*n*t_k)));
    z_traj(k) = real(coeff_opt(5)*exp(1i*n*t_k) + coeff_opt(6)*exp(-1i*n*t_k));
end

% Punto di partenza e di arrivo
pos_start = [x_traj(1);   y_traj(1);   z_traj(1)];
pos_end   = [x_traj(end); y_traj(end); z_traj(end)];

% ── 7. PLOT ──────────────────────────────────────────────────────────────────


traj = struct(); % conterrà i campi A e D

for i = 1:size(x0, 2)
    ODE_obj            = ode;
    ODE_obj.ODEFcn     = @(t,x) proximityP_f(t, x, proximityP);
    ODE_obj.InitialValue = x0(:,i) .* 1e3;
    ODE_obj.Solver     = 'ode45';
    ODEResults_obj     = solve(ODE_obj, t_0, t_f);
    tt = ODEResults_obj.Time';
    xx = ODEResults_obj.Solution';

    % ── Salva traiettoria per A (i=1) e D (i=4) ──────────────────────────
    if i == 1
        traj.A.x = xx(:,1); traj.A.y = xx(:,3); traj.A.z = xx(:,5);
    elseif i == 4
        traj.D.x = xx(:,1); traj.D.y = xx(:,3); traj.D.z = xx(:,5);
    end

    % ... il resto del tuo loop rimane identico ...
end

figure('Name', 'Orbite A e D – Confronto 3D', 'NumberTitle', 'off');

plot3(traj.A.x, traj.A.y, traj.A.z, 'b-', 'LineWidth', 2, 'DisplayName', 'Orbita A');
hold on;
plot3(traj.D.x, traj.D.y, traj.D.z, 'r-', 'LineWidth', 2, 'DisplayName', 'Orbita D');

% Punti di partenza
plot3(traj.A.x(1), traj.A.y(1), traj.A.z(1), 'bo', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Start A');
plot3(traj.D.x(1), traj.D.y(1), traj.D.z(1), 'ro', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Start D');

% Origine (target / Chief)
plot3(0, 0, 0, 'k*', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', 'Chief (origine)');

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Orbite Relative – Scenari A e D');
legend('Location', 'best');
grid on; axis equal;
view(35, 25);
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
    n  = proximityP.n;

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    x6 = x(6);

    rc = sqrt((R0+x1)^2 + x3^2 + x5^2);

    xdot_1 = x2;
    xdot_2 = 2*n*x4 + n^2*(R0+x1) - (mu*(R0+x1))/rc^3;
    xdot_3 = x4;
    xdot_4 = -2*n*x2 + n^2*x3-(mu*x3)/(rc^3);
    xdot_5 = x6;
    xdot_6 = -(mu*x5)/rc^3;

    xdot = [xdot_1, xdot_2, xdot_3, xdot_4, xdot_5, xdot_6]';
end