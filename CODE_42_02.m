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

if k==1
    Ax=x_vect;
    Ay=y_vect;
    Az=z_vect;
elseif k==4
    Dx=x_vect;
    Dy=y_vect;
    Dz=z_vect;
end
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

% %% Task 5
% % 5.1 - Linear Optimization and Nonlinear Verification
% % Homing : avviciniamo il chaser all'obbiettivo ponendolo su un'orbita di
% % avvicinamento
% % Partiamo dalle condizioni iniziali D, con il moto che inizia a t = 0
% % In questa sezione cerchiamo l'impulso ottimale DeltaV1 per arrivare al target
% 
% % 1. Impostazioni del problema
% scenario_idx = 4;            % Scegliamo lo Scenario D 
% T=proximityP.T;
% t_transfer = T;       % Tempo di trasferimento massimo
% xD_start = x0(:, scenario_idx) * 10^3;
% % Stato iniziale in metri [x; vx; y; vy; z; vz]
% target_pos = [Ax(1:t_transfer); Ay(1:t_transfer); Az(1:t_transfer)]; % Obiettivo: un punto sull'orbita di A 
% 
% 
% % 2. Definizione delle variabili decisionali
% % Cerchiamo le 3 componenti dell'impulso iniziale: u = [dv1x, dv1y, dv1z]
% u0 = [0; 0; 0; T/2]; % Punto di partenza per l'ottimizzatore (m/s), A CAZZO 
% % Fare algoritmo per cercare la guess-iniziale perfetta
% 
% % 3. Opzioni dell'algoritmo fmincon
% options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
%     'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-9);
% 
% % Limiti Inferiori (Lower Bounds)
% % dv_x >= -10, dv_y >= -10, dv_z >= -10, t_inj >= T/3, t_tof >= 0 (o un valore minimo > 0)
% lb = [-10, -10, -10, T/3];
% 
% % Limiti Superiori (Upper Bounds)
% % dv_x <= 10, dv_y <= 10, dv_z <= 10, t_inj <= Inf, t_tof <= T
% ub = [10, 10, 10, T];
% 
% % 4. Chiamata a fmincon
% % La funzione obiettivo minimizza la somma dei moduli dei due impulsi (partenza + arrivo)
% [u_opt, J_min] = fmincon(@(u) objective_dv(u, xD_start, coeff(1,:), S, proximityP), u0, [], [], [], [], lb, ub, ...
%                          @(u) constraints_pos(u, xD_start, coeff(1,:), S, proximityP), options);
% 
% % 5. Visualizzazione dei risultati
% fprintf('\n--- RISULTATI OTTIMIZZAZIONE TASK 5 ---\n');
% fprintf('DeltaV1 Ottimale (m/s): [%.4f, %.4f, %.4f]\n', u_opt(1:3));
% fprintf('Costo Totale (DeltaV1 + DeltaV2): %.4f m/s\n', J_min);
% 
% dv1_ottimo = u_opt(1:3);
% tempo_di_volo_ottimo = u_opt(4);
% 
% fprintf('Costo totale minimo: %.4f m/s\n', J_min);
% fprintf('Tempo di volo ottimale: %.2f sec (che equivale a %.2f T)\n', tempo_di_volo_ottimo, tempo_di_volo_ottimo/T);

%% Task 5
% 1. Impostazioni del problema
scenario_idx = 4;
T = proximityP.T; 
n = proximityP.n;
xD_start = x0(:, scenario_idx) * 10^3;

% 2. Variabili decisionali
% x = [dv1x, dv1y, dv1z, t_wait, t_tof, tau]
% t_wait: tempo di permanenza sull'orbita D prima di accendere i motori
% t_tof:  tempo di volo dell'arco di trasferimento
% tau:    offset di fase (sceglie il punto ottimale sull'orbita A in cui agganciarsi)

% Punto di partenza per l'ottimizzatore 
% (T/2 di attesa, T/2 di volo, 0 shift di fase)
x0_guess = [3; 2.5; 0; T/2; T/2; 0]; 

options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-9);

% 3. Limiti (Lower e Upper Bounds)
% dv compresi tra -10 e 10
% t_wait >= T/3
% t_tof compreso tra 0 e T
% tau compreso tra -T e T (libertà totale di agganciarsi ovunque sull'orbita)
lb = [-10; -10; -10; T/3; 0.1; -T];
ub = [ 10;  10;  10; Inf;   T;  T];

% 4. Chiamata fmincon
xA_start = x0(:, 1) * 10^3; 
coeff_A = S \ xA_start; % Coefficienti dell'orbita bersaglio A
[x_opt, J_min] = fmincon(@(u) objective_dv(u, xD_start, coeff_A, S, n), ...
                         x0_guess, [], [], [], [], lb, ub, ...
                         @(u) constraints_pos(u, xD_start, coeff_A, S, n), options);

% 5. Estrazione Risultati
dv1_opt    = x_opt(1:3);
t_wait_opt = x_opt(4);
t_tof_opt  = x_opt(5);
tau_opt    = x_opt(6);

fprintf('\n--- RISULTATI OTTIMIZZAZIONE TASK 5 ---\n');
fprintf('DeltaV1 Ottimale (m/s) : [%.4f, %.4f, %.4f]\n', dv1_opt);
fprintf('Tempo di attesa su D   : %.2f sec (%.2f T) [Requisito: >= 0.33 T]\n', t_wait_opt, t_wait_opt/T);
fprintf('Tempo di trasferimento : %.2f sec (%.2f T) [Requisito: <= 1.00 T]\n', t_tof_opt, t_tof_opt/T);
fprintf('Costo totale minimo    : %.4f m/s\n', J_min);

% 6. Calcolo traiettorie pulite per il plot

% 1. Fase di attesa sull'orbita D (da t=0 a t_wait)
coeff_D = S \ xD_start;
t_wait_vec = linspace(0, t_wait_opt, 500);
stati_wait = getState(coeff_D, t_wait_vec, n);

% 2. Fase di trasferimento (da t_wait a t_wait + t_tof)
stato_dopo_dv1 = stati_wait(:, end);
stato_dopo_dv1([2,4,6]) = stato_dopo_dv1([2,4,6]) + dv1_opt(:); % Applico DV1
coeff_T = S \ stato_dopo_dv1;

t_tof_vec = linspace(0, t_tof_opt, 500);
stati_trans = getState(coeff_T, t_tof_vec, n);

% 3. Orbita A di background (Disegno esattamente 1 periodo per chiarezza)
t_A_vec = linspace(0, T, 1000);
stati_A = getState(coeff_A, t_A_vec, n);

% Punto effettivo di aggancio previsto sull'orbita A
t_arrivo_assoluto = t_wait_opt + t_tof_opt + tau_opt;
stato_A_target = getState(coeff_A, t_arrivo_assoluto, n);

% 7. PLOT (Solo Vista 3D)
col_A = [0.00, 0.45, 0.74];   % blu (Orbita A)
col_D = [0.85, 0.33, 0.10];   % arancio-rosso (Orbita D)
col_T = [0.47, 0.67, 0.19];   % verde (Trasferimento)

figure('Name','D->A: Vista 3D','NumberTitle','off');

% Usiamo direttamente Ax, Ay, Az dal tuo workspace per ripristinare l'Orbita A
plot3(Ax, Ay, Az, '--', 'Color', col_A, 'LineWidth', 1.5, 'DisplayName', 'Orbita A (Target)'); hold on;

% Plotto l'arco di attesa su D
plot3(stati_wait(1,:), stati_wait(3,:), stati_wait(5,:), '--', 'Color', col_D, 'LineWidth', 2, 'DisplayName', 'Attesa su Orbita D');

% Plotto la traiettoria di trasferimento
plot3(stati_trans(1,:), stati_trans(3,:), stati_trans(5,:), '-', 'Color', col_T, 'LineWidth', 2.5, 'DisplayName', 'Arco trasferimento');

% Punti notevoli
plot3(xD_start(1), xD_start(3), xD_start(5), 'd', 'Color', col_D, 'MarkerSize',8, 'MarkerFaceColor', col_D, 'DisplayName', 't=0');
plot3(stati_wait(1,end), stati_wait(3,end), stati_wait(5,end), 'o', 'Color', col_D, 'MarkerSize',10, 'MarkerFaceColor', col_D, 'DisplayName', 'Partenza (dv1)');
plot3(stati_trans(1,end), stati_trans(3,end), stati_trans(5,end), 's', 'Color', col_T, 'MarkerSize',10, 'MarkerFaceColor', col_T, 'DisplayName', 'Aggancio (dv2)');
plot3(stato_A_target(1), stato_A_target(3), stato_A_target(5), 'k*', 'MarkerSize', 14, 'LineWidth', 1.5, 'DisplayName', 'Punto di Inserimento Ottimo su A');
plot3(0, 0, 0, 'k+', 'MarkerSize',10, 'LineWidth',2, 'DisplayName', 'Chief (origine)');

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('Inserimento Ottimale D \\rightarrow A (T_{attesa} = %.0fs, T_{volo} = %.0fs)', t_wait_opt, t_tof_opt));
legend('Location','best'); grid on; axis equal; view(35,25);

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

% % FUNZIONE OBIETTIVO: Minimizza il DeltaV totale (J = |dv1| + |dv2|)
% function J = objective_dv(x, x0_start, coeff_A, S, proximityP)
%     % x(1:3) è il vettore DeltaV1 [dvx, dvy, dvz]
%     % x(4) è il tempo di trasferimento t
%     n=proximityP.n;
%     u = x(1:3); % Estraggo il deltaV
%     t = x(4);   % Estraggo il tempo scelto da fmincon per questa iterazione
% 
%     dv1_mag = norm(u);
% 
%     % Calcoliamo dove arriveremmo applicando questo impulso
%     x0_plus = x0_start;
%     x0_plus([2, 4, 6]) = x0_plus([2, 4, 6]) + u(:); 
%     coeff = S \ x0_plus;                            
% 
%     % Calcolo la velocità finale al tempo t 
%     vx_f = real(1i*n*coeff(1)*exp(1i*n*t) - 1i*n*coeff(2)*exp(-1i*n*t));
%     vy_f = real(-3*coeff(4) - 2*n*coeff(1)*exp(1i*n*t) - 2*n*coeff(2)*exp(-1i*n*t));
%     vz_f = real(1i*n*coeff(5)*exp(1i*n*t) - 1i*n*coeff(6)*exp(-1i*n*t));
% 
%     vx_Af = real(1i*n*coeff_A(1)*exp(1i*n*t) - 1i*n*coeff_A(2)*exp(-1i*n*t));
%     vy_Af = real(-3*coeff_A(4) - 2*n*coeff_A(1)*exp(1i*n*t) - 2*n*coeff_A(2)*exp(-1i*n*t));
%     vz_Af = real(1i*n*coeff_A(5)*exp(1i*n*t) - 1i*n*coeff_A(6)*exp(-1i*n*t));
% 
%     % Vettore differenza di velocità (già corretto!)
%     dv2_mag = norm([vx_f - vx_Af; vy_f - vy_Af; vz_f - vz_Af]);
% 
%     J = dv1_mag + dv2_mag;
% end
% 
% % VINCOLI: Garantisce che al tempo t_transfer la posizione sia quella target
% function [c, ceq] = constraints_pos(u, x0_start, coeff_A, S, proximityP)
%     u_in = u(1:3);
%     t = u(4);
%     n=proximityP.n;
%     % Propagazione Chaser
%     x0_plus = x0_start;
%     x0_plus([2, 4, 6]) = x0_plus([2, 4, 6]) + u_in(:);
%     coeff = S \ x0_plus;
% 
%     xf = real(coeff(1)*exp(1i*n*t) + coeff(2)*exp(-1i*n*t) + 2*(coeff(4)/n));
%     yf = real(coeff(3) - 3*coeff(4)*t + 2*1i*(coeff(1)*exp(1i*n*t) - coeff(2)*exp(-1i*n*t)));
%     zf = real(coeff(5)*exp(1i*n*t) + coeff(6)*exp(-1i*n*t));
% 
%     % Calcolo posizione Target A al tempo t
%     x_Af = real(coeff_A(1)*exp(1i*n*t) + coeff_A(2)*exp(-1i*n*t) + 2*(coeff_A(4)/n));
%     y_Af = real(coeff_A(3) - 3*coeff_A(4)*t + 2*1i*(coeff_A(1)*exp(1i*n*t) - coeff_A(2)*exp(-1i*n*t)));
%     z_Af = real(coeff_A(5)*exp(1i*n*t) + coeff_A(6)*exp(-1i*n*t));
% 
%     % Vincolo: Posizione Chaser == Posizione Target
%     ceq = [xf - x_Af; yf - y_Af; zf - z_Af];
%     c = [];
% end

function J = objective_dv(u, xD_start, coeff_A, S, n)
    dv1    = u(1:3);
    t_wait = u(4);
    t_tof  = u(5);
    tau    = u(6);
    
    % 1. Propago l'orbita D fino al tempo di accensione
    coeff_D = S \ xD_start;
    stato_partenza = getState(coeff_D, t_wait, n);
    
    % 2. Applico DV1 alle velocità
    stato_dopo_dv1 = stato_partenza;
    stato_dopo_dv1([2,4,6]) = stato_dopo_dv1([2,4,6]) + dv1(:);
    
    % 3. Trovo coefficienti trasferimento e propago per il tempo di volo
    coeff_T = S \ stato_dopo_dv1;
    stato_arrivo = getState(coeff_T, t_tof, n);
    
    % 4. Valuto la velocità dell'orbita A nel punto di inserimento
    stato_A = getState(coeff_A, t_wait + t_tof + tau, n);
    
    % 5. Il DV2 è la differenza di velocità per restare su A
    dv2_mag = norm(stato_arrivo([2,4,6]) - stato_A([2,4,6]));
    
    J = norm(dv1) + dv2_mag;
end

function [c, ceq] = constraints_pos(u, xD_start, coeff_A, S, n)
    dv1    = u(1:3);
    t_wait = u(4);
    t_tof  = u(5);
    tau    = u(6);
    
    % Stessa identica propagazione
    coeff_D = S \ xD_start;
    stato_partenza = getState(coeff_D, t_wait, n);
    
    stato_dopo_dv1 = stato_partenza;
    stato_dopo_dv1([2,4,6]) = stato_dopo_dv1([2,4,6]) + dv1(:);
    
    coeff_T = S \ stato_dopo_dv1;
    stato_arrivo = getState(coeff_T, t_tof, n);
    
    % Valuto la posizione dell'orbita A nel punto di inserimento
    stato_A = getState(coeff_A, t_wait + t_tof + tau, n);
    
    % VINCOLO: La posizione del chaser e quella del punto bersaglio devono coincidere
    ceq = stato_arrivo([1,3,5]) - stato_A([1,3,5]);
    c = [];
end

function state = getState(coeff, t, n)
    % Questa utility valuta analiticamente stato e velocità
    % Funziona anche se 't' è un vettore (per il plot!)
    x = real(coeff(1)*exp(1i*n*t) + coeff(2)*exp(-1i*n*t) + 2*(coeff(4)/n));
    y = real(coeff(3) - 3*coeff(4)*t + 2*1i*(coeff(1)*exp(1i*n*t) - coeff(2)*exp(-1i*n*t)));
    z = real(coeff(5)*exp(1i*n*t) + coeff(6)*exp(-1i*n*t));
    
    vx = real(1i*n*coeff(1)*exp(1i*n*t) - 1i*n*coeff(2)*exp(-1i*n*t));
    vy = real(-3*coeff(4) - 2*n*coeff(1)*exp(1i*n*t) - 2*n*coeff(2)*exp(-1i*n*t));
    vz = real(1i*n*coeff(5)*exp(1i*n*t) - 1i*n*coeff(6)*exp(-1i*n*t));
    
    state = [x; vx; y; vy; z; vz];
end