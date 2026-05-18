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
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];
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
x0_guess = [0.1; 0.1; 0; T/2; T/2; 0]; 

options = optimoptions('fmincon','Algorithm', 'interior-point', 'OptimalityTolerance', 1e-9);

% 3. Limiti (Lower e Upper Bounds)
% dv compresi tra -10 e 10
% t_wait >= T/3
% t_tof compreso tra 0 e T
% tau compreso tra -T e T (libertà totale di agganciarsi ovunque sull'orbita)
lb = [-10; -10; -10; T/3; 0.1; -T];
ub = [ 10;  10;  10; Inf;   T;  T];

% 4. Chiamata fmincon
xA_start = x0(:, 1) * 10^3; 
coeff_A = S \ xA_start;     % Coefficienti dell'orbita bersaglio A

% Creazione del problema per l'ottimizzatore globale
problem = createOptimProblem('fmincon', ...
    'objective', ...
    @(u) objective_dv(u, xD_start, coeff_A, S, n), ...
    'x0', x0_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'nonlcon', @(u) constraints_pos(u, xD_start, coeff_A, S, n), ...
    'options', options);

% Configurazione e Lancio di GlobalSearch
% (NumTrialPoints = 200 è un buon compromesso per non rallentare troppo il pc)
gs = GlobalSearch('Display', 'iter', 'NumTrialPoints', 200);
[x_opt, J_min] = run(gs, problem);

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

%applichiamo i risultati al sistema non lineare
%- diciamo a ODE di integrare normale fino all'istante del primo impulso
%- usiamo le condizioni iniziali in quel punto+nostro DV per calcolare la
%traiettoria fino a t_of_flight
%- calcoliamo l'errore con il confronto con il punto di arrivo lineare vs
%non lineare.

ODE_obj = ode;  
ODE_obj.ODEFcn = @(t,x) proximityP_f(t,x, proximityP);
ODE_obj.InitialValue = x0(:,4).*10^3;
ODE_obj.Solver = 'ode45';
ODEResults_obj = solve(ODE_obj, t_0, t_wait_opt);
tt = ODEResults_obj.Time'; 
D_orbit = (ODEResults_obj.Solution');

Dvx=x_opt(1);
Dvy=x_opt(2);
Dvz=x_opt(3);

ODE_obj = ode;  
ODE_obj.ODEFcn = @(t,x) proximityP_f(t,x, proximityP);
ODE_obj.InitialValue = (D_orbit(end,:)+[0;Dvx;0;Dvy;0;Dvz]')';
ODE_obj.Solver = 'ode45';
ODEResults_obj = solve(ODE_obj, 0, t_tof_opt);
t_trasferimento = ODEResults_obj.Time'; 
transfer_orbit = (ODEResults_obj.Solution');

[D_orbit_linear]=initial(sys, x0(:,4).*10^3, linspace(0, t_wait_opt, 1000));

% 1. Fase di attesa sull'orbita D (da t=0 a t_wait)
coeff_D = S \ xD_start;
t_wait_vec = linspace(0, t_wait_opt, 500);
stati_wait = getState(coeff_D, t_wait_vec, n);

%2. Fase di trasferimento (da t_wait a t_wait + t_tof)
stato_dopo_dv1 = stati_wait(:, end);
stato_dopo_dv1([2,4,6]) = stato_dopo_dv1([2,4,6]) + [Dvx;Dvy;Dvz]; % Applico DV1
coeff_T = S \ stato_dopo_dv1;
t_tof_vec = linspace(0, t_tof_opt, 500);
stati_trans = getState(coeff_T, t_tof_vec, n);

delta_pos=transfer_orbit(end,[1,3,5])-stati_trans([1,3,5],end)';
delta_vel=transfer_orbit(end,[2,4,6])-stati_trans([2,4,6],end)';
disp(delta_pos);
disp(delta_vel);

% 6. Preparazione Dati per il Plot
% Estraiamo le posizioni (x, y, z) dai risultati calcolati al punto 5.
% Ricordiamo che le posizioni sono le colonne 1, 3, 5.

% -- Modello NON Lineare --
pos_D_NL = D_orbit(:, [1, 3, 5]);
pos_T_NL = transfer_orbit(:, [1, 3, 5]);

% -- Modello Lineare --
pos_D_L = D_orbit_linear(:, [1, 3, 5]);
pos_T_L = stati_trans([1, 3, 5],:);

% -- Orbita Target (A) di riferimento --
% Ricalcoliamo 1 periodo dell'orbita bersaglio per contesto visivo
t_A_vec = linspace(0, T, 1000);
stati_A = getState(coeff_A, t_A_vec, n);
pos_A_L = stati_A([1, 3, 5], :)'; % Trasposta per avere una matrice N x 3

% 7. PLOT 3D: Confronto Lineare vs Non Lineare
figure('Name','Comparison between Linear and Nonlinear Models','NumberTitle','off');
hold on; grid on;
% Definizione colori per chiarezza
col_A  = [0.50, 0.50, 0.50]; % Grigio per l'orbita bersaglio
col_NL = [0.85, 0.33, 0.10]; % Rosso-Arancio per il Non Lineare
col_L  = [0.00, 0.45, 0.74]; % Blu per il Lineare
% 1. Plot Orbita Target A (Sfondo)
plot3(pos_A_L(:,1), pos_A_L(:,2), pos_A_L(:,3), ':', 'Color', col_A, 'LineWidth', 1.5, 'DisplayName', 'Orbit A Target');
% 2. Plot Fase di Attesa sull'Orbita D
plot3(pos_D_NL(:,1), pos_D_NL(:,2), pos_D_NL(:,3), '-', 'Color', col_NL, 'LineWidth', 1.5, 'DisplayName', 'Waiting time D (Non Linear)');
plot3(pos_D_L(:,1), pos_D_L(:,2), pos_D_L(:,3), '--', 'Color', col_L, 'LineWidth', 1.5, 'DisplayName', 'Waiting time D (Linear)');
% 3. Plot Fase di Trasferimento
plot3(pos_T_NL(:,1), pos_T_NL(:,2), pos_T_NL(:,3), '-', 'Color', col_NL, 'LineWidth', 2.5, 'DisplayName', 'Transfer (Non Linear)');
plot3(pos_T_L(1,:), pos_T_L(2,:), pos_T_L(3,:), '--', 'Color', col_L, 'LineWidth', 2.5, 'DisplayName', 'Transfer (Linear)');
% 4. MARKERS: Punti Notevoli
% Punto di partenza comune (t=0)
plot3(pos_D_L(1,1), pos_D_L(1,2), pos_D_L(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Start (t=0)');
% Punto di applicazione del Delta-V 1
plot3(pos_D_NL(end,1), pos_D_NL(end,2), pos_D_NL(end,3), 'o', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 7, 'DisplayName', 'Impulse NL');
plot3(pos_D_L(end,1), pos_D_L(end,2), pos_D_L(end,3), 'o', 'Color', col_L, 'MarkerFaceColor', col_L, 'MarkerSize', 7, 'DisplayName', 'Impulse Linear');
% Punto di Arrivo (Dove si vede l'errore calcolato!)
plot3(pos_T_NL(end,1), pos_T_NL(end,2), pos_T_NL(end,3), 's', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 10, 'DisplayName', 'Arrival NL');
plot3(pos_T_L(1,end), pos_T_L(2,end), pos_T_L(3,end), 's', 'Color', col_L, 'MarkerFaceColor', col_L, 'MarkerSize', 10, 'DisplayName', 'Arrival Linear');
% Origine (Target)
plot3(0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Origin (Chief)');
% 5. Formattazione del grafico
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Comparison between Linear and Nonlinear Models');
legend('Location', 'best');
axis equal; 
% Io farei vista (70, 50), vista (120,30) e (60,10)
view(70,50); % Angolo per visualizzare bene la differenza lungo l'asse Y
hold off;
ylim([-900, 300]);
xlim([-500, 300]);

% Grafico della sola traiettoria prevista dal modello non lineare
figure('Name','Nonlinear transfer trajectory','NumberTitle','off');
hold on; grid on;
% 1. Plot Orbita A (sfondo)
plot3(pos_A_L(:,1), pos_A_L(:,2), pos_A_L(:,3), ':', 'Color', col_A, 'LineWidth', 1.5, 'DisplayName', 'Orbit A Target');
% 2. Plot Fase di Attesa sull'Orbita D
plot3(pos_D_NL(:,1), pos_D_NL(:,2), pos_D_NL(:,3), '-', 'Color', col_NL, 'LineWidth', 1.5, 'DisplayName', 'Waiting Time D');
% 3. Plot Fase di Trasferimento
plot3(pos_T_NL(:,1), pos_T_NL(:,2), pos_T_NL(:,3), '-', 'Color', col_NL, 'LineWidth', 2.5, 'DisplayName', ['Transfer']);
% 4. MARKERS: Punti Notevoli
% Punto di partenza (t=0)
plot3(pos_D_NL(1,1), pos_D_NL(1,2), pos_D_NL(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Start (t=0)');
% Punto di applicazione del Delta-V 1
plot3(pos_D_NL(end,1), pos_D_NL(end,2), pos_D_NL(end,3), 'o', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 7, 'DisplayName', 'Impulse');
% Punto di Arrivo (Dove si vede l'errore calcolato!)
plot3(pos_T_NL(end,1), pos_T_NL(end,2), pos_T_NL(end,3), 's', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 10, 'DisplayName', ['Arrival']);
% Origine (Target)
plot3(0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Origin (Chief)');
% 5. Formattazione del grafico
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Nonlinear transfer trajectory');
legend('Location', 'best');
axis equal; 
view(30, 40); % Angolo per visualizzare bene la differenza lungo l'asse Y
hold off;
ylim([-900, 300]);
xlim([-500, 300]);

% figure('Name','3D Real-Time Trajectory Animation','NumberTitle','off');
% hold on; grid on; axis equal; 
% view(70, 50); % Stesso angolo di visualizzazione del tuo primo plot
% xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
% title('3D Animation: Linear vs Nonlinear Proximity Operations');
% 
% % 1. Disegno dello sfondo statico (Orbita Target A e Origine)
% plot3(pos_A_L(:,1), pos_A_L(:,2), pos_A_L(:,3), ':', 'Color', col_A, 'LineWidth', 1.5, 'Handlevisibility', 'off');
% plot3(0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Origin (Chief)');
% 
% % 2. Inizializzazione delle linee animate (Scia dei satelliti)
% an_D_NL = animatedline('Color', col_NL, 'LineWidth', 1.5, 'LineStyle', '-',  'DisplayName', 'Nonlinear Path');
% an_T_NL = animatedline('Color', col_NL, 'LineWidth', 2.5, 'LineStyle', '-');
% an_D_L  = animatedline('Color', col_L,  'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Linear Path');
% an_T_L  = animatedline('Color', col_L,  'LineWidth', 2.5, 'LineStyle', '--');
% 
% % 3. Inizializzazione dei Marker per i satelliti in movimento
% sat_NL = plot3(NaN, NaN, NaN, 'o', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 6, 'DisplayName', 'Sat NL');
% sat_L  = plot3(NaN, NaN, NaN, 'o', 'Color', col_L,  'MarkerFaceColor', col_L,  'MarkerSize', 6, 'DisplayName', 'Sat Linear');
% 
% legend('Location', 'best');
% 
% % Configurazione campionamento animazione (Numero di fotogrammi per fase)
% N_frames = 200; 
% pause_time = 0.01; % Minore = più veloce, Maggiore = più lento
% 
% t_wait_common = linspace(0, t_wait_opt, N_frames);
% pos_D_NL_interp = interp1(tt(:), pos_D_NL, t_wait_common(:));
% 
% % SOLUZIONE: Genera il vettore dei tempi basandoti sul numero reale di righe di pos_D_L (1000 punti)
% t_D_L_real = linspace(0, t_wait_opt, size(pos_D_L, 1));
% pos_D_L_interp  = interp1(t_D_L_real(:), pos_D_L, t_wait_common(:));
% 
% % Punto di partenza iniziale (t=0)
% plot3(pos_D_L(1,1), pos_D_L(1,2), pos_D_L(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'Handlevisibility', 'off');
% 
% for i = 1:N_frames
%     % Aggiorna le scie (Fase Attesa)
%     addpoints(an_D_NL, pos_D_NL_interp(i,1), pos_D_NL_interp(i,2), pos_D_NL_interp(i,3));
%     addpoints(an_D_L,  pos_D_L_interp(i,1),  pos_D_L_interp(i,2),  pos_D_L_interp(i,3));
% 
%     % Muovi i satelliti
%     set(sat_NL, 'XData', pos_D_NL_interp(i,1), 'YData', pos_D_NL_interp(i,2), 'ZData', pos_D_NL_interp(i,3));
%     set(sat_L,  'XData', pos_D_L_interp(i,1),  'YData', pos_D_L_interp(i,2),  'ZData', pos_D_L_interp(i,3));
% 
%     drawnow;
%     pause(pause_time);
% end
% % --- FASE 2: ANIMAZIONE TRASFERIMENTO (Transfer) ---
% t_tof_common = linspace(0, t_tof_opt, N_frames);
% pos_T_NL_interp = interp1(t_trasferimento(:), pos_T_NL, t_tof_common(:));
% pos_T_L_interp  = interp1(t_tof_vec(:), pos_T_L', t_tof_common(:)); % Trasposta per match dimensioni
% 
% for i = 1:N_frames
%     % Aggiorna le scie (Fase Trasferimento)
%     addpoints(an_T_NL, pos_T_NL_interp(i,1), pos_T_NL_interp(i,2), pos_T_NL_interp(i,3));
%     addpoints(an_T_L,  pos_T_L_interp(i,1),  pos_T_L_interp(i,2),  pos_T_L_interp(i,3));
% 
%     % Muovi i satelliti
%     set(sat_NL, 'XData', pos_T_NL_interp(i,1), 'YData', pos_T_NL_interp(i,2), 'ZData', pos_T_NL_interp(i,3));
%     set(sat_L,  'XData', pos_T_L_interp(i,1),  'YData', pos_T_L_interp(i,2),  'ZData', pos_T_L_interp(i,3));
% 
%     drawnow;
%     pause(pause_time);
% end
% 
% % Marker finali di Arrivo
% plot3(pos_T_NL(end,1), pos_T_NL(end,2), pos_T_NL(end,3), 's', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 10, 'Handlevisibility', 'off');
% plot3(pos_T_L(1,end),  pos_T_L(2,end),  pos_T_L(3,end),  's', 'Color', col_L,  'MarkerFaceColor', col_L,  'MarkerSize', 10, 'Handlevisibility', 'off');
% 
% hold off;


% 5.3 Advanced Trajectory Design.

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

function J = objective_dv(u, xD_start, coeff_A, S, n) 
%minimizza il DV totale considerando il boost iniziale e il boost finale per entrare in traiettoria
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
% VINCOLO: La posizione del chaser e quella del punto bersaglio devono coincidere all'arrivo
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

    ceq = stato_arrivo([1,3,5]) - stato_A([1,3,5]);
    c = [];
end

function state = getState(coeff, t, n) %funzione che calcola velocità e posizione a partire da t (vettore o scalare)
    x = real(coeff(1)*exp(1i*n*t) + coeff(2)*exp(-1i*n*t) + 2*(coeff(4)/n));
    y = real(coeff(3) - 3*coeff(4)*t + 2*1i*(coeff(1)*exp(1i*n*t) - coeff(2)*exp(-1i*n*t)));
    z = real(coeff(5)*exp(1i*n*t) + coeff(6)*exp(-1i*n*t));
    
    vx = real(1i*n*coeff(1)*exp(1i*n*t) - 1i*n*coeff(2)*exp(-1i*n*t));
    vy = real(-3*coeff(4) - 2*n*coeff(1)*exp(1i*n*t) - 2*n*coeff(2)*exp(-1i*n*t));
    vz = real(1i*n*coeff(5)*exp(1i*n*t) - 1i*n*coeff(6)*exp(-1i*n*t));
    
    state = [x; vx; y; vy; z; vz];
end