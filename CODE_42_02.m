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
J_best=inf;

for i=1:10
    if i == 1
        [x_opt, J_min] = fmincon(@(u) objective_dv(u, xD_start, coeff_A, S, n), ...
                         x0_guess, [], [], [], [], lb, ub, ...
                         @(u) constraints_pos(u, xD_start, coeff_A, S, n), options);
    %J_min è il costo totale della missione
    end
    if  J_min >= norm (x0_guess(1:3)) && i > 1
        norma  = J_min;
        dir_casuale = randn(3, 1);
    
    % 2. Normalizzo il vettore dividendolo per la sua stessa norma (ora è lungo esattamente 1)
    dir_normalizzata = dir_casuale / norm(dir_casuale);
    
    % 3. Ti riduce rispetto a J_min questo! quindi ottengo delle guess che
    % in modulo sono sempre inferiori a J_MIN
    nuovo_modulo = rand() * norma;
    
    % Moltiplico la direzione di raggio 1 per il nuovo modulo
    nuovo_vettore_3d = dir_normalizzata * nuovo_modulo;
    
    % 4. Ricostruisco il guess iniziale concatenando le parti
    x0_guess = [nuovo_vettore_3d; x_opt(4); x_opt(5); x_opt(6)];

    [x_opt, J_min] = fmincon(@(u) objective_dv(u, xD_start, coeff_A, S, n), ...
                         x0_guess, [], [], [], [], lb, ub, ...
                         @(u) constraints_pos(u, xD_start, coeff_A, S, n), options);
    end
disp(J_min);

    if J_min < J_best %sta roba non succede mai dopo la prima iterazione, succede solo una volta. Questo conferma che il migliore risultato possibile, per 10 test di guess
        %iniziali minori del deltaV totale della prima guess, è quello dato
        %dalla prima guess. tutte queste guess sono quindi guess iniziali
        %ottimali!
        J_best = J_min;
        x_best = x_opt;
    end
end

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

    %[transfer_orbit_linear]=initial(sys, (D_orbit_linear(end,:)+[0;Dvx;0;Dvy;0;Dvz]')', linspace(0, t_tof_opt, 1000));

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
figure('Name','Confronto Modelli 3D: Lineare vs Non Lineare','NumberTitle','off');
hold on; grid on;

% Definizione colori per chiarezza
col_A  = [0.50, 0.50, 0.50]; % Grigio per l'orbita bersaglio
col_NL = [0.85, 0.33, 0.10]; % Rosso-Arancio per il Non Lineare
col_L  = [0.00, 0.45, 0.74]; % Blu per il Lineare

% 1. Plot Orbita Target A (Sfondo)
plot3(pos_A_L(:,1), pos_A_L(:,2), pos_A_L(:,3), ':', 'Color', col_A, 'LineWidth', 1.5, 'DisplayName', 'Orbita A Target');

% 2. Plot Fase di Attesa sull'Orbita D
plot3(pos_D_NL(:,1), pos_D_NL(:,2), pos_D_NL(:,3), '-', 'Color', col_NL, 'LineWidth', 1.5, 'DisplayName', 'Attesa D (Non Lineare)');
plot3(pos_D_L(:,1), pos_D_L(:,2), pos_D_L(:,3), '--', 'Color', col_L, 'LineWidth', 1.5, 'DisplayName', 'Attesa D (Lineare)');

% 3. Plot Fase di Trasferimento
plot3(pos_T_NL(:,1), pos_T_NL(:,2), pos_T_NL(:,3), '-', 'Color', col_NL, 'LineWidth', 2.5, 'DisplayName', 'Transfer (Non Lineare)');
plot3(pos_T_L(1,:), pos_T_L(2,:), pos_T_L(3,:), '--', 'Color', col_L, 'LineWidth', 2.5, 'DisplayName', 'Transfer (Lineare)');

% 4. MARKERS: Punti Notevoli
% Punto di partenza comune (t=0)
plot3(pos_D_L(1,1), pos_D_L(1,2), pos_D_L(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Start (t=0)');

% Punto di applicazione del Delta-V 1
plot3(pos_D_NL(end,1), pos_D_NL(end,2), pos_D_NL(end,3), 'o', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 7, 'DisplayName', 'Impulso NL');
plot3(pos_D_L(end,1), pos_D_L(end,2), pos_D_L(end,3), 'o', 'Color', col_L, 'MarkerFaceColor', col_L, 'MarkerSize', 7, 'DisplayName', 'Impulso Lineare');

% Punto di Arrivo (Dove si vede l'errore calcolato!)
plot3(pos_T_NL(end,1), pos_T_NL(end,2), pos_T_NL(end,3), 's', 'Color', col_NL, 'MarkerFaceColor', col_NL, 'MarkerSize', 10, 'DisplayName', 'Arrivo NL');
plot3(pos_T_L(1,end), pos_T_L(2,end), pos_T_L(3,end), 's', 'Color', col_L, 'MarkerFaceColor', col_L, 'MarkerSize', 10, 'DisplayName', 'Arrivo Lineare');

% Origine (Target)
plot3(0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Origine (Chief)');

% 5. Formattazione del grafico
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Traiettorie di Avvicinamento: Lineare vs Non Lineare');
legend('Location', 'best');
axis equal; 
view(35, 25); % Angolo per visualizzare bene la differenza lungo l'asse Y
hold off;

%%
% % 6. Calcolo traiettorie pulite per il plot
% 
% % 1. Fase di attesa sull'orbita D (da t=0 a t_wait)
% coeff_D = S \ xD_start;
% t_wait_vec = linspace(0, t_wait_opt, 500);
% stati_wait = getState(coeff_D, t_wait_vec, n);
% 
% % 2. Fase di trasferimento (da t_wait a t_wait + t_tof)
% stato_dopo_dv1 = stati_wait(:, end);
% stato_dopo_dv1([2,4,6]) = stato_dopo_dv1([2,4,6]) + dv1_opt(:); % Applico DV1
% coeff_T = S \ stato_dopo_dv1;
% 
% t_tof_vec = linspace(0, t_tof_opt, 500);
% stati_trans = getState(coeff_T, t_tof_vec, n);
% 
% % 3. Orbita A di background (Disegno esattamente 1 periodo per chiarezza)
% t_A_vec = linspace(0, T, 1000);
% stati_A = getState(coeff_A, t_A_vec, n);
% 
% % Punto effettivo di aggancio previsto sull'orbita A
% t_arrivo_assoluto = t_wait_opt + t_tof_opt + tau_opt;
% stato_A_target = getState(coeff_A, t_arrivo_assoluto, n);

% % 7. PLOT (Solo Vista 3D)
% col_A = [0.00, 0.45, 0.74];   % blu (Orbita A)
% col_D = [0.85, 0.33, 0.10];   % arancio-rosso (Orbita D)
% col_T = [0.47, 0.67, 0.19];   % verde (Trasferimento)
% 
% figure('Name','D->A: Vista 3D','NumberTitle','off');
% 
% % Usiamo direttamente Ax, Ay, Az dal tuo workspace per ripristinare l'Orbita A
% plot3(Ax, Ay, Az, '--', 'Color', col_A, 'LineWidth', 1.5, 'DisplayName', 'Orbita A (Target)'); hold on;
% 
% % Plotto l'arco di attesa su D
% plot3(stati_wait(1,:), stati_wait(3,:), stati_wait(5,:), '--', 'Color', col_D, 'LineWidth', 2, 'DisplayName', 'Attesa su Orbita D');
% 
% % Plotto la traiettoria di trasferimento
% plot3(stati_trans(1,:), stati_trans(3,:), stati_trans(5,:), '-', 'Color', col_T, 'LineWidth', 2.5, 'DisplayName', 'Arco trasferimento');
% 
% % Punti notevoli
% plot3(xD_start(1), xD_start(3), xD_start(5), 'd', 'Color', col_D, 'MarkerSize',8, 'MarkerFaceColor', col_D, 'DisplayName', 't=0');
% plot3(stati_wait(1,end), stati_wait(3,end), stati_wait(5,end), 'o', 'Color', col_D, 'MarkerSize',10, 'MarkerFaceColor', col_D, 'DisplayName', 'Partenza (dv1)');
% plot3(stati_trans(1,end), stati_trans(3,end), stati_trans(5,end), 's', 'Color', col_T, 'MarkerSize',10, 'MarkerFaceColor', col_T, 'DisplayName', 'Aggancio (dv2)');
% plot3(stato_A_target(1), stato_A_target(3), stato_A_target(5), 'k*', 'MarkerSize', 14, 'LineWidth', 1.5, 'DisplayName', 'Punto di Inserimento Ottimo su A');
% plot3(0, 0, 0, 'k+', 'MarkerSize',10, 'LineWidth',2, 'DisplayName', 'Chief (origine)');
% 
% xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
% title(sprintf('Inserimento Ottimale D \\rightarrow A (T_{attesa} = %.0fs, T_{volo} = %.0fs)', t_wait_opt, t_tof_opt));
% legend('Location','best'); grid on; axis equal; view(35,25);

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