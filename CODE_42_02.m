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
idx_D = 4;  idx_A = 1;
x0_D  = x0(:, idx_D) * 1e3;   % condizioni iniziali orbita D [m, m/s]
x0_A  = x0(:, idx_A) * 1e3;   % condizioni iniziali orbita A [m, m/s]

t_transfer = T;                % tempo di trasferimento = 1 periodo

% Propagazione ANALITICA di A fino a t_transfer
% (serve il target di posizione E velocità all'arrivo)
coeff_A = S \ x0_A;

xA_f  = real(coeff_A(1)*exp(1i*n*t_transfer) + coeff_A(2)*exp(-1i*n*t_transfer) ...
           + 2*(coeff_A(4)/n));
yA_f  = real(coeff_A(3) - 3*coeff_A(4)*t_transfer ...
           + 2*1i*(coeff_A(1)*exp(1i*n*t_transfer) - coeff_A(2)*exp(-1i*n*t_transfer)));
zA_f  = real(coeff_A(5)*exp(1i*n*t_transfer) + coeff_A(6)*exp(-1i*n*t_transfer));

vxA_f = real( 1i*n*coeff_A(1)*exp(1i*n*t_transfer) - 1i*n*coeff_A(2)*exp(-1i*n*t_transfer));
vyA_f = real(-3*coeff_A(4) - 2*n*coeff_A(1)*exp(1i*n*t_transfer) ...
           - 2*n*coeff_A(2)*exp(-1i*n*t_transfer));
vzA_f = real( 1i*n*coeff_A(5)*exp(1i*n*t_transfer) - 1i*n*coeff_A(6)*exp(-1i*n*t_transfer));

target_pos = [xA_f; yA_f; zA_f];
target_vel = [vxA_f; vyA_f; vzA_f];

% ── 3. fmincon ───────────────────────────────────────────────────────────
u0      = [2; 2; 2];
options = optimoptions('fmincon', ...
    'Display',              'iter-detailed', ...
    'Algorithm',            'interior-point', ...
    'OptimalityTolerance',  1e-9);

[u_opt, J_min] = fmincon( ...
    @(u) objective_transfer(u, x0_D, S, n, t_transfer, target_vel), ...
    u0, [], [], [], [], [], [], ...
    @(u) constraints_transfer(u, x0_D, S, n, t_transfer, target_pos), ...
    options);

% ── 4. Risultati a schermo ────────────────────────────────────────────────
fprintf('\n=== TRASFERIMENTO D → A ===\n');
fprintf('DeltaV1      [m/s]: [%+.4f, %+.4f, %+.4f]\n', u_opt);
fprintf('|dv1|        [m/s]: %.4f\n', norm(u_opt));
fprintf('Costo totale [m/s]: %.4f\n', J_min);

% ── 5. Propagazione arco di trasferimento ────────────────────────────────
x0_D_plus           = x0_D;
x0_D_plus([2,4,6])  = x0_D_plus([2,4,6]) + u_opt(:);
coeff_T             = S \ x0_D_plus;

t_vec = linspace(0, t_transfer, 1000);
x_T = zeros(size(t_vec));
y_T = x_T;  z_T = x_T;

for k = 1:length(t_vec)
    tk    = t_vec(k);
    x_T(k) = real(coeff_T(1)*exp(1i*n*tk) + coeff_T(2)*exp(-1i*n*tk) + 2*(coeff_T(4)/n));
    y_T(k) = real(coeff_T(3) - 3*coeff_T(4)*tk ...
               + 2*1i*(coeff_T(1)*exp(1i*n*tk) - coeff_T(2)*exp(-1i*n*tk)));
    z_T(k) = real(coeff_T(5)*exp(1i*n*tk) + coeff_T(6)*exp(-1i*n*tk));
end

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

% ── 6. PLOT ───────────────────────────────────────────────────────────────

% Colori
col_A = [0.00, 0.45, 0.74];   % blu
col_D = [0.85, 0.33, 0.10];   % arancio-rosso
col_T = [0.47, 0.67, 0.19];   % verde

% --- Figura 1: 3D completo ---
figure('Name','D→A: Vista 3D','NumberTitle','off');

plot3(traj.A.x, traj.A.y, traj.A.z, '--', 'Color', col_A, ...
      'LineWidth', 1.5, 'DisplayName', 'Orbita A');
hold on;
plot3(traj.D.x, traj.D.y, traj.D.z, '--', 'Color', col_D, ...
      'LineWidth', 1.5, 'DisplayName', 'Orbita D');
plot3(x_T, y_T, z_T, '-', 'Color', col_T, ...
      'LineWidth', 2.5, 'DisplayName', 'Arco trasferimento');

% Punti notevoli
plot3(traj.D.x(1), traj.D.y(1), traj.D.z(1), 'o', ...
      'Color', col_D, 'MarkerSize',10, 'MarkerFaceColor', col_D, ...
      'DisplayName', 'Partenza D (dv1)');
plot3(x_T(end), y_T(end), z_T(end), 's', ...
      'Color', col_T, 'MarkerSize',10, 'MarkerFaceColor', col_T, ...
      'DisplayName', 'Aggancio su A (dv2)');
plot3(0, 0, 0, 'k*', 'MarkerSize',14, 'LineWidth',2, ...
      'DisplayName', 'Chief (origine)');

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Trasferimento Ottimale D \rightarrow A (3D)');
legend('Location','best'); grid on; axis equal; view(35,25);

% --- Figura 2: Proiezioni affiancate ---
figure('Name','D→A: Proiezioni','NumberTitle','off');

% Dati per il loop sui piani
plans = { {traj.A.x, traj.A.y, traj.D.x, traj.D.y, x_T, y_T, 'x [m]','y [m]','Piano x-y'}, ...
          {traj.A.x, traj.A.z, traj.D.x, traj.D.z, x_T, z_T, 'x [m]','z [m]','Piano x-z'}, ...
          {traj.A.y, traj.A.z, traj.D.y, traj.D.z, y_T, z_T, 'y [m]','z [m]','Piano y-z'} };


for p = 1:3
    pl = plans{p};
    subplot(1,3,p);
    plot(pl{1}, pl{2}, '--', 'Color', col_A, 'LineWidth',1.5); hold on;
    plot(pl{3}, pl{4}, '--', 'Color', col_D, 'LineWidth',1.5);
    plot(pl{5}, pl{6}, '-',  'Color', col_T, 'LineWidth',2.5);
    plot(pl{3}(1), pl{4}(1), 'o','Color',col_D,'MarkerSize',8,'MarkerFaceColor',col_D);
    plot(pl{5}(end),pl{6}(end),'s','Color',col_T,'MarkerSize',8,'MarkerFaceColor',col_T);
    plot(0, 0, 'k*','MarkerSize',12,'LineWidth',2);
    xlabel(pl{7}); ylabel(pl{8}); title(pl{9});
    grid on; axis equal;
end
legend('Orbita A','Orbita D','Trasferimento','Partenza D','Aggancio','Chief', ...
       'Location','best');
sgtitle('Proiezioni – Trasferimento D \rightarrow A');

% --- Figura 3: Posizione vs tempo ---
figure('Name','D→A: Posizione nel tempo','NumberTitle','off');
t_min = t_vec / 60;

comps   = {'x','y','z'};
data_T_all = [x_T; y_T; z_T];

for p = 1:3
    subplot(3,1,p);
    plot(t_min, data_T_all(p,:), '-', 'Color', col_T, 'LineWidth', 2);
    yline(target_pos(p), '--k', 'LineWidth', 1.2, ...
          'Label', sprintf('Target %s = %.1f m', comps{p}, target_pos(p)));
    xlabel('Tempo [min]');
    ylabel(sprintf('%s [m]', comps{p}));
    title(sprintf('Componente %s – Arco di trasferimento', comps{p}));
    grid on;
end
sgtitle('Evoluzione Temporale – Trasferimento D \rightarrow A');

% =========================================================================
% FUNZIONI LOCALI
% =========================================================================

function J = objective_transfer(u, x0_D, S, n, t, target_vel)
    dv1_mag = norm(u);
    x0_plus = x0_D;
    x0_plus([2,4,6]) = x0_plus([2,4,6]) + u(:);
    coeff = S \ x0_plus;

    vx_f = real( 1i*n*coeff(1)*exp(1i*n*t) - 1i*n*coeff(2)*exp(-1i*n*t));
    vy_f = real(-3*coeff(4) - 2*n*coeff(1)*exp(1i*n*t) - 2*n*coeff(2)*exp(-1i*n*t));
    vz_f = real( 1i*n*coeff(5)*exp(1i*n*t) - 1i*n*coeff(6)*exp(-1i*n*t));

    % dv2 = scarto rispetto alla velocità di A all'arrivo
    dv2_mag = norm([vx_f; vy_f; vz_f] - target_vel);
    J = dv1_mag + dv2_mag;
end

function [c, ceq] = constraints_transfer(u, x0_D, S, n, t, target_pos)
    x0_plus = x0_D;
    x0_plus([2,4,6]) = x0_plus([2,4,6]) + u(:);
    coeff = S \ x0_plus;

    xf = real(coeff(1)*exp(1i*n*t) + coeff(2)*exp(-1i*n*t) + 2*(coeff(4)/n));
    yf = real(coeff(3) - 3*coeff(4)*t ...
            + 2*1i*(coeff(1)*exp(1i*n*t) - coeff(2)*exp(-1i*n*t)));
    zf = real(coeff(5)*exp(1i*n*t) + coeff(6)*exp(-1i*n*t));

    ceq = [xf; yf; zf] - target_pos;
    c   = [];
end


% ── 7. PLOT ──────────────────────────────────────────────────────────────────


traj = struct(); % conterrà i campi A e D




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