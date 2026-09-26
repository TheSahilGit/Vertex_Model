%% ==========================================================
%  ODE system with ΔA(t) being:
%     (1) random in time   OR
%     (2) sinusoidal in time
%  Choose mode using the switch "mode"
% ==========================================================

clear; clc;

%% -------------------------------
% Parameters
% -------------------------------
A_rho = 0.015;
n     = 10;
K     = 0.10;
A_R   = 0.015;
A_M   = 0.015;
D_rho = 0.015;
D_R   = 0.004;
D_M   = 0.015;

h = @(x) double(x >= 0);

%% -------------------------------
% Time settings
% -------------------------------
t0 = 0; 
tf = 500;

%% ==========================================================
% Choose mode:
% mode = "random";      % noisy ΔA(t)
 mode = "random";         % sinusoidal ΔA(t)
% mode = "random";        % <--- change this
% ==========================================================

%% -------------------------------
% Define ΔA(t)
% -------------------------------
switch mode
    case "random"
        % Random time series for ΔA(t)
        Npts = 2001;
        t_rand = linspace(t0, tf, Npts);

        % Gaussian random noise
        DeltaA_rand = randn(size(t_rand)) ;

        % Smooth it for biological realism
        DeltaA_rand = smoothdata(DeltaA_rand, 'gaussian', 40);

        % Interpolant
        DeltaA_fun = @(t) interp1(t_rand, DeltaA_rand, t, 'linear', 'extrap');

    case "sin"
        % Sinusoidal forcing: amplitude & frequency
        A0 = 0.0;     % baseline
        A1 = 0.2;     % oscillation amplitude
        omega = 2*pi/50;  % period = 50

        DeltaA_fun = @(t) A0 + A1 * sin(omega * t);

    otherwise
        error("Unknown mode")
end

%% -------------------------------
% ODE system
% -------------------------------
odefun = @(t, y) [
    A_rho * h(DeltaA_fun(t)) * (DeltaA_fun(t)^n) / (K^n + DeltaA_fun(t)^n) * (1 - y(1)) ...
         - D_rho * y(1);

    A_R * y(1) * (1 - y(2)) - D_R * y(2);

    A_M * y(2) * (1 - y(3)) - D_M * y(3)
];

%% Initial conditions
y0 = [0.0; 0.0; 0.0];

%% Solve
[t, y] = ode45(odefun, [t0 tf], y0);

rho = y(:,1);
R   = y(:,2);
M   = y(:,3);

%% -------------------------------
% Plot ΔA(t)
% -------------------------------
figure;

subplot(2,1,1)
tt = linspace(t0, tf, 2000);
plot(tt, DeltaA_fun(tt), 'k', 'LineWidth', 2)
ylabel('\DeltaA(t)')
title(['Mode: ' mode])
grid on

subplot(2,1,2)
hold on
plot(t, rho, 'LineWidth', 2)
plot(t, R,   'LineWidth', 2)
plot(t, M,   'LineWidth', 2)
xlabel('Time')
ylabel('Values')
legend('\rho','ROCK','Myosin')
grid on
set(gca,'FontSize',16);



figure; 

plot(M, rho, 'LineWidth', 2);     % trajectory
hold on;
plot(M(1), rho(1), 'go', 'MarkerSize', 10, 'LineWidth', 2);  % start point
plot(M(end), rho(end), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % end point

xlabel('M(t)');
ylabel('\rho (t)');

title('Phase space:  M vs  \rho');
grid on;
set(gca,'FontSize',16);
axis tight;




