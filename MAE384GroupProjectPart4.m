% Main script
% Parameters for transmission rate
beta_0 = 0.3;
Lambda = 0.5;
omega1 = 2 * pi; % 1 day periodicity
omega2 = 2 * pi * 100 / 365; % ~3.65 days periodicity

% Time settings
h = 0.1; % Time step (days)
t = 0:h:30; % Time vector

% Initial conditions
S0 = 990;
I0 = 10;
R0 = 0;
N = S0 + I0 + R0; % Total population

% Simulate for omega1 (1 day periodicity)
[S1, I1, R1] = runge_kutta_SIR_periodic(beta_0, Lambda, omega1, 0.1, S0, I0, R0, h, t, N);

% Plot the signals for omega1
figure;
hold on;
plot(t, S1, 'b-', 'DisplayName', 'S(t)');
plot(t, I1, 'r-', 'DisplayName', 'I(t)');
plot(t, R1, 'g-', 'DisplayName', 'R(t)');
xlabel('Time (days)');
ylabel('Population');
title('SIR Model with Periodic Transmission Rate (\omega = 2\pi)');
legend;
grid on;
hold off;

% Perform FFT for omega1
I_fft1 = fft(I1);
f1 = (0:length(I_fft1)-1)*(1/(h*length(I_fft1)));

% Plot spectrum for infected cases I1
figure;
plot(f1, abs(I_fft1));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Infected Cases (\omega = 2\pi)');
xlim([0 0.5]);
grid on;

% Simulate for omega2 (~3.65 days periodicity)
[S2, I2, R2] = runge_kutta_SIR_periodic(beta_0, Lambda, omega2, 0.1, S0, I0, R0, h, t, N);

% Plot the signals for omega2
figure;
hold on;
plot(t, S2, 'b-', 'DisplayName', 'S(t)');
plot(t, I2, 'r-', 'DisplayName', 'I(t)');
plot(t, R2, 'g-', 'DisplayName', 'R(t)');
xlabel('Time (days)');
ylabel('Population');
title('SIR Model with Periodic Transmission Rate (\omega = 2\pi \times 100/365)');
legend;
grid on;
hold off;

% Perform FFT for omega2
I_fft2 = fft(I2);
f2 = (0:length(I_fft2)-1)*(1/(h*length(I_fft2)));

% Plot spectrum for infected cases I2
figure;
plot(f2, abs(I_fft2));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Infected Cases (\omega = 2\pi \times 100/365)');
xlim([0 0.5]);
grid on;

% Function to perform Runge-Kutta 4th Order Method with periodic beta
function [S, I, R] = runge_kutta_SIR_periodic(beta_0, Lambda, omega, gamma, S0, I0, R0, h, t, N)
    % Initialize arrays
    S = zeros(1, length(t));
    I = zeros(1, length(t));
    R = zeros(1, length(t));
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;

    for i = 1:(length(t) - 1)
        % Calculate beta(t)
        beta = beta_0 + Lambda * sin(omega * t(i));
        
        % Define ODE functions
        fS = @(S, I) -(beta/N) * S * I;
        fI = @(S, I) (beta/N) * S * I - gamma * I;
        fR = @(I) gamma * I;

        % Runge-Kutta 4th Order Method
        k1_S = fS(S(i), I(i));
        k1_I = fI(S(i), I(i));
        k1_R = fR(I(i));

        k2_S = fS(S(i) + 0.5 * k1_S * h, I(i) + 0.5 * k1_I * h);
        k2_I = fI(S(i) + 0.5 * k1_S * h, I(i) + 0.5 * k1_I * h);
        k2_R = fR(I(i) + 0.5 * k1_I * h);

        k3_S = fS(S(i) + 0.5 * k2_S * h, I(i) + 0.5 * k2_I * h);
        k3_I = fI(S(i) + 0.5 * k2_S * h, I(i) + 0.5 * k2_I * h);
        k3_R = fR(I(i) + 0.5 * k2_I * h);

        k4_S = fS(S(i) + k3_S * h, I(i) + k3_I * h);
        k4_I = fI(S(i) + k3_S * h, I(i) + k3_I * h);
        k4_R = fR(I(i) + k3_I * h);

        S(i+1) = S(i) + (1/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S) * h;
        I(i+1) = I(i) + (1/6) * (k1_I + 2*k2_I + 2*k3_I + k4_I) * h;
        R(i+1) = R(i) + (1/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R) * h;
    end
end
