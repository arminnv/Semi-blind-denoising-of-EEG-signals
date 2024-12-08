close all;

% Parameters
fs = 1000; % Sampling frequency
t = 0:1/fs:1-1/fs; % Time vector

% Triangular train signal
n_channels = 10;
signal = zeros(n_channels, length(t));
refrence = zeros(n_channels, length(t));
noise = zeros(n_channels, length(t));

for ch=1:n_channels
    signal(ch, :) = ch^2 * 0.5 * sin(2 * pi * 20 * t) .* cos(2 * pi * 400 * t + 0.5);
    refrence(ch, :) = ch^3 * sawtooth(2 * pi * 10 * t, 0.5);
    noise(ch, :) = 2 * refrence(ch, :) .* (1+ 0.3 * randn(1, size(signal, 2))) .* (1+pulstran(t, [0:0.1:1], 'rectpuls', 0.03));
end

% Sinusoidal noise
%noise = 0.2 * sin(2 * pi * 200 * t);
%noise = 0.1 * randn(size(signal));

% Noisy signal
noisy_signal = signal + noise;

mu = 0.01; % Step size 
filter_order = 64; % Number of filter coefficients
X_denoised = ANC_DSS(noisy_signal, refrence, mu, filter_order);

RRMSE = sqrt(sumsqr(X_denoised - signal))/sqrt(sumsqr(signal));
disp(RRMSE)

channel = 1;
% Plot results
figure;
subplot(5, 1, 1);
plot(t, signal(channel, :));
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 2);
plot(t, noisy_signal(channel, :));
title('Noisy Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 3);
plot(t, X_denoised(channel, :));
title('Filtered Signal using LMS');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 4);
plot(t, noisy_signal(channel, :)-X_denoised(channel, :));
title('Estimated noise');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 5);
plot(t, noise(channel, :));
title('Estimated noise');
xlabel('Time (s)');
ylabel('Amplitude');


scale = var(noisy_signal, [], 2).^0.5;
noisy_signal = noisy_signal ./ scale;
refrence = refrence ./ var(refrence, [], 2).^0.5;

X_denoised = zeros(size(noisy_signal));
for ch=1:n_channels
    X_denoised(ch, :) = noisy_signal(ch, :) - anc(noisy_signal(ch, :), refrence(ch, :), mu, filter_order);
end

X_denoised = X_denoised .* scale;
noisy_signal = noisy_signal .* scale;

RRMSE = sqrt(sumsqr(X_denoised - signal))/sqrt(sumsqr(signal));
disp(RRMSE)

channel = 1;
% Plot results
figure;
subplot(5, 1, 1);
plot(t, signal(channel, :));
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 2);
plot(t, noisy_signal(channel, :));
title('Noisy Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 3);
plot(t, X_denoised(channel, :));
title('Filtered Signal using LMS');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 4);
plot(t, noisy_signal(channel, :)-X_denoised(channel, :));
title('Estimated noise');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5, 1, 5);
plot(t, noise(channel, :));
title('Estimated noise');
xlabel('Time (s)');
ylabel('Amplitude');
