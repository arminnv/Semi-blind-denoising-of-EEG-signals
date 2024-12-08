% Parameters
close all
figure 
plot(pulstran(t, [0:0.1:1], 'rectpuls', 0.03))
fs = 1000; % Sampling frequency
t = 0:1/fs:1-1/fs; % Time vector
L = 32; % Length of the filter
mu = 0.01; % Learning rate
sigma = 0.1; % Standard deviation'

% Triangular train signal
signal = 0.5 * sin(2 * pi * 20 * t) .* cos(2 * pi * 400 * t + 0.5);

% Sinusoidal noise
refrence = sawtooth(2 * pi * 10 * t, 0.5);
noise = 2 * refrence .* (1+ 0.3 * randn(size(signal))) .* (1+pulstran(t, [0:0.1:1], 'rectpuls', 0.03));
%noise = 0.1 * randn(size(signal));

% Noisy signal
noisy_signal = signal + noise;

% Initialize filter weights
W = zeros(L, 1);

% Initialize variables
n_iterations = length(noisy_signal);
y = zeros(n_iterations, 1); % Filter output
e = zeros(n_iterations, 1); % Error signal

X = zeros(n_iterations, L);


for n = L+1:n_iterations
     % Input vector
    %x = noisy_signal(n-L+1:n)';
    %X(n, :) = x(n-L+1, n);
    X(n, :) = refrence(n-L+1: n);

    % Filter output
    y(n) = W' * X(n, :)';

    % Error signal
    e(n) = noisy_signal(n) - y(n);
    
    grad = calc_gradient(e, X, L, sigma, n);

    % Update filter coefficients
    W = W + mu * grad;
    %W = W +  mu * x * e(n);

    
    %y = W' * X(:, n); % Filter output
    %e(n) = x(n) - y; % Error signal
    %grad = calc_gradient(e, X, L, sigma);
    %W = W + mu * grad; % Update weights
end

% Output filtered signal
%filtered_signal = filter(W, 1, x);

figure;
subplot(5,1,1);
plot(t, signal);
title('Original Signal');

subplot(5,1,2);
plot(t, noisy_signal);
title('Noisy Signal');

subplot(5,1,3);
plot(t, noise);
title('Noise');

subplot(5,1,4);
plot(t, e');
title('Filtered Signal');

subplot(5,1,5);
plot(t, noisy_signal-e');
title('Estimated noise');


% Gradient calculation function
function grad = calc_gradient(e, X, L, sigma, n)
    
    grad = zeros(size(X, 2), 1); % Initialize gradient with correct size
    for i =n-L:n-1
        diff_e = e(n) - e(i); % Ensure diff_e is a scalar
        diff_X = X(n, :) - X(i, :); % Ensure diff_X is a column vector
        G_sigma = exp(-(diff_e^2) / (4*sigma^2)) / (sqrt(4*pi) * sigma);
        grad = grad + G_sigma * diff_e * diff_X'; % Scalar * Vector multiplication
    end
    grad = grad / (2*sigma^2 * L);
end


