subplot(4, 1, 4)
Xf = abs(fft(data(channel, t_start:t_end)));
L = size(Xf, 2);
f = (0:L/2-1)/L*fs;
plot(f, Xf(1:L/2))

%%

fs = EEG.srate;
raw_data = EEG.data;

data = raw_data;

% Perform PCA
[coeff, score, ~, ~, explained] = pca(data');

% Use only the first principal component to reconstruct the data
first_pc = coeff(:, 1);
first_component = (score(:, 1) * first_pc')';

data_filt = data - first_component;
channel = 20;

%%
figure
subplot(6, 1, 1)
plot(EEG_Sig(channel, t_start:t_end))
title("EEG")
subplot(6, 1, 2)
plot(GA_Sig(channel, t_start:t_end))
title("GA")
subplot(6, 1, 3)
plot(GVS_Sig(channel, t_start:t_end))
title("GVS")
subplot(6, 1, 4)
plot(raw_data(channel, t_start:t_end))
title("Sim")
subplot(6, 1, 5)
plot(data_filt(channel, t_start:t_end))
title("PCA")
subplot(6, 1, 6)
plot(first_component(channel, t_start:t_end))
title("refrence")

%%
EEG.data = data_filt;
% Pefrom fMRIb
%%
data_noga = EEG.data;
ga = data_filt - data_noga;

%%
scale = var(data, [], 2).^0.5;
data = data ./ scale;


denoised = zeros(size(data));
for ch=1:size(data, 1)
    refrence = [first_component(ch, t_start: t_end)];
    refrence = refrence ./ var(refrence, [], 2).^0.5;
    denoised(ch, :) = data(ch, :) - anc(data(ch, :), refrence, mu, filter_order);
end

denoised = denoised .* scale;
data = data .* scale;


%%
types = {EEG.event.type};
index = find(strcmp(types, '101')); %105

t_start = EEG.event(index).latency;
t_end = t_start + 5*fs;

channel = 10;
data = data_noga(:, t_start: t_end);

%eeg = double(EEG_Sig(:, t_start:t_end));

%a = 1;
%b = 1;
%ngft = 30;

%A = gsp_learn_graph_log_degrees(squareform(pdist(double(eeg), "correlation")), a, b);
%[X_hat, X_LPF] = gft_coef(data, A, ngft);


% LMS Parameters
mu = 0.02; % Step size
gamma = 1000000;

filter_order = 30; % 30  36
refrence = [first_component(channel, t_start: t_end); ga(channel, t_start: t_end)];
[denoised, R, r] = ANC_DSS(data, refrence, mu, gamma, filter_order, 5); %5

r = r ./ var(r, [], 2).^0.5; 

scale = var(data, [], 2).^0.5;
data = data ./ scale;
denoised = denoised ./ scale;

for ch=1:size(data, 1)
    R = R ./ var(R, [], 2).^0.5;
    %denoised(ch, :) = data(ch, :) - anc(data(ch, :), R(ch,:), mu, gamma, filter_order);
end


denoised = denoised .* scale;
data = data .* scale;

%%

r = r ./ var(r, [], 2).^0.5;

scale = var(data, [], 2).^0.5;
data = data ./ scale;

denoised = data;
for i=1:size(r, 5)
    denoised = sb_infomax(denoised, r(i,:), 1.5, 1);
end

denoised = denoised .* scale;
data = data .* scale;
%%
noise = data - denoised;

nfft = 10000;
window = 5000;
n_ovelap = 2000;

shame_clean = EEG_Sig;

[p_denoised, f] = pwelch(denoised(channel, :), [], [], [], fs);
[p_noisy, f] = pwelch(data(channel, :), [], [], [], fs);
[p_sham, f_sham] = pwelch(shame_clean(channel, :), [], [], [], fs);
%[p_original, f] = pwelch(eeg(channel, :), [], [], [], fs);

ind = f<70;

figure 
subplot(5, 1, 1)
hold on
plot(data(channel, :))
%plot(eeg(channel, :))
hold off
%legend('noisy', 'denoised', 'original');
title('Noisy')

subplot(5, 1, 2)
hold on
plot(denoised(channel, :))
%plot(eeg(channel, :))
hold off
title('Denoised')

subplot(5, 1, 3)
hold on
%plot(f(ind), p_denoised(ind))
ind_sham = f_sham<70;
plot(f_sham(ind_sham), p_sham(ind_sham))
%plot(f(ind), p_noisy(ind))
%plot(f(ind), p_original(ind))
hold off
legend('denoised', 'sham');
title('PSD')
sgtitle('Before ICA')

subplot(5, 1, 4)
plot(EEG_Sig(channel, t_start: t_end))

subplot(5, 1, 5)
plot(noise(channel, :))

%subplot(5, 1, 4)
%plot(f(ind), p_original(ind))
%plot(abs(fft(eeg(channel,:), 2000)))
%plot(X_LPF(channel, :))
%title("Extracted Noise")

%subplot(5, 1, 5)
%plot(f(ind), p_original(ind))
%plot(abs(fft(eeg(channel,:), 2000)))
%plot(data(channel, :)-eeg(channel, :))
%title("Real Noise")
%plot(data(channel, :)-denoised(channel, :))

%subplot(4, 1, 4)

%plot(refrence(2, :))
%title("GA noise")

%figure

%plot(f, Xf(1:L/2))

%%
EEG.data = denoised;

%%
IC_removed = EEG.data;
ic = denoised - IC_removed;

noise = data - IC_removed;

nfft = 10000;
window = 5000;
n_ovelap = 2000;

[p_denoised, f] = pwelch(IC_removed(channel, :), [], [], [], fs);
[p_noisy, f] = pwelch(data(channel, :), [], [], [], fs);
[p_sham, f_sham] = pwelch(shame_clean(channel, :), [], [], [], fs);
%[p_original, f] = pwelch(eeg(channel, :), [], [], [], fs);

ind = f<70;

figure 
subplot(3, 1, 1)
hold on
plot(data(channel, :))
%plot(eeg(channel, :))
hold off
title('Noisy')
%legend('noisy', 'denoised', 'original');

subplot(3, 1, 2)
hold on
plot(IC_removed(channel, :))
%plot(eeg(channel, :))
title('Denoised')
hold off

subplot(3, 1, 3)
hold on
plot(f(ind), p_denoised(ind))
ind_sham = f_sham<70;
plot(f_sham(ind_sham), p_sham(ind_sham))
%plot(f(ind), p_noisy(ind))
%plot(f(ind), p_original(ind))
hold off
title('PSD')
%legend('denoised', 'noisy', 'original');
sgtitle("After ICA")
%%
EEG.data = noise;


%%
[p_ica, f] = pwelch(IC_removed(channel, :), [], [], [], fs);
[p_noisy, f] = pwelch(data(channel, :), [], [], [], fs);
[p_sham, f_sham] = pwelch(shame_clean(channel, :), [], [], [], fs);
%[p_original, f] = pwelch(eeg(channel, :), [], [], [], fs);

ind = f<70;

figure 
subplot(3, 1, 1)
hold on
plot(data(channel, :))
%plot(eeg(channel, :))
hold off
%legend('noisy', 'denoised', 'original');

subplot(3, 1, 2)
hold on
plot(IC_removed(channel, :))
%plot(eeg(channel, :))
hold off

subplot(3, 1, 3)
hold on
plot(f(ind), p_ica(ind))
ind_sham = f_sham<70;
plot(f_sham(ind_sham), p_sham(ind_sham))
%plot(abs(fft(shame_clean(channel, :))))
%plot(f(ind), p_noisy(ind))
%plot(f(ind), p_original(ind))
hold off
%legend('denoised', 'noisy', 'original');


%%
channel = 20;
shame_clean = EEG.data;

figure 
subplot(1, 1, 1)
hold on
plot(sham_clean(channel, :))
hold off
%legend('noisy', 'denoised', 'original');




function [X_hat, X_LPF] = gft_coef(X, A, ngft)
    D = diag (sum (A, 1)); % degree matrix
    L = D - A; % laplacian matrix
    
    [V,D] = eig(L);
    [D,I] = sort(diag(D), 'ascend');
    V = V(:, I(1:ngft));
    
    X_hat = V' * X;
    X_LPF = V * X_hat;
    %X_hat = normalize(X_hat, 1, "zscore");
end

%%


