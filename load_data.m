% Loading data using EEGLAB and applying PCA

fs = EEG.srate;
raw_data = EEG.data;
data = raw_data;

% Perform PCA
[coeff, score, ~, ~, explained] = pca(data');

% Use only the first principal component to reconstruct the data
first_pc = coeff(:, 1);
first_component = (score(:, 1) * first_pc')';

% Remove first PC
data_filt = data - first_component;

%%

EEG.data = data_filt;
% Pefrom fMRIb in EEGLAB to remove Gradient Artifact

%%

data_noga = EEG.data;
ga = data_filt - data_noga;

%%

% Plotting the signals

types = {EEG.event.type}; % Events
index = find(strcmp(types, '105')); % Time index of stimulation 105
stim_105 = types(index:end);

t_start = EEG.event(index).latency; 
t_end = t_start + 10*fs; 
channel = 20; % Signal to plot

figure
subplot(5, 1, 1)
plot(raw_data(channel, t_start:t_end))
title("Raw data")
subplot(5, 1, 2)
plot(first_component(channel, t_start:t_end))
title("First Principal Component")
subplot(5, 1, 3)
plot(data_filt(channel, t_start:t_end))
title("First PC Removed")
subplot(5, 1, 4)
plot(ga(channel, t_start:t_end))
title("Extracted GA")
subplot(5, 1, 5)
plot(data_noga(channel, t_start:t_end))
title("GA and First PC removed")

%%
% Finding order of stimulations

types = {EEG.event.type};
stimulation_list = [];

for i=1:length(types)
    if strcmp(types{i}(1:2), '10')
        if str2double(types{i}(3)) > 1
            stimulation_list = [stimulation_list, types(i)];
        end
    end
end

disp(stimulation_list)

%%
% Extracting trials

index_start = find(strcmp(types, '105'));
index_end = find(strcmp(types, '103'));
trial_starts = [];
names = {};
flags = {'11', '12', '13', '21', '31', '41', '42', '51'};

for i=index_start:index_end
    for flag=flags
        if strcmp(types{i}, flag)
            trial_starts = [trial_starts, EEG.event(i).latency];
            names{end+1} = flag;
        end
    end
end

%%
% Initializing hyper-parameters
mu = 0.05; % Step size
% mu = 0.02
n_components = 5; % Number of noise components to estimate
filter_order1 = 40; % Must be an even Integer
filter_order2 = 40; % Must be an even Integer
trial_length = 5*fs;
initial_length = 1*fs;

exp_num = 0;

% Denoising all trials using ANC_DSS
for trial=1:length(trial_starts)
    t_start = trial_starts(trial)-initial_length;
    t_end = t_start + trial_length - 1;
    channel = 20;
    if strcmp(names{trial},'11') || strcmp(names{trial},'12') || strcmp(names{trial},'13')
        exp_num = exp_num + 1;
    end
    
    refrence = [first_component(channel, t_start: t_end); ga(channel, t_start: t_end)];
    trial_signal = data_noga(:, t_start: t_end);
    denoised = two_step_ANC_DSS(trial_signal, refrence, mu, filter_order1, filter_order2, n_components); 
    denoised = denoised(:, initial_length:end);
    trial_signal = trial_signal(:, initial_length:end);

    name = "Exp" + num2str(exp_num)  + "-flag" + names{trial}; 
    path = "Results\ANC_DSS\";
    plot_results(trial_signal, denoised, channel, name, path);
    save(path + name + ".mat", 'denoised');

end 

%%

alpha = 0.5;

exp_num = 0;
for trial=1:1
    t_start = trial_starts(trial)-initial_length;
    t_end = t_start + trial_length - 1;
    channel = 20;
    if strcmp(names{trial},'11') || strcmp(names{trial},'12') || strcmp(names{trial},'13')
        exp_num = exp_num + 1;
    end

    trial_signal = data_noga(:, t_start: t_end);

    refrence = [first_component(channel, t_start: t_end); ga(channel, t_start: t_end)];
    trial_signal = data_noga(:, t_start: t_end);
    denoised = two_step_ANC_DSS(trial_signal, refrence, mu, filter_order1, filter_order2, n_components); 
    
    scale = var(trial_signal, [], 2).^0.5;
    refrence = refrence./var(refrence, [], 2).^0.5;

    [~, denoised] = sb_infomax(denoised./scale, refrence(1,:), alpha, 1);
    
    denoised = denoised .* scale;

    denoised = denoised(:, initial_length:end);
    trial_signal = trial_signal(:, initial_length:end);

    name = "Exp" + num2str(exp_num)  + "-flag" + names{trial}; 
    path = "Results\sb_ICA\";
    plot_results(trial_signal, denoised, channel, name, path);
    %save(path + name + ".mat", 'denoised');

end 



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





%%

function plot_results(data, denoised, channel, name, path)
    figure 
    
    subplot(3, 1, 1)
    plot(data(channel, :))
    title('Noisy') 
    subplot(3, 1, 2)
    plot(denoised(channel, :))
    title('Denoised')
    subplot(3, 1, 3)
    plot(data(channel, :)-denoised(channel, :))
    title('Noise')
    sgtitle(name)

    exportgraphics(gcf, path + "Plots\" + name + ".png", 'Resolution', 300);
end
