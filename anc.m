function noise  = anc(noisy_signal, refrence, mu, gamma, filter_order)
    
    n_iterations = length(noisy_signal); % Number of iterations
    
    % Initialize variables
    h = zeros(filter_order, size(refrence, 1)); % Adaptive filter coefficients
    y = zeros(n_iterations, 1); % Filter output
    e = zeros(n_iterations, 1); % Error signal
    
    % Adaptive filtering
    for n = filter_order/2:n_iterations-filter_order/2
        
        % Input vector
        x = refrence(:, n-filter_order/2+1:n+filter_order/2)';
    
        % Filter output
        y(n) = sum(h .* x, "all");
    
        % Error signal
        e(n) = noisy_signal(n) - y(n);
    
        % Update filter coefficients
        h = h + mu * x * e(n) ./ (1); %+ sum(x.^2, 1)/gamma);
    end
    e = e';
    noise = noisy_signal - e;

end