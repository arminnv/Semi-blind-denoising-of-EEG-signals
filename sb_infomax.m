function [noise, denoised] = sb_infomax(X, refrence, alpha, n_components)

    [num_sources, num_samples] = size(X);
    eta = 0.01;
    max_iter = 60000; 
    
    W = randn(num_sources); 
    
    D = zeros(num_sources, num_samples);
    D(1:n_components, :) = refrence;
    D(1:n_components, :) = D(1:n_components, :)/norm(D(1:n_components, :));

    
    for iter = 1:max_iter   
        disp(iter)
        U = W * X;                        % Estimated sources
        f_U = tanh(U);                    % Non-linearity function (tanh)
        
        % Natural gradient update for W
        delta_W = eta * ((eye(num_sources) - (f_U * U')/num_samples) * W - alpha * D * X');

        W = W + delta_W;
        
        if norm(delta_W) < 1e-10
            disp(iter)
            break
        end
    end
    
    % Output the estimated sources
    estimated_sources = W * X;
    noise = estimated_sources(1:n_components, :);
    estimated_sources(1:n_components, :) = 0;
    denoised = pinv(W) * estimated_sources;
end
