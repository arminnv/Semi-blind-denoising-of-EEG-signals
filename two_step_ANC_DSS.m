% Applying ANC-DSS method and feeding the estimated noise components to
% another ANC algorithm

function denoised = two_step_ANC_DSS(data, refrence, mu, filter_order1, filter_order2, n_components)

    gamma = 1000000;

    scale = max(abs(data), [], 2);
    data = data ./ scale;
    
  
    [denoised, R, r] = ANC_DSS(data, refrence, mu, gamma, filter_order1, n_components);
    
    for ch=1:size(data, 1)
        R = R ./ max(R, [], 2);
        denoised(ch, :) = data(ch, :) - anc(data(ch, :), R(ch,:), mu, gamma, filter_order2);
    end
     
    denoised = denoised .* scale;
    data = data .* scale;
    
end