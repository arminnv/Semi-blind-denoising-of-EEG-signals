function [X_den, X_noise, R]  = ANC_DSS(X_noisy, refrence, mu, gamma, filter_order, M)
    refrence = refrence ./ var(refrence, [], 2).^0.5;
    disp(size(refrence))
    n_channels = size(X_noisy,1);
    L = size(X_noisy, 2);

    [U, S, V] = svd(X_noisy*X_noisy');
    Z = S^(-1/2) * U' * X_noisy;
    disp(size(Z))
    
    %M = 2;
    W = zeros(36, M);
    R = zeros(M, L);
    disp(L)

    w = randn(n_channels,1);
    for update=1:1
        for p=1:M
            for iteration=1:100
                r = w'*Z;
                r_new = anc(r, refrence, mu, gamma, filter_order);
                %refrence = r_new;
        
                w_new = Z*r_new';
                w_new = (eye(n_channels)-W(:, 1:p-1)*W(:, 1:p-1)')*w_new;
                w = w_new/norm(w_new);
            end
            R(p, :) = r_new;
            W(:, p) = w_new;
    
            Z = Z - W(:, p) * R(p, :);
            %disp(w_new)
        end
        %refrence = R;
        %refrence = refrence ./ var(refrence, [], 2).^0.5;
    end
    disp(size(R))
    X_noise = U * S^(1/2) * W * R;
    X_den =  X_noisy - X_noise;
end