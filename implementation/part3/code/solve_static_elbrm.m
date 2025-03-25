function [w, y, dif_tol] = solve_static_elbrm(T, d, alpha, theta, N, w)

    % inner loop parameters
    dif_tol   = 1;          % initial tolerance
    tol       = 1e-6;      % tolerance parameter
    count     = 1;          % count
    max_count = 10;         % max count
    adj_exp   = 0.5;        % update parameter exp (helpful with extreme values)
    adj_sum   = 0.5;        % update parameter sum
    while dif_tol > tol && count < max_count
        
        % income
        w       = w./w(1);  % ensure first country is numeraire    
        E       = w.*N;
    
        % expenditure shares
        T_tilde = T.*(N.^alpha);       
        pi_nom  = T_tilde.*(d.*w).^(-theta); % dim 1 i origin - dim 2 is destination
        Phi     = sum(pi_nom,1);       % sum over dim 1
        pi      = pi_nom./Phi;
        P       = Phi'.^(1/theta);
    
        % total sales
        X       = pi.*E';
        Y       = sum(X,2);
    
        % demand for labor
        N_dem     = Y./w;
        dif_adj_w = N_dem./N;
        dif_tol_w = max(abs(1-dif_adj_w));
        w_new     = w.*dif_adj_w.^adj_exp;
        w_new     = w_new./w_new(1); 
        w         = w*adj_sum + w_new*(1-adj_sum); 
    
        % update
        dif_tol = max(dif_tol_w);
        count   = count + 1;
    
    end

    % export real income
    y = w./P;

end

