clc;
clear;

% THE PURPOSE OF THIS CODE IS TO PRODUCE TFPS AND AMENITIES THAT CAN BE
% USED IN OUR DYNAMIC MODEL. THIS IS NOT AN INVERSION OF A DYNAMIC
% MODEL.

% bring data
N_i_data   = readmatrix('../output/pop.csv')/1000; % in thousands
Y_i_data   = readmatrix('../output/gdp.csv');
micro_code = readmatrix('../output/micro_code.csv');
dist       = readmatrix('../output/mat_dist.csv');
dist       = dist - diag(diag(dist)); % enforce zero distance with itself

% normalize GDP
Y_i_data = Y_i_data./mean(Y_i_data);

% make trade cost
d0  = 1;
d1  = 0.0005;
d   = d0*exp(d1*dist); % not trying to be realistic!

% make migration cost
mu0  = 1;
mu1  = 0.0010;
mu   = mu0*exp(mu1*dist); % not trying to be realistic!

% parameters
N_total = sum(N_i_data);
theta   = 4;
alpha   = 0.10;
beta    = -0.20;
I       = length(N_i_data);

% guesses
w    = ones(I,1);
N    = ones(I,1);
Tbar = ones(I,1);
Bbar = ones(I,1);

% outer loop parameters
dif_tol_outer   = 1;          % initial tolerance
tol_outer       = 1e-6;       % tolerance parameter
count_outer     = 1;          % count
max_count_outer = 1000;       % max count
adj_exp_outer   = 0.2;        % update parameter exp (helpful with extreme values)
adj_sum_outer   = 0.2;        % update parameter sum

% construct productivities and amenities
while dif_tol_outer > tol_outer && count_outer < max_count_outer

    % inner loop parameters
    dif_tol   = 1;          % initial tolerance
    tol       = 1e-6;       % tolerance parameter
    count     = 1;          % count
    max_count = 10;       % max count
    adj_exp   = 0.8;        % update parameter exp (helpful with extreme values)
    adj_sum   = 0.8;        % update parameter sum
    while dif_tol > tol && count < max_count
        
        % construct productivities and amenities
        T = Tbar.*(N.^alpha);
        B = Bbar.*(N.^beta);
    
        % income
        w       = w./w(1);  % ensure first country is numeraire    
        E       = w.*N;
    
        % expenditure shares
        pi_nom  = T.*(d.*w).^(-theta); % dim 1 i origin - dim 2 is destination
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
    
        % utility
        U         = w./P.*B;
        U_mean    = mean(U);
        dif_adj_U = U./U_mean;
        dif_tol_U = max(abs(1-dif_adj_U));
        N_new     = N.*dif_adj_U.^adj_exp;
        N         = N*adj_sum + N_new*(1-adj_sum);
        N         = (N./sum(N)).*N_total;
    
        % update
        dif_tol = max([dif_tol_w,dif_tol_U]);
        count   = count + 1;
    
    end

    % display results
    disp(['Tolerance in wage = ', num2str(dif_tol_w)])
    disp(['Tolerance in labor = ', num2str(dif_tol_U)])
    disp(['Iteration = ', num2str(count)]) 

    % Adjust Tbar
    Y_norm       = Y./mean(Y);
    dif_adj_Tbar = Y_i_data./Y_norm;
    dif_tol_Tbar = max(abs(1-dif_adj_Tbar));
    Tbar_new     = Tbar.*dif_adj_Tbar.^adj_exp_outer;
    Tbar         = Tbar*adj_sum_outer + Tbar_new*(1-adj_sum_outer);
    Tbar         = Tbar./mean(Tbar);

    % Adjust Bbar
    dif_adj_Bbar = N_i_data./N;
    dif_tol_Bbar = max(abs(1-dif_adj_Bbar));
    Bbar_new     = Bbar.*dif_adj_Bbar.^adj_exp_outer;
    Bbar         = Bbar*adj_sum_outer + Bbar_new*(1-adj_sum_outer);
    Bbar         = Bbar./mean(Bbar);

    % update
    dif_tol_outer = max([dif_tol_Tbar,dif_tol_Bbar]);
    count_outer   = count_outer + 1;

    % display results
    disp(['Tolerance in Tbar = ', num2str(dif_tol_Tbar)])
    disp(['Tolerance in Bbar = ', num2str(dif_tol_Bbar)])
    disp(['Iteration = ', num2str(count_outer)]) 

end

% export data to make some maps!
fundamentals = struct(...
    'micro_code', micro_code, ...
    'T', T, ...
    'B', B, ...
    'd', d, ...
    'mu', mu, ...
    'N_0', N);
save('../output/fundamentals.mat', 'fundamentals');
