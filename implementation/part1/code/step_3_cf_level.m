clc;
clear;

% load data
load('../output/calibration.mat', 'calibration');
d     = calibration.d;
T     = calibration.T;
N     = calibration.N;
theta = calibration.theta;
I     = calibration.I;

%% SIMULATE BASELINE
w               = ones(I,1);  % guess for w
adj_exp_N       = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum_N       = 0.5;        % update parameter sum

% we will run this for two Ts. We will increase T(1) by 20 percent for the
% second round
T_cases = [T, T];
T_cases(1,2) = T_cases(1,2)*1.2;
W_cases = zeros(size(T_cases));
w_cases = zeros(size(T_cases));

% loop related parameters
for tx = 1:2

    % pick T    
    T = T_cases(:,tx);
    
    % initiate loop
    dif_tol_N   = 1;          % initial tolerance
    tol_N       = 1e-6;       % tolerance parameter
    count_N     = 1;          % count
    max_count_N = 1000;       % max count
    
    % solve the model in the inner loop
    while dif_tol_N > tol_N && count_N < max_count_N
    
        % income
        w       = w./w(1);  % ensure first country is numeraire    
        E       = w.*N;
    
        % expenditure shares
        pi_nom  = T.*(d.*w).^(-theta); % dim 1 i origin - dim 2 is destination
        Phi     = sum(pi_nom,1);       % sum over dim 1
        pi      = pi_nom./Phi;
    
        % total sales
        X       = pi.*E';
    
        % demand for labor
        N_dem     = sum(X,2)./w;
        dif_adj_N = N_dem./N;
        dif_tol_N = max(abs(1-dif_adj_N));
        w_new     = w.*dif_adj_N.^adj_exp_N;
        w_new     = w_new./w_new(1); 
        w         = w*adj_sum_N + w_new*(1-adj_sum_N); 
    
        % update
        count_N   = count_N + 1;
    
    end

    P_index = Phi.^(-1/theta);
    w_cases(:,tx) = w;
    W_cases(:,tx) = w./P_index';
end

% see relative wages
w_hat = w_cases(:,2)./w_cases(1);

% make some figure
D_rincome = 100*(W_cases(:,2)-W_cases(:,1))./W_cases(:,1);
D_rincome(1) = max(D_rincome(2:end)); % FOR VISUALIZATION WE WILL REMOVE THE GAINS FOR 1
I_side  = sqrt(I);
D_rincome_mat = reshape(D_rincome, [I_side, I_side]);
imagesc(D_rincome_mat);
colorbar
saveas(gcf, '../output/map_D_real_income_levels.png');

