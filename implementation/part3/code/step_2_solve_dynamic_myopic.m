clc;
clear;

% bring fundamentals
load ('../output/fundamentals.mat', 'fundamentals');

% unpack variables
[micro_code, T, B, d, mu, N_0] = deal(...
    fundamentals.micro_code, ...
    fundamentals.T, ...
    fundamentals.B, ...
    fundamentals.d, ...
    fundamentals.mu, ...
    fundamentals.N_0);

% parameters
theta   = 4;
kappa   = 2;
delta   = 0;
beta    = -0.2;
alpha   = 0.1;
I       = length(N_0);
N_total = sum(N_0);
maxT    = 50; % number of periods

% vector of wages and labor allocations and guess of real income
N_t         = repmat(ones(I,1), [1 maxT]);
N_t_lag     = repmat(N_0, [1 maxT]); % important to start with correct initial values
w_t         = repmat(ones(I,1), [1 maxT]);
y_t         = repmat(ones(I,1), [1 maxT]);

% solve the model sequentially
tic
for tx = 1:(maxT-1)    

    % guess of real income and wages
    y_guess = y_t(:,tx);
    N_guess = N_t(:,tx);
    w_guess = w_t(:,tx);
    N_lag   = N_t_lag(:,tx);

    % start with solution to N
    dif_tol_outer   = 1;          % initial tolerance
    tol_outer       = 1e-4;       % tolerance parameter
    count_outer     = 1;          % count
    max_count_outer = 1000;       % max count
    while dif_tol_outer > tol_outer && count_outer < max_count_outer 
    
        % compute N_t given N_t_lag and y_t (guess)               
        Btilde       = B.*N_guess.^beta;
        lambda_nom   = Btilde.*(mu.*y_guess').^kappa;
        Omega        = sum(lambda_nom,2).^(1/kappa);
        lambda       = lambda_nom./Omega.^kappa;
        N_new        = lambda'*N_lag;
        N_new        = N_new./sum(N_new)*N_total;

        % solve static equilibrium given implied N_t
        [w_guess, y_new, dif_tol]  = solve_static_elbrm(T, d, alpha, theta, N_guess, w_guess);
        
        % compare implied income and population
        adj_y = y_new./y_guess;
        dif_tol_y = max(abs(1-adj_y));
        y_guess = y_guess.*adj_y.^(0.8);

        adj_N = N_new./N_guess;
        dif_tol_N = max(abs(1-adj_N));
        N_guess = N_guess.*adj_N.^(0.8);
        N_guess = N_guess/sum(N_guess)*N_total;

        % finish
        count_outer = count_outer + 1;
        dif_tol_outer = max([dif_tol_y, dif_tol_N, dif_tol]);

    end

    % update the lagged variable
    N_t_lag(:,tx+1)     = N_guess;
    N_t(:,tx)           = N_guess;
    y_t(:,tx)           = y_guess;
    w_t(:,tx)           = w_guess;    

    N_t(:,tx+1)         = N_guess; % use as a guess for next period
    y_t(:,tx+1)         = y_guess; % use as a guess for next period
    w_t(:,tx+1)         = w_guess; % use as a guess for next period

    % display results
    disp(['Tolerance in inner = ', num2str(dif_tol)])   
    disp(['Tolerance in outer = ', num2str(dif_tol_outer)])   
    disp(['Iterations to solve = ', num2str(count_outer)]) 
    disp(['Period = ', num2str(tx)]) 
end
toc

% export steady state and t+3 to make some maps!
N_ss = N_t(:,end);
N_3  = N_t(:,3);
N_changes_ss = (N_t(:,end) - N_t(:,1))./N_t(:,1);
N_changes_3 = (N_t(:,3) - N_t(:,1))./N_t(:,1);
results = [micro_code, N_ss, 100*N_changes_ss, N_3, 100*N_changes_3];
writematrix(results, '../output/results_myopic.csv');