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
delta   = 0.95^30; % assume 30 years lag between discrete periods
beta    = -0.2;
alpha   = 0.1;
I       = length(N_0);
N_total = sum(N_0);
maxT    = 50; % number of periods

%% SOLVE FOR THE STEADY STATE
N_guess     = ones(I,1)/sum(ones(I,1))*N_total;
y_guess     = ones(I,1);
N_lag       = N_guess;
w_guess     = ones(I,1);
Omega_guess = ones(I,1);

% start with solution to N
dif_tol_outer   = 1;          % initial tolerance
tol_outer       = 1e-4;       % tolerance parameter
count_outer     = 1;          % count
max_count_outer = 1000;       % max count
while dif_tol_outer > tol_outer && count_outer < max_count_outer 

    % compute N_t given N_t_lag and y_t (guess)               
    Btilde       = B.*N_guess.^beta;
    lambda_nom   = Btilde.*(mu.*y_guess'.*Omega_guess'.^delta).^kappa;
    Omega_new    = sum(lambda_nom,2).^(1/kappa);
    lambda       = lambda_nom./Omega_new.^(kappa*(1-delta));
    N_new        = lambda'*N_lag;
    N_new        = N_new./sum(N_new)*N_total;

    % solve static equilibrium given implied N_t
    [w_guess, y_new, dif_tol]  = solve_static_elbrm(T, d, alpha, theta, N_guess, w_guess);
    
    % compare implied to guess
    adj_y = y_new./y_guess;
    dif_tol_y = max(abs(1-adj_y));
    y_guess = y_guess.*adj_y.^(0.8);

    adj_N = N_new./N_guess;
    dif_tol_N = max(abs(1-adj_N));
    N_guess = N_guess.*adj_N.^(0.8);
    N_lag = N_guess;

    adj_Omega = Omega_new./Omega_guess;
    dif_tol_Omega = max(abs(1-adj_Omega));
    Omega_guess = Omega_guess.*adj_Omega.^(0.8);    

    % finish
    count_outer = count_outer + 1;
    dif_tol_outer = max([dif_tol_y, dif_tol_N, dif_tol_Omega, dif_tol]);
    
end

% display results
disp(['Solved SS - Tolerance in inner = ', num2str(dif_tol)])   
disp(['Solved SS - Tolerance in outer = ', num2str(dif_tol_outer)])   
disp(['Solved SS - Iterations to solve = ', num2str(count_outer)]) 

% SS results
Omega_ss   = Omega_guess;
N_ss       = N_guess;
y_ss       = y_guess;
w_ss       = w_guess;
B_tilde_ss = Btilde;
lambda_ss  = lambda;

%% SOLVE THE DYNAMIC EQUILIBRIUM
y_t      = repmat(y_ss, [1 maxT]);
w_t      = repmat(w_ss, [1 maxT]);
y_t_new  = repmat(y_ss, [1 maxT]);
w_t_new  = repmat(w_ss, [1 maxT]);
N_t      = repmat(N_ss, [1 maxT]);
N_t_lag  = repmat(N_0, [1 maxT]);
Omega_t  = repmat(Omega_ss, [1 maxT]);
Btilde_t = repmat(B_tilde_ss, [1 maxT]);
lambta_t = repmat(lambda_ss, [1 1 maxT]);

% start with solution to N
tic;
dif_tol_inner   = zeros(maxT,1);          % to ensure inner loop solution
dif_tol_outer   = 1;          % initial tolerance
tol_outer       = 1e-4;       % tolerance parameter
count_outer     = 1;          % count
max_count_outer = 1000;       % max count
while dif_tol_outer > tol_outer && count_outer < max_count_outer 

    
    % backward induction and migration choices
    for tx = maxT-1:-1:1     
        lambda_nom       = Btilde_t(:,tx).*(mu.*y_t(:,tx)'.*Omega_t(:,tx+1)'.^delta).^kappa;
        Omega_t(:,tx)    = sum(lambda_nom,2).^(1/kappa);
        lambda_t(:,:,tx) = lambda_nom./Omega_t(:,tx).^kappa;
    end

    % forward induction
    for tx = 1:(maxT-1)                
        N_t(:,tx) = lambda_t(:,:,tx)'*N_t_lag(:,tx);
        N_t(:,tx) = N_t(:,tx)./sum(N_t(:,tx))*N_total;
        N_t_lag(:,tx+1) = N_t(:,tx);
    end

    % solve market clearing for every period given sequence of N_t
    for tx = 1:(maxT-1)
        w_guess = w_t(:,tx);
        N_guess = N_t(:,tx);
        [w_guess, y_guess, dif_tol]  = solve_static_elbrm(T, d, alpha, theta, N_guess, w_guess);

        w_t_new(:,tx) = w_guess;
        y_t_new(:,tx) = y_guess;
        dif_tol_inner(tx) = dif_tol;
    end

    % adjust guesses
    adj_y = y_t_new./y_t;
    dif_tol_y = max(max(abs(1-adj_y)));
    y_t = y_t.*adj_y.^(0.8);
    
    adj_w = w_t_new./w_t;
    dif_tol_w = max(max(abs(1-adj_w)));
    w_t = w_t.*adj_w.^(0.8);

    % finish
    count_outer = count_outer + 1;
    dif_tol_outer = max([dif_tol_y, dif_tol_w, max(dif_tol_inner)]);

    % display results
    disp(['Tolerance in inner = ', num2str(dif_tol)])   
    disp(['Tolerance in outer = ', num2str(dif_tol_outer)])   
    disp(['Iterations to solve = ', num2str(count_outer)])     
end
toc

% export steady state and t+3 to make some maps!
N_ss = N_t(:,end);
N_3  = N_t(:,3);
N_changes_ss = (N_t(:,end) - N_t(:,1))./N_t(:,1);
N_changes_3 = (N_t(:,3) - N_t(:,1))./N_t(:,1);
results = [micro_code, N_ss, 100*N_changes_ss, N_3, 100*N_changes_3];
writematrix(results, '../output/results_forward_looking.csv');
