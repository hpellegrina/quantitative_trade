clc;
clear;

% load data
load('../output/data_square_world.mat');
E_data    = data.E;
X_data    = data.X;
N         = data.N;

% Comment: Notice that this is how the data will come to you. You observe
% total expenditure (you might need to build that value, because we
% often do not observe sales of a country to itself). And you also observe
% trade. 
% Let us first recover trade cost first using gravity equations. I always 
% prefer to do that in stata, because we might have missing
% values, or because you might want to use PPMLHDE. But to keep things more
% parsimonious, for this course, I'll run the gravity here in matlab.

%% RECOVER TRADE COSTS
% recover residuals.
I         = length(E_data);                     % recover number of countries
FE_o      = kron(eye(I),ones(I,1));             % make origin FE
FE_d      = kron(ones(I,1),eye(I));             % make destination FE
y_var     = log(X_data(:));                     % DV    
X_var     = [FE_o, FE_d];                       % EV
beta      = pinv(X_var'*X_var)*(X_var'*y_var);  % plim to avoid multicollinearity
res       = y_var - X_var*beta;

% with residuals, we can use an assumption on trade elasticity
theta     = 4;
tcost     = exp(res).^(-(1/theta));
tcost_mat = reshape(tcost, [I, I]);
d         = tcost_mat./diag(tcost_mat);

%% RECOVER TFPS

% prepare targets in the data
E_target = E_data/E_data(1);

% guesses and tolerances
T               = ones(I,1);  % guess for T
w               = ones(I,1);  % guess for w
adj_exp_N       = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum_N       = 0.5;        % update parameter sum
adj_exp_T       = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum_T       = 0.5;        % update parameter sum

% loop related parameters
dif_tol_T       = 1;          % initial tolerance
tol_T           = 1e-6;       % tolerance parameter
count_T         = 1;          % count
max_count_T     = 1000;       % max count

% search for values of T in the outer loop
while dif_tol_T > tol_T && count_T < max_count_T

    % loop related parameters
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
        Pi      = sum(pi_nom,1);       % sum over dim 1
        pi      = pi_nom./Pi;
    
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

    % compare model implied with actual gross output
    E_model   = sum(X,2);
    E_model   = E_model./E_model(1);
    dif_adj_T = E_target./E_model;
    dif_tol_T = max(abs(1-dif_adj_T));
    T_new     = T.*dif_adj_T.^(adj_exp_T);
    T_new     = T_new./T_new(1);    
    T         = T*adj_sum_T + T_new*(1-adj_sum_T);
    T         = T/T(1); % ensure normalization

    % display results
    disp(['Tolerance in N (inner) = ', num2str(dif_tol_N)]); 
    disp(['Tolerance in T (outer) = ', num2str(dif_tol_T)]); 
    disp(['count (outer) = ', num2str(count_T)]);
end

% export calibration
calibration = struct(...
        'd', d, ...
        'T', T, ...
        'N', N, ...
        'theta', theta, ...
        'I', I);
save('../output/calibration.mat', 'calibration');

% notice that trade flows are perfectly rationalized
X_check = X./X_data;

% let's check results if they are the same as the fundamentals
load('../output/fundamental_for_simulation.mat', 'fundamentals');
T_fund  = fundamentals.T;
d_fund  = fundamentals.d;
T_check = T./T_fund;
d_check = d./d_fund;

% We perfectly rationalize the data, but fundamentals are not exactly the
% same. that's because we are under-identified. Hat-algebra and levels will
% still give the same results from counterfactuals. We will get there.

close all;