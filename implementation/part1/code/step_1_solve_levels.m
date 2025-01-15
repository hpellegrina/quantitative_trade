clc;
clear;

% load fundamentals
load('../output/fundamental_for_simulation.mat');
d     = fundamentals.d; 
I     = fundamentals.I;
T     = fundamentals.T;
N     = fundamentals.N;
coord = fundamentals.coord;

% trade elasticity
theta = 4;

% solve the model
w         = ones(I,1);  % guess
dif_tol   = 1;          % initial tolerance
tol       = 1e-6;       % tolerance parameter
count     = 1;          % count
max_count = 1000;       % max count
adj_exp   = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum   = 0.5;        % update parameter sum

while dif_tol > tol && count < max_count
    
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
    w_new     = w.*dif_adj_N.^adj_exp;
    w_new     = w_new./w_new(1); 
    w         = w*adj_sum + w_new*(1-adj_sum); 

    % update
    dif_tol = dif_tol_N;
    count   = count + 1;

    % display results
    disp(['Tolerance in N = ', num2str(dif_tol_N)]) 
    disp(['Iteration = ', num2str(count)]) 

end

% export data
data = struct(...
        'E', E, ...
        'X', X, ...
        'N', N, ...
        'coord', coord);
save('../output/data_square_world.mat', 'data');

% make real income figure
figure;
Pindex  = Pi.^(-(1/theta));
rincome = w./Pindex';
I_side  = sqrt(I);
rincome_mat = reshape(rincome, [I_side, I_side]);
imagesc(rincome_mat);
colorbar
saveas(gcf, '../output/map_real_income.png');

% make tfp figure
figure;
T_mat = reshape(T, [I_side, I_side]);
imagesc(T_mat);
colorbar
saveas(gcf, '../output/map_tfp.png');

close all;
