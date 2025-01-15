clc;
clear;

% load data
load('../output/data_square_world.mat');
E_data    = data.E;
X_data    = data.X;
N         = data.N;
I         = length(E_data);

% trade elasticity
theta     = 4;

% let's construct the pi_0 matrix and E_0 first
pi_0 = X_data./sum(X_data,1);
E_0  = E_data;

% counterfactual changes in exogenous parameters
T_hat = ones(length(E_data),1);
T_hat(1) = 1.2;
d_hat = ones(size(X_data));

% initiate loop
w_hat       = ones(I,1);
adj_exp     = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum     = 0.5;        % update parameter sum
dif_tol     = 1;          % initial tolerance
tol         = 1e-6;       % tolerance parameter
count       = 1;          % count
max_count   = 1000;       % max count

while dif_tol > tol && count < max_count

    % ensure first country is a numeraire
    w_hat     = w_hat./w_hat(1);
    E_hat     = w_hat;

    % change in multilateral resistance term
    Phi_hat   = sum(pi_0.*T_hat.*(d_hat.*w_hat).^(-theta),1);

    % change in trade
    X_hat     = ((T_hat.*(d_hat.*w_hat).^(-theta))./Phi_hat).*E_hat';

    % implied change in wages
    w_hat_new = sum(pi_0.*E_0'.*X_hat,2)./E_0;
    dif_adj   = w_hat_new./w_hat;
    dif_tol   = max(abs(1-dif_adj));
    w_hat_new = w_hat.*dif_adj.^adj_exp;    
    w_hat_new = w_hat_new./w_hat_new(1);
    w_hat     = w_hat_new*adj_sum + w_hat*(1-adj_sum);     

    % update
    count     = count + 1;

end

% change in price index
Pindex_hat  = Phi_hat.^(-1/theta);
D_rincome = 100*(w_hat./Pindex_hat'-1);

% make some figure
D_rincome(1) = max(D_rincome(2:end)); % FOR VISUALIZATION WE WILL REMOVE THE GAINS FOR 1
I_side  = sqrt(I);
D_rincome_mat = reshape(D_rincome, [I_side, I_side]);
imagesc(D_rincome_mat);
colorbar
saveas(gcf, '../output/map_D_real_income_hat.png');

% YOU CAN CHECK THE RESULTS ARE EXACTLY THE SAME AS WHEN WE SOLVED IN
% LEVELS

