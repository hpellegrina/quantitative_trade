clc;
clear;

%% GEOGRAPHY
% we will make an artificial world in which we know exactly the
% fundamentals. This will allow us to test if our algorithms are working
% by the end. I recommend that you always do that in your papers before 
% bringing actual data. The world here will have a square shape.

% basic characteristics
I_side      = 10;                          % number of countries in each side of the square
I           = I_side*I_side;              % total number of countries
coord       = zeros(I_side*I_side,2);     % coordinates
matrix_ref  = zeros(I_side, I_side);      % reference matrix
N           = ones(I,1);                  % population

% make tfps
for i = 1:I_side
    for j = 1:I_side
        matrix_ref(i, j) = abs(i - 1) + abs(j - 1) + 1;
    end
end
T_sq = matrix_ref.^(-0.1);               % tfp of countries in a square
T    = T_sq(:);                          % tfp of countries in a vector

% make coordinates
row = 1;
for i = 1:I_side
    for j = 1:I_side
        coord(row,1) = i;
        coord(row,2) = j;
        row = row + 1;
    end
end

% make iceberg trade cost
dist = pdist2(coord, coord, 'euclidean');
d0   = 0.01;
d1   = 0.5;
d    = exp(d0.*dist.^d1);

% export data
fundamentals = struct(...
        'I', I, ...
        'd', d, ...
        'T', T, ...
        'N', N, ...
        'coord', coord);
save('../output/fundamental_for_simulation.mat', 'fundamentals');