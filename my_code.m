%% CLEAN THE WORK SPACE

clear; close all; clc;

%% LOAD DATA

my_data = load('assignment_graph.txt');

%% Question 1

% create edge list
m = max(unique(my_data(:, 1)));
n = max(unique(my_data(:, 2)));

my_graph = sparse(my_data(:, 1), my_data(:, 2), 1, m, n);

% spy plot
spy(my_graph)

%close all;

%% Question 2

% degree per node
degree_per_node = full(sum(my_graph, 2)) .* 2;

% plot deg distribution
uniq_deg = unique(degree_per_node);
count_per_deg = zeros(length(uniq_deg), 1);

for i = 1 : length(uniq_deg)
    count_per_deg(i) = nnz(find(degree_per_node == uniq_deg(i)));
end

figure
scatter(uniq_deg, count_per_deg)
set(gca,'xscale','log', 'yscale', 'log')
xlabel('Degree')
ylabel('Node count')
title('Degree distribution')

%% CLEAR AGAIN

clear my_data;

clear count_per_deg uniq_deg degree_per_node

%% Question 3

[U, S, V] = svds(my_graph, 40);
new_adj_matrix = U * S * V';

diff_mat = 0.0;
ori_mat = 0.0;
for i = 1 : m
    i
    for j = 1 : n
        diff_mat = diff_mat + ( new_adj_matrix(i, j) - ...
            my_graph(i, j) ) ^ 2;
        ori_mat = ori_mat + my_graph(i, j) ^ 2;
    end
end
