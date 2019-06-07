% author: Yujia Zhai
% This is for CS 235 at UC, Riverside in 19 spring.

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


large_uniq_deg = uniq_deg(uniq_deg > 100);
abnormal_node = [];
for i = 1 : length(large_uniq_deg)
    curr_deg = large_uniq_deg(i);
    curr_deg_set = find(degree_per_node == curr_deg);
    if nnz(curr_deg_set) > 5
        abnormal_node = [abnormal_node; curr_deg_set];
    end
end

sort_abn_set = sort(abnormal_node);
tmp_new_graph = my_graph;
sort_normal_set = setdiff(1:length(degree_per_node), sort_abn_set);
tmp_new_graph(sort_normal_set, :) = 0;
tmp_new_graph(:, sort_normal_set) = 0;
figure; spy(tmp_new_graph)
title('extracted abnormal blocks according to degree distribution')
%% CLEAR AGAIN

%clear my_data;

%clear count_per_deg uniq_deg degree_per_node

%% Question 3

% a
k = 40;
[U, S, V] = svds(my_graph, k);

ori_svd = diag(S);
total_norm = norm(ori_svd, 2);

for i = 1 : k
    curr_norm = norm(ori_svd(1:i), 2);
    if (norm(ori_svd(i+1, :), 2)) / total_norm <= 0.1
        break;
    end
end

msg = sprintf('# of eigs is %d, \n loss rate is %f\n', i, ...
    (norm(ori_svd(i+1, :), 2)) / total_norm);

disp(msg)


% b
graph_length = length(U);

figure
hold on
color_box = ['r', 'b', 'g', 'k', 'm'];

for i = 1 : 5
    plot(1:graph_length, U(:,i), color_box(i))
end
xlim([1 graph_length])

legend('sing vec 1', 'sing vec 2', 'sing vec 3',...
    'sing vec 4', 'sing vec 5')

xlabel('node index')
ylabel('singular vector')
title('plot for first 5 singular vectors')

% d
idx_max = zeros(5, 100);
for i = 1 : 5
    [~, idx_max(i, :)] = maxk(U(:,i), 100);
end

row_nnz = unique(my_data(:, 1));

figure;
for i = 1 : 5
    new_graph = my_graph;
    tmp = setdiff(row_nnz, idx_max(i, :));
    new_graph(tmp, :) = 0;
    new_graph(:, tmp) = 0;
    subplot(2, 3, i)
    spy(new_graph)
    msg = sprintf('singular value %d', i);
    title(msg)
end





