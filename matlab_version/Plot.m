% 拆分位移向量
d_m_x = d_m(1:2:end)';
d_m_y = d_m(2:2:end)';
node_matrix = zeros(2,8);
for i = 1:8
    node_matrix(1,i) = opt.nodes{i}.x;
    node_matrix(2,i) = opt.nodes{i}.y;
end
% 计算变形后的坐标
nodes_def = node_matrix + [d_m_x; d_m_y];

% 插值并画图
hold on;
for i = 1:8
    next_i = mod(i, 8) + 1;
    x = linspace(node_matrix(1,i), node_matrix(1,next_i), 100);
    y = linspace(node_matrix(2,i), node_matrix(2,next_i), 100);
    x_def = linspace(nodes_def(1,i), nodes_def(1,next_i), 100);
    y_def = linspace(nodes_def(2,i), nodes_def(2,next_i), 100);
    plot(x_def, y_def, 'k-');
end
axis equal;