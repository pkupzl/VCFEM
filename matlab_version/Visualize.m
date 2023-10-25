function Visualize(color,mesh,nodes, topPoints, bottomPoints, leftPoints, rightPoints,displacement,particle_nodes)
    x_values = zeros(1, length(nodes));
    y_values = zeros(1, length(nodes));
    for i = 1:length(nodes)
        x_values(i) = nodes{i}.x;
        y_values(i) = nodes{i}.y;
        if nargin>=8
            x_values(i) = x_values(i) + displacement(2*i-1);
            y_values(i) = y_values(i) + displacement(2*i);
        end
    end

    % 绘制初始点
    %scatter(x_values, y_values, 'k', 'filled');
    hold on;

    % 绘制上边界
    for i = 1:(length(topPoints)-1)
        plot([x_values(topPoints(i)), x_values(topPoints(i+1))], [y_values(topPoints(i)), y_values(topPoints(i+1))], color);
    end

    % 绘制右边界
    for i = 1:(length(rightPoints)-1)
        plot([x_values(rightPoints(i)), x_values(rightPoints(i+1))], [y_values(rightPoints(i)), y_values(rightPoints(i+1))], color);
    end

    % 绘制下边界
    for i = 1:(length(bottomPoints)-1)
        plot([x_values(bottomPoints(i)), x_values(bottomPoints(i+1))], [y_values(bottomPoints(i)), y_values(bottomPoints(i+1))], color);
    end

    % 绘制左边界
    for i = 1:(length(leftPoints)-1)
        plot([x_values(leftPoints(i)), x_values(leftPoints(i+1))], [y_values(leftPoints(i)), y_values(leftPoints(i+1))], color);
    end
    if nargin>=9
        %位移前颗粒
        for i = 1:mesh.element_num
            element = mesh.elements{i};
            len = length(element.node_c_id);
            for j = 1:(len-1)
                plot([particle_nodes{element.node_c_id(j)}.x, particle_nodes{element.node_c_id(j+1)}.x], [particle_nodes{element.node_c_id(j)}.y, particle_nodes{element.node_c_id(j+1)}.y], 'g--');
            end
            plot([particle_nodes{element.node_c_id(len)}.x, particle_nodes{element.node_c_id(1)}.x], [particle_nodes{element.node_c_id(len)}.y, particle_nodes{element.node_c_id(1)}.y], 'g--');
        end
        %位移后颗粒
        for i = 1:mesh.element_num
            element = mesh.elements{i};
            len = length(element.node_c_id);
            d_c = element.d_c;
            for j = 1:(len-1)
                plot([particle_nodes{element.node_c_id(j)}.x+d_c(2*j-1), particle_nodes{element.node_c_id(j+1)}.x+d_c(2*j+1)], [particle_nodes{element.node_c_id(j)}.y+d_c(2*j), particle_nodes{element.node_c_id(j+1)}.y+d_c(2*j+2)], 'b-');
            end
            plot([particle_nodes{element.node_c_id(len)}.x+d_c(2*len-1), particle_nodes{element.node_c_id(1)}.x+d_c(1)], [particle_nodes{element.node_c_id(len)}.y+d_c(2*len), particle_nodes{element.node_c_id(1)}.y+d_c(2)], 'b-');
        end

    end

    

    % 绘制闭合的正方形边界
    %plot([x_values(topPoints(1)), x_values(rightPoints(1)), x_values(bottomPoints(1)), x_values(leftPoints(1)), x_values(topPoints(1))], ...
    %     [y_values(topPoints(1)), y_values(rightPoints(1)), y_values(bottomPoints(1)), y_values(leftPoints(1)), y_values(topPoints(1))], 'b-');

    % 标记用到的点
    used_points = unique([topPoints, bottomPoints, leftPoints, rightPoints]);
    scatter(x_values(used_points), y_values(used_points), 'r', 's', 'filled');

    % 设置坐标轴范围
    %xlim([min(x_values)-1, max(x_values)+1]);
    %ylim([min(y_values)-1, max(y_values)+1]);

    axis equal;

    hold off;
end
