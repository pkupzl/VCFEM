function Visualize(nodes, topPoints, bottomPoints, leftPoints, rightPoints,displacement)
    x_values = zeros(1, length(nodes));
    y_values = zeros(1, length(nodes));
    for i = 1:length(nodes)
        x_values(i) = nodes{i}.x;
        
        y_values(i) = nodes{i}.y;
        if nargin==6
            x_values(i) = x_values(i) + displacement(2*i-1);
            y_values(i) = y_values(i) + displacement(2*i);
        end
    end
    
    % 绘制初始点
    %scatter(x_values, y_values, 'k', 'filled');
    hold on;

    % 绘制上边界
    for i = 1:(length(topPoints)-1)
        plot([x_values(topPoints(i)), x_values(topPoints(i+1))], [y_values(topPoints(i)), y_values(topPoints(i+1))], 'b-');
    end

    % 绘制右边界
    for i = 1:(length(rightPoints)-1)
        plot([x_values(rightPoints(i)), x_values(rightPoints(i+1))], [y_values(rightPoints(i)), y_values(rightPoints(i+1))], 'b-');
    end

    % 绘制下边界
    for i = 1:(length(bottomPoints)-1)
        plot([x_values(bottomPoints(i)), x_values(bottomPoints(i+1))], [y_values(bottomPoints(i)), y_values(bottomPoints(i+1))], 'b-');
    end

    % 绘制左边界
    for i = 1:(length(leftPoints)-1)
        plot([x_values(leftPoints(i)), x_values(leftPoints(i+1))], [y_values(leftPoints(i)), y_values(leftPoints(i+1))], 'b-');
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
