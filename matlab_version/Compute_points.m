function points = Compute_points(x0, y0, r)
    points = zeros(8,2);  % 初始化points矩阵
    for i = 0:7
        angle = 2 * pi * i / 8;  % 计算每个点对应的角度
        x = x0 + r * cos(angle);  % 计算x坐标
        y = y0 + r * sin(angle);  % 计算y坐标
        points(i+1, :) = [x, y];  % 将计算的坐标添加到points矩阵中
    end
end