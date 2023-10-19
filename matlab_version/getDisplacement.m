function displacement = getDisplacement(x, y, d_m)
    % 检测是否在正方形边界上
    if (x == 0 && y >= 0 && y <= 2.5) || (x == 2.5 && y >= 0 && y <= 2.5) || (y == 0 && x >= 0 && x <= 2.5) || (y == 2.5 && x >= 0 && x <= 2.5)
        % 在边界上，进行插值计算
        if x == 0
            % 左边界
            if y <= 1.25
                u1 = d_m(1:2);
                u2 = d_m(15:16);
                displacement = (y / 1.25) * u2 + (1 - y / 1.25) * u1;
            else
                u1 = d_m(15:16);
                u2 = d_m(13:14);
                displacement = ((y - 1.25) / 1.25) * u2 + (1 - (y - 1.25) / 1.25) * u1;
            end
        elseif x == 2.5
            % 右边界
            if y <= 1.25
                u1 = d_m(5:6);
                u2 = d_m(7:8);
                displacement = (y / 1.25) * u2 + (1 - y / 1.25) * u1;
            else
                u1 = d_m(7:8);
                u2 = d_m(9:10);
                displacement = ((y - 1.25) / 1.25) * u2 + (1 - (y - 1.25) / 1.25) * u1;
            end
        elseif y == 0
            % 下边界
            if x <= 1.25
                u1 = d_m(1:2);
                u2 = d_m(13:14);
                displacement = (x / 1.25) * u2 + (1 - x / 1.25) * u1;
            else
                u1 = d_m(13:14);
                u2 = d_m(7:8);
                displacement = ((x - 1.25) / 1.25) * u2 + (1 - (x - 1.25) / 1.25) * u1;
            end
        elseif y == 2.5
            % 上边界
            if x <= 1.25
                u1 = d_m(13:14);
                u2 = d_m(11:12);
                displacement = (x / 1.25) * u2 + (1 - x / 1.25) * u1;
            else
                u1 = d_m(11:12);
                u2 = d_m(9:10);
                displacement = ((x - 1.25) / 1.25) * u2 + (1 - (x - 1.25) / 1.25) * u1;
            end
        end
    else
        % 不在边界上
        displacement = NaN;
        disp('点不在正方形的边上');
    end
end