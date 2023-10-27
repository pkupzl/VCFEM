function opt = Initialize_from_csv(datapath,scale_factor,type)
    % 读取保存的数据
    matrix_cleaned_nodes_filename = strcat(datapath,'matrix_cleaned_nodes.csv');
    matrix_polygon_node_indices_filename = strcat(datapath,'matrix_polygon_node_indices.csv');
    matrix_cleaned_nodes = csvread(matrix_cleaned_nodes_filename);
    fid = fopen(matrix_polygon_node_indices_filename, 'r');
    matrix_dataArray = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    matrix_polygon_node_indices = matrix_dataArray{1};
    particle_cleaned_nodes_filename = strcat(datapath,'particle_cleaned_nodes.csv');
    particle_polygon_node_indices_filename = strcat(datapath,'particle_polygon_node_indices.csv');
    particle_cleaned_nodes = csvread(particle_cleaned_nodes_filename);
    fid = fopen(particle_polygon_node_indices_filename, 'r');
    particle_dataArray = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    particle_polygon_node_indices = particle_dataArray{1};
    [opt.node_num,~] = size(matrix_cleaned_nodes);
    [opt.particle_node_num,~] = size(particle_cleaned_nodes);
    [opt.element_num,~] = size(matrix_polygon_node_indices);
    
    % 创建nodes,基体节点坐标总表
    opt.nodes = cell(1, opt.node_num);
    for i=1:opt.node_num
        opt.nodes{i} = Node(matrix_cleaned_nodes(i,1)*scale_factor,(1024-matrix_cleaned_nodes(i,2))*scale_factor);
    end
    % 创建particle_nodes，颗粒节点坐标总表
    opt.particle_nodes = cell(1,opt.particle_node_num);
    for i=1:opt.particle_node_num
        opt.particle_nodes{i} = Node(particle_cleaned_nodes(i,1)*scale_factor,(1024-particle_cleaned_nodes(i,2))*scale_factor);
    end

    opt.edge_m_nums = zeros(1,opt.element_num);
    opt.edge_c_nums = zeros(1,opt.element_num);
    % 创建node_ms
    opt.node_ms = cell(1,opt.element_num);
    opt.node_cs = cell(1,opt.element_num);
    opt.node_m_ids = cell(1,opt.element_num);
    for k = 1:opt.element_num
        matrix_str = matrix_polygon_node_indices{k};
        particle_str = particle_polygon_node_indices{k};
        matrix_strArray = strsplit(matrix_str, ',');
        particle_strArray = strsplit(particle_str,',');
        matrix_indice = cellfun(@str2num, matrix_strArray)+1;
        opt.node_m_ids{k} = matrix_indice;
        particle_indice = cellfun(@str2num, particle_strArray)+1;
        opt.node_c_ids{k} = particle_indice;
        [~,edge_m_num] = size(matrix_indice);
        [~,edge_c_num] = size(particle_indice);
        opt.edge_m_nums(k) = edge_m_num;
        opt.edge_c_nums(k) = edge_c_num;
        opt.node_ms{k} = zeros(1,2*edge_m_num);
        opt.node_cs{k} = zeros(1,2*edge_c_num);
        %disp(k);
        for i = 1:edge_m_num
            opt.node_ms{k}(2 * i - 1:2 * i) = [opt.nodes{matrix_indice(i)}.x, opt.nodes{matrix_indice(i)}.y];
        end
        for i = 1:edge_c_num
            opt.node_cs{k}(2 * i - 1:2 * i) = [opt.particle_nodes{particle_indice(i)}.x, opt.particle_nodes{particle_indice(i)}.y];
        end
        
    end
    opt.leftPoints = [];
    opt.rightPoints = [];
    opt.topPoints = [];
    opt.bottomPoints = [];
    for i = 1:length(opt.nodes)
        point = opt.nodes{i};
        x = point.x;
        y = point.y;
        % 判断点是否接近边界，可以自定义阈值
        threshold = 10*scale_factor;  % 阈值用于判断接近边界的距离
        if x < threshold
            opt.leftPoints = [opt.leftPoints, i];
        end
        if x > 1-threshold
            opt.rightPoints = [opt.rightPoints, i];
        end
        if y < threshold
            opt.bottomPoints = [opt.bottomPoints, i];
        end
        if y > 1-threshold
            opt.topPoints = [opt.topPoints, i];
        end
    end
    % 获取topPoints中每个点的x坐标
    xCoords = cellfun(@(i) opt.nodes{i}.x, num2cell(opt.topPoints));
    % 根据x坐标进行排序，并获取排序后的序号
    [~, sortedIndices] = sort(xCoords,'descend');
    % 根据排序后的序号重新排列topPoints数组
    sortedTopPoints = opt.topPoints(sortedIndices);
    opt.topPoints = sortedTopPoints;
    % 获取bottomPoints中每个点的x坐标
    xCoords = cellfun(@(i) opt.nodes{i}.x, num2cell(opt.bottomPoints));
    % 根据x坐标进行排序，并获取排序后的序号
    [~, sortedIndices] = sort(xCoords);
    % 根据排序后的序号重新排列bottomPoints数组
    sortedBottomPoints = opt.bottomPoints(sortedIndices);
    opt.bottomPoints = sortedBottomPoints;
    % 获取leftPoints中每个点的y坐标
    yCoords = cellfun(@(i) opt.nodes{i}.y, num2cell(opt.leftPoints));
    % 根据y坐标进行排序，并获取排序后的序号
    [~, sortedIndices] = sort(yCoords,'descend');
    % 根据排序后的序号重新排列leftPoints数组
    sortedLeftPoints = opt.leftPoints(sortedIndices);
    opt.leftPoints = sortedLeftPoints;
    % 获取rightPoints中每个点的y坐标
    yCoords = cellfun(@(i) opt.nodes{i}.y, num2cell(opt.rightPoints));
    % 根据y坐标进行排序，并获取排序后的序号
    [~, sortedIndices] = sort(yCoords);
    % 根据排序后的序号重新排列rightPoints数组
    sortedRightPoints = opt.rightPoints(sortedIndices);
    opt.rightPoints = sortedRightPoints;
    if strcmp(type,'上下均布')
        %假设在上边界施加载荷
        opt.load_num = length(opt.topPoints)-1;
        load_type_str = '作用均布载荷';
        opt.load_type = cellstr(repmat(load_type_str, opt.load_num, 1));
        q_value = -10;
        opt.q = repmat(q_value, 1, opt.load_num);
        opt.node_i = opt.topPoints(1:end-1);
        opt.node_j = opt.topPoints(2:end);
        
        opt.dbc_num = 2*length(opt.bottomPoints);
        opt.dbc_node = repelem(opt.bottomPoints, 2);
        opt.dbc_type = cell(1, opt.dbc_num);
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 0) = {'x位移'};
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 1) = {'y位移'};
        %opt.d = -1e-5 * ones(1, opt.dbc_num);
        opt.d = zeros(1,opt.dbc_num);
    elseif strcmp(type,'左右均布')
        %假设在右边界施加载荷
        opt.load_num = length(opt.rightPoints)-1;
        load_type_str = '作用均布载荷';
        opt.load_type = cellstr(repmat(load_type_str, opt.load_num, 1));
        q_value = -10;
        opt.q = repmat(q_value, 1, opt.load_num);
        opt.node_i = opt.rightPoints(1:end-1);
        opt.node_j = opt.rightPoints(2:end);
        
        opt.dbc_num = 2*length(opt.leftPoints);
        opt.dbc_node = repelem(opt.leftPoints, 2);
        opt.dbc_type = cell(1, opt.dbc_num);
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 0) = {'x位移'};
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 1) = {'y位移'};
        %opt.d = -1e-5 * ones(1, opt.dbc_num);
        opt.d = zeros(1,opt.dbc_num);
    elseif strcmp(type,'双轴拉伸')
        opt.load_num = length(opt.rightPoints)+length(opt.topPoints)+length(opt.leftPoints)+length(opt.bottomPoints)-4;
        load_type_str = '作用均布载荷';
        opt.load_type = cellstr(repmat(load_type_str, opt.load_num, 1));
        q_value = -10;
        opt.q = repmat(q_value, 1, opt.load_num);
        opt.node_i = [opt.rightPoints(1:end-1),opt.topPoints(1:end-1),opt.leftPoints(1:end-1),opt.bottomPoints(1:end-1)];
        opt.node_j = [opt.rightPoints(2:end),opt.topPoints(2:end),opt.leftPoints(2:end),opt.bottomPoints(2:end)];
        opt.dbc_num = 2*length(opt.bottomPoints);
        opt.dbc_node = repelem(opt.bottomPoints, 2);
        opt.dbc_type = cell(1, opt.dbc_num);
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 0) = {'x位移'};
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 1) = {'y位移'};
        %opt.d = -1e-5 * ones(1, opt.dbc_num);
        opt.d = zeros(1,opt.dbc_num);
    elseif strcmp(type,'上边左切')
        %假设在上边界施加载荷
        opt.load_num = length(opt.topPoints)-1;
        load_type_str = '作用均布切向载荷';
        opt.load_type = cellstr(repmat(load_type_str, opt.load_num, 1));
        q_value = 10;
        opt.q = repmat(q_value, 1, opt.load_num);
        opt.node_i = opt.topPoints(1:end-1);
        opt.node_j = opt.topPoints(2:end);
        
        opt.dbc_num = 2*length(opt.bottomPoints);
        opt.dbc_node = repelem(opt.bottomPoints, 2);
        opt.dbc_type = cell(1, opt.dbc_num);
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 0) = {'x位移'};
        opt.dbc_type(mod(0:opt.dbc_num-1, 2) == 1) = {'y位移'};
        %opt.d = -1e-5 * ones(1, opt.dbc_num);
        opt.d = zeros(1,opt.dbc_num);
    end
    opt.E_m = 1000;
    opt.E_c = 3000;
    opt.pr_m = 0.2;
    opt.pr_c = 0.2;
end