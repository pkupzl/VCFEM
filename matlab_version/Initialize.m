function opt = Initialize()
    opt.node_num = 8;
    opt.element_num = 1;
    
    % 创建nodes
    opt.nodes = cell(1, opt.node_num);
    opt.nodes{1} = Node(0, 0);
    opt.nodes{2} = Node(1.25, 0);
    opt.nodes{3} = Node(2.5, 0);
    opt.nodes{4} = Node(2.5, 1.25);
    opt.nodes{5} = Node(2.5, 2.5);
    opt.nodes{6} = Node(1.25, 2.5);
    opt.nodes{7} = Node(0, 2.5);
    opt.nodes{8} = Node(0, 1.25);
    opt.topPoints = [5,6,7];
    opt.bottomPoints = [1,2,3];
    opt.rightPoints = [3,4,5];
    opt.leftPoints = [7,8,1];
    % 创建edge_ms
    opt.edge_ms = cell(1, opt.node_num);
    for i = 1:opt.node_num
        opt.edge_ms{i} = Edge(opt.nodes{i}.x, opt.nodes{i}.y, opt.nodes{mod(i, opt.node_num) + 1}.x, opt.nodes{mod(i, opt.node_num) + 1}.y);
    end
    
    opt.edge_m_nums = [8];
    opt.edge_c_nums = [8];
    
    % 创建edge_cs
    points = Compute_points(1.25, 1.25, 0.7);
    opt.edge_cs = cell(1, length(points));
    opt.particle_nodes = cell(1,length(points));
    for i = 1:length(points)
        opt.edge_cs{i} = Edge(points(i, 1), points(i, 2), points(mod(i - 1, length(points)) + 1, 1), points(mod(i - 1, length(points)) + 1, 2));
        opt.particle_nodes{i} = Node(points(i,1),points(i,2));
    end
    
    % 创建node_ms
    opt.node_ms = cell(1,opt.element_num);
    opt.node_ms{1} = zeros(1, 2 * opt.node_num);
    for k = 1:opt.element_num
        for i = 1:opt.node_num
            opt.node_ms{k}(2 * i - 1:2 * i) = [opt.nodes{i}.x, opt.nodes{i}.y];
        end
    end
    
    % 创建node_cs
    opt.node_cs = cell(1,opt.element_num);
    opt.node_cs{1} = zeros(1, 2 * length(points));
    for k = 1:opt.element_num
        for i = 1:length(points)
            opt.node_cs{k}(2 * i - 1:2 * i) = [points(i, 1), points(i, 2)];
        end
    end
    
    %opt.node_m_ids = 0:opt.node_num-1;
    opt.node_m_ids = cell(1,opt.element_num);
    opt.node_m_ids{1} = 1:opt.node_num;
    opt.node_c_ids{1} = 1:opt.node_num;
    opt.load_num = 2;
    opt.load_type = {'作用均布载荷', '作用均布载荷'};
    opt.q = [-10.0, -10.0];
    %opt.node_i = [4, 5];
    %opt.node_j = [5, 6];
    opt.node_i = [5,6];
    opt.node_j = [6,7];
    opt.dbc_num = 6;
    opt.dbc_node = [1,1,2,2,3,3];
    opt.dbc_type = cell(1, opt.dbc_num);
    for i = 1:opt.dbc_num
        if mod(i-1, 2) == 0
            opt.dbc_type{i} = 'x位移';
        else
            opt.dbc_type{i} = 'y位移';
        end
    end
    %opt.d = -1e-5 * ones(1, opt.dbc_num);
    opt.d = zeros(1,opt.dbc_num);
    opt.E_m = 1000;
    opt.E_c = 3000;
    opt.pr_m = 0.2;
    opt.pr_c = 0.2;
end
