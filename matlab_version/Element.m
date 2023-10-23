classdef Element < handle
    properties
        edge_m_num
        edge_c_num
        node_m
        node_c
        node_m_id
        edge_m
        edge_c
        element_number
        n_beta_M
        n_beta_C
    end
    
    methods
        function obj = Element(edge_m_num, edge_c_num, node_m, node_c, node_m_id, element_number)
            if nargin < 6, element_number = []; end  % 默认值为[]
            obj.edge_m_num = edge_m_num;
            obj.edge_c_num = edge_c_num;
            obj.node_m = node_m;
            obj.node_c = node_c;
            obj.node_m_id = node_m_id;
            obj.element_number = element_number;
            obj.edge_m = cell(1,edge_m_num);
            obj.edge_c = cell(1,edge_c_num);
            for i = 1:edge_m_num
                edge = Edge(node_m(2*i-1),node_m(2*i),node_m(2*mod(i,edge_m_num)+1),node_m(2*mod(i,edge_m_num)+2));
                obj.edge_m{i} = edge;
            end
            for i = 1:edge_c_num
                edge = Edge(node_c(2*i-1),node_c(2*i),node_c(2*mod(i,edge_c_num)+1),node_c(2*mod(i,edge_c_num)+2));
                obj.edge_c{i} = edge;
            end
            if obj.edge_m_num <=5
                obj.n_beta_M = 7;
            elseif obj.edge_m_num<=7
                obj.n_beta_M = 12;
            elseif obj.edge_m_num>7 
                obj.n_beta_M = 18;
            end
            if obj.edge_c_num <=5
                obj.n_beta_C = 7;
            elseif obj.edge_c_num<=7
                obj.n_beta_C = 12;
            elseif obj.edge_c_num>7 
                obj.n_beta_C = 18;
            end
            %obj.n_beta_C = obj.n_beta_M;
        end
    end
end

