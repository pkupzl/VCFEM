classdef Mesh < handle
    properties
        node_num
        element_num
        nodes
        elements
    end
    
    methods
        function obj = Mesh(node_num, element_num, nodes, edge_m_nums, edge_c_nums, node_ms, node_cs, node_m_ids,node_c_ids)
            obj.node_num = node_num;
            obj.element_num = element_num;
            obj.nodes = nodes;
            obj.elements = cell(1,element_num);
            for i = 1:element_num
                element = Element(edge_m_nums(i),edge_c_nums(i),node_ms{i},node_cs{i},node_m_ids{i},node_c_ids{i},i);
                obj.elements{i} = element;
            end
        end
    end
end

