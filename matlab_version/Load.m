classdef Load < handle
    properties
        load_num
        type
        q
        node_i
        node_j
    end
    
    methods
        function obj = Load(load_num, type, q, node_i, node_j)
            obj.load_num = load_num;
            obj.type = type;
            obj.q = q;
            obj.node_i = node_i;
            obj.node_j = node_j;
        end
    end
end

