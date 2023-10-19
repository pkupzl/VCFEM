classdef Displacement_BC < handle
    properties
        num
        node
        type
        d
    end
    
    methods
        function obj = Displacement_BC(num, node, type, d)
            obj.num = num;
            obj.node = node;
            obj.type = type;
            obj.d = d;
        end
    end
end

