classdef Node < handle
    properties
        x
        y
    end
    
    methods
        function obj = Node(x, y)
            obj.x = x;
            obj.y = y;
        end
    end
end

