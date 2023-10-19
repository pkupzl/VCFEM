classdef Edge < handle
    properties
        x1
        y1
        x2
        y2
        length
        n1
        n2
    end
    
    methods
        function obj = Edge(x1, y1, x2, y2)
            obj.x1 = x1;
            obj.y1 = y1;
            obj.x2 = x2;
            obj.y2 = y2;
            dx = x2 - x1;
            dy = y2 - y1;
            obj.length = sqrt(dx^2 + dy^2);
            obj.n1 = dy / obj.length;
            obj.n2 = -dx / obj.length;
        end
    end
end

