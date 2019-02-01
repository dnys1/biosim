classdef soil < handle
    %SOIL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        As      % Area of soil
        type    % Type of soil
    end
    
    methods
        function obj = soil(As, type)
            %SOIL Construct an instance of this class
            %   Detailed explanation goes here
            obj.As = As;
            obj.type = type;
        end
    end
end

