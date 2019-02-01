classdef block
    %BLOCK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x       % The x-coordinate of the block
        y       % The y-coordinate of the block
        w       % The width (x-dir) of the block
        h       % The height (y-dir) of the block
        data    % The block data (i.e. waterBody, crop, etc)
    end
    
    methods
        function obj = block(x, y, w, h, data)
            %BLOCK Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin > 0
                obj.x=x;
                obj.y=y;
                obj.w=w;
                obj.h=h;
                obj.data=data;
            end
        end
    end
end

