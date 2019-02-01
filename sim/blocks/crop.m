classdef crop < handle
    %CROP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type        % Type of crop being grown (i.e. wheat)
        GDD         % Growing degree days (oC-d)
    end
    
    properties (Dependent)
        Kc          % Crop coefficient
        LAI         % Leaf area index
    end
    
    methods
        function obj = crop(type)
            %CROP Construct an instance of this class
            %   Detailed explanation goes here
            if (strcmpi(type, 'Wheat') || ...
                strcmpi(type, 'Corn') || ...
                strcmpi(type, 'Soy'))
                obj.type = type;
            else
                error('Invalid crop type')
            end
            obj.GDD = 0;
        end
        
        function Kc = get.Kc(obj)
            %TODO
           Kc = obj.GDD; 
        end
        
        function LAI = get.LAI(obj)
            %TODO
           LAI = obj.GDD; 
        end
    end
end

