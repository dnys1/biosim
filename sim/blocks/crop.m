classdef crop < handle
    %CROP Summary of this class goes here
    %   https://aces.nmsu.edu/aes/irrigation/wheat.html
    
    properties
        type        % Type of crop being grown (i.e. wheat)
        Ac          % Area of the crop (m2)
        
        model       % A handle for the model
    end
    
    properties (Dependent)
        Kc          % Crop coefficient
        yield       % Yield if harvested now (kg)
        LAI         % Leaf area index
        Tbase       % Base temperature
        Tcutoff     % Cutoff temperature
        GDD         % Growing degree days (oC-d)
    end
    
    methods
        function obj = crop(model, type, Ac)
            %CROP Construct an instance of this class
            %   Detailed explanation goes here
            if (strcmpi(type, 'Wheat') || ...
                strcmpi(type, 'Corn') || ...
                strcmpi(type, 'Soy'))
                obj.type = type;
                obj.Ac = Ac;
                obj.model = model;
            else
                error('Invalid crop type')
            end
        end
        
        function GDD = get.GDD(obj)
            T = min(obj.model.trackTa, obj.Tcutoff) - obj.Tbase;
            T = T(T>0);
            GDD = max(trapz(T) / 24, 0);
        end
        
        function Kc = get.Kc(obj)
           % https://aces.nmsu.edu/aes/irrigation/wheat.html
           switch(obj.type)
               case 'Wheat'
                   Kc = 2.7e-1 - 4.8e-4*obj.GDD + 6.27e-7*obj.GDD - 1.3e-10*obj.GDD^3;
               case 'Corn'
                   Kc = 1.2e-1 + 1.68e-3*obj.GDD - 2.46e-7*obj.GDD^2 - 4.37e-10*obj.GDD^3;
               case 'Soy'
                   Kc = 4.35e-2 + 1.37e-3*obj.GDD - 5.3e-7*obj.GDD^2 + 5.43e-11*obj.GDD^3;
           end
        end
        
        function Tbase = get.Tbase(obj)
            switch(obj.type)
                case 'Wheat'
                    Tbase = 0;
                case 'Corn'
                    Tbase = 10;
                case 'Soy'
                    Tbase = 10;
            end
        end
        
        function Tcutoff = get.Tcutoff(obj)
            switch(obj.type)
                case 'Wheat'
                    Tcutoff = 30;
                case 'Corn'
                    Tcutoff = 30;
                case 'Soy'
                    Tcutoff = 30;
            end
        end
        
        function LAI = get.LAI(obj)
            %TODO
           LAI = obj.GDD; 
        end
    end
end

