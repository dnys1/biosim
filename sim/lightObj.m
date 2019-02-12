classdef lightObj < handle
    %LIGHT Summary of this class goes here
    %   Detailed explanation goes here
    % Lighting requirements for living space
    % https://www.engineeringtoolbox.com/light-level-rooms-d_708.html
    % Bulbs used:
    % Living: http://www.alconlighting.com/specsheets/alcon/13124%20-%20Sinch.pdf
    % Growing: https://www.lumigrow.com/wp-content/uploads/2018/11/LumiGrow-TopLight-Specifications_PN-770-00016-C.pdf
    
    properties
        A       % Area of space illuminated by light (m2)
        W       % Heat wattage per bulb (W)
        L       % Lumens per bulb (lumen)
        Lmax    % Desired lumens per m2 (lumen)
        desc    % Textual description of light
        
        onT     % Hour of day to turn on lights (1-23)
        offT    % Hour of day to turn off lights (1-23)
        
        no_bulbs% Number of lightbulbs for this light fixture
        
        model   % Model controlling light
    end
    
    properties (Dependent)
        Q       % Current heat output (W/hr)
        Qperm2  % Heat per m2 (W/m2)
    end
    
    methods
        function obj = lightObj(A,W,L,Lmax,model,onT,offT,desc)
            %LIGHT Construct an instance of this class
            %   Detailed explanation goes here
            obj.A = A;
            obj.W = W;
            obj.L = L;
            obj.Lmax = Lmax;
            obj.model = model;
            obj.onT = onT;
            obj.offT = offT;
            obj.desc = desc;
            
            obj.no_bulbs = ceil(Lmax * A / L);
        end
        
        function Qperm2 = get.Qperm2(obj)
            if obj.model.hour >= obj.onT && obj.model.hour <= obj.offT
                Qperm2 = obj.no_bulbs * obj.W / obj.A * obj.model.dt;
            else
                Qperm2 = 0;
            end
        end
    end
end

