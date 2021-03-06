classdef waterBody < handle
    %WATERBODY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Aw      % Area of the water body (m2)
        hw      % Height of the water body (m)
        Tw_C    % Temperature of the water body (oC)
    end
    
    properties (Dependent)
        Tw_K    % Temperature of the water body (K)
        RHOw    % Water body density (kg/m3)
        Vw      % Volume of water (m3)
        Mw      % Mass of waterbody (kg)
    end
    
    methods
        function obj = waterBody(Aw, hw, Tw_C)
            %WATERBODY Construct an instance of this class
            %   Detailed explanation goes here
            obj.Aw = Aw;
            obj.hw = hw;
            obj.Tw_C = Tw_C;
        end
        
        function Tw_K = get.Tw_K(obj)
            Tw_K = obj.Tw_C + 273.15;
        end
        
        function RHOw = get.RHOw(obj)
            RHOw = refEQ.RHOw(obj.Tw_C);
        end
        
        function Vw = get.Vw(obj)
            Vw = obj.Aw * obj.hw;
        end
        
        function Mw = get.Mw(obj)
            Mw = obj.Vw * obj.RHOw;
        end
    end
end

