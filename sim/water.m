classdef water < handle
    %WATER Tracks the flow of water throughout the system
    %   Water is compartmentalized in this class into
    %   different areas throughout the system.
    
    properties
        M           % Total mass of water in system (kg)
        Mtank       % Mass of water in storage tank (kg)
        Msoil       % Mass of water in soil (kg)
        Mvfw        % Mass of water in vertical wetlands (kg)
        Mvapor      % Mass of water in air as vapor (kg)
        Mbody       % Mass of water in water bodies (kg)
    end
    
    methods
        function obj = water(Mtank, Mbody)
            %WATER Construct an instance of this class
            %   Mtank - initial water amount in storage
            %   Mbody - initial water amount in open water bodies
           obj.M = Mtank + Mbody;
           
           % Make some initial assumptions about the distribution of water
           % in the system
           obj.Mtank = Mtank;
           obj.Mbody = Mbody;
           obj.Mvapor = 0;
           obj.Msoil = 0;
           obj.Mvfw = 0;
        end
    end
end

