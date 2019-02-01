classdef model < handle
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Length      % Length of BioSim (rows, m)
        Width       % Width of BioSim (cols, m)
        Height      % Height of BioSim (m)
        Ta_C        % Air temperature (C)
        To_C        % Outside temperature (C)
        Ma          % Dry air mass (kg)
        Mv          % Water vapor mass (kg)
        RHmax       % Maximum indoor humidity
        wv          % Wind velocity (m/s)
        dt          % Model time step (s)
        
        blocks      % Array of all block objects
        trackTa     % Air temperature monitoring array
        trackTw     % Average water body temperature array
        trackETo    % Reference evapotranspiration (mm/h)
        
        hour        % Tracks the hour of the day
        day         % Tracks the day since starting
        date        % Tracks the actual date
        
        H           % Heat equation handler
        heatperm2   % Heat delivered per m2
    end
    
    properties (Constant)
        Ra = 0.288;     % gas constant for air kJ/(kg-K)
        Rv = 0.463;     % gas constant for water vapor kJ/(kg-K)
        epsilon = 0.622;% epsilon in gas law equations (=Ra/Rv)
        ca = 1.000;     % air specific heat (kJ/kg-deg C)
        cw = 4.186;     % water specific heat (kJ/kg-deg C)
        cv = 1.996;     % water vapor specific heat (kJ/kg-deg C)
        kWall = 0.001;  % thermal conductivity of building wall material (kW/m-K)
        sWall = 0.7;    % wall thickness (m)
        
    end
    
    properties (Dependent)
        Area        % Area of BioSim (m2)
        Volume      % Volume of BioSim (m3)
        P           % Total pressure (kPa)
        Pa          % Dry air pressure (kPa)
        Pv          % Water vapor pressure (kPa)
        Pvs         % Saturation water vapor pressure (kPa)
        M           % Total moist air mass (kg)
        RHO         % Total moist air density (kg/m3)
        RHOa        % Dry air density (kg/m3)
        RHOv        % Water vapor density (kg/m3)
        Ta_K        % Air temperature (K)
        To_K        % Outside air temperature (K)
        RH          % Current relative humidity
        
        waterBlocks % Array of all water blocks
    end
    
    methods
        function obj = model(L, W, H, P, Ta, To, RH, dt)
            %MODEL Construct an instance of this class
            %   A   = Area of BioSim (m2)
            %   H   = Height of BioSim space (m)
            %   P   = Initial pressure (kPa)
            %   Ta  = Initial air temperature (oC)
            %   To  = Outside air temperature (oC)
            %   RH  = Desired indoor humidity
            %   dt  = Model time step (time enum)
            
            if nargin > 0
                % Initialize model
                obj.Length = L;
                obj.Width = W;
                obj.Height = H;
                obj.Ta_C = Ta;
                obj.To_C = To;
                obj.Mv = 0; 
                RHOa = P / ((Ta+273.15) * obj.Ra);
                obj.Ma = RHOa * obj.Volume;
                obj.RHmax = RH;
                obj.wv = 0;
                
                obj.dt = dt;
                
                obj.hour = 1; obj.day = 1;
                obj.date = datetime('today');

                % Create empty (soil) array of blocks to initialize area
                % Blocks are 1x1 soil blocks
                for j = 1:1:W
                    for i = 1:1:L
                        blocks(j, i) = block(j, i, 1, 1, soil(1, 'soil'));
                    end
                end
                obj.blocks = blocks;
            end
        end
        
        function A = get.Area(obj)
            A = obj.Length * obj.Width;
        end
        
        function V = get.Volume(obj)
            V = obj.Area * obj.Height;
        end
        
        function M = get.M(obj)
            M = obj.Mv + obj.Ma; 
        end
        
        function RHO = get.RHO(obj)
            RHO = obj.RHOa + obj.RHOv;
        end
        
        function RHOa = get.RHOa(obj)
            RHOa = obj.Ma / obj.Volume;
        end
        
        function RHOv = get.RHOv(obj)
            RHOv = obj.Mv / obj.Volume;
        end
        
        function P = get.P(obj)
            P = obj.Pa + obj.Pv;
        end
        
        function Pa = get.Pa(obj)
            Pa = obj.RHOa * obj.Ra * obj.Ta_K;
        end
        
        function Pv = get.Pv(obj)
            Pv = obj.RHOv * obj.Rv * obj.Ta_K;
        end
        
        function Pvs = get.Pvs(obj)
            Pvs = refEQ.Pvs(obj.Ta_C);
        end
        
        function RH = get.RH(obj)
            RH = obj.Pv / obj.Pvs;
        end
        
        function Ta_K = get.Ta_K(obj)
            Ta_K = obj.Ta_C + 273.15;
        end
        
        function To_K = get.To_K(obj)
            To_K = obj.To_C + 273.15;
        end
        
        function waterBlocks = get.waterBlocks(obj)
            c = arrayfun(@(b) class(b.data), obj.blocks, 'UniformOutput', false);
            mask = strcmp(c, 'water');
            w = arrayfun(@(b) b.data, obj.blocks(mask));
            waterBlocks = unique(w);
        end
        
        function setup(obj)
            % Create necessary handlers & variables
            obj.H = heat(obj.dt);
            
            % TODO: Put lighting in separate object or area
            % Lighting requirements for living space
            % https://www.engineeringtoolbox.com/light-level-rooms-d_708.html
            lumensperm2 = 1000;
            % Bulbs used: http://www.alconlighting.com/specsheets/alcon/13124%20-%20Sinch.pdf
            lumensperbulb = 1000;
            heatperbulb = 12; % (W)
            % No. of bulbs needed
            no_bulbs = ceil(lumensperm2 * obj.Area / lumensperbulb);
            obj.heatperm2 = no_bulbs * heatperbulb / obj.Area;
            
            % Create a lake object
            w.x=1;
            w.y=1;
            w.w=3;
            w.h=10;
            waterBody = water(w.w*w.h, 1, 20);
            for i = w.x:1:(w.x+w.w-1)
                for j = w.y:1:(w.y+w.h-1)
                obj.blocks(i, j).data = waterBody;
                end
            end
        end
        
        function stepTime(obj)
            %% TODO: Update to change based off dt
            if obj.hour == 24
                obj.hour = 0;
                obj.day = obj.day + 1;
            else
                obj.hour = obj.hour + 1;
            end
            obj.date = obj.date + 1/24;
        end
        
        function output = run(obj, hours)
            for i = hours
                obj.step;
            end
            
            output.Ta = obj.trackTa;
            output.Tw = obj.trackTw;
            output.ETo = obj.trackETo;
        end
        
        function output = step(obj)
        %STEP Run the model for one dt step
            if isempty(obj.H)
                obj.setup;
            end
        
            %% Calculate the amount of heat exchange
            Tw = [];
            QQsn = 0;
            for j = 1:1:(obj.dt*1)
                % Calculate amount of water released into air from
                % all water blocks
                dQ = 0; dQsn = 0; dQan = 0; dQe = 0; dQc = 0; dQbr = 0;
                dMv = 0;
                for w = obj.waterBlocks
                    % Calculate the heat absorbed by water due to lighting
                    % Valid for "daylight" hours (5AM - 8PM, for example)
                    if obj.hour > 5 && obj.hour < 20
                       Qsn = obj.heatperm2 * w.Aw * obj.dt; % (W/m^2) * m^2
                    else
                       Qsn = 0;
                    end

                    % Calculate saturation vapor pressure and density
                    RHOvs = obj.Pvs / (obj.Rv * obj.Ta_K);

                    % Evaporation coefficient (kg/m2-s)
                    ae = (25 + 19*obj.wv) / 3600;
                    % Hypothetical saturated vapor mass per total mass (kg/kg)
                    % i.e. saturation humidity ratio
                    x_s = RHOvs * obj.Volume / obj.Ma;
                    % Actual (current) vapor mass per total mass (kg/kg)
                    % i.e. humidity ratio
                    x = obj.Mv / obj.Ma;
                    % Amount of water vapor evaporated (kg/dt)
                    dmv = w.Aw * (obj.RHmax*x_s - x) * ae * obj.dt;

                    % Solve nonlinear equation for equilibrium temp. of water
                    fun = @(T) obj.H.water(T, obj.Ta_C, obj.Pv, dmv, Qsn, obj.wv, w.Aw);
                    [Ts_C, Q] = fzero(fun, w.Tw_C);

                    % Calculate change in water temp based off surface
                    % temp change
                    ds = 0.001; % Surface depth = 1mm
                    RHOs = refEQ.RHOw(Ts_C);
                    Vs = w.Aw * ds;
                    w.Tw_C = (w.RHOw*(w.Vw-Vs)*w.Tw_C + RHOs*Vs*Ts_C)/(w.RHOw*(w.Vw-Vs) + RHOs*Vs);

                    % Update mass of water
                    w.Mw = w.Mw - dmv;

                    % Track water temp
                    Tw = [Tw w.Tw_C];

                    dQ = dQ + Q/1000;
                    dQsn = dQsn + Qsn/1000;
                    dQan = dQan + obj.H.Qan/1000;
                    dQbr = dQbr + obj.H.Qbr/1000;
                    dQe = dQe + obj.H.Qe/1000;
                    dQc = dQc + obj.H.Qc/1000;
                    dMv = dMv + dmv;
                end
                % Add short-wave radiation to accumulator
                QQsn = QQsn + dQsn;

                % Calculate heat exchange due to walls of BioSim
                qWall = (obj.kWall / obj.sWall) * (2*obj.Length*obj.Height + 2*obj.Width*obj.Height) * (obj.To_C - obj.Ta_C) * obj.dt;

                % Update temperature of air
                obj.Ta_C = (qWall+(dQbr+dQc+dQe-dQan))/(obj.Mv*obj.cv + obj.Ma*obj.ca) + obj.Ta_C;

                % Update mass of water vapor in air
                obj.Mv = obj.Mv + dMv;
            end
            
            ETo = refEQ.ETo(obj.Ta_C, QQsn, obj.wv, obj.P, obj.Pv);
            obj.trackETo = [obj.trackETo ETo];
            obj.trackTa = [obj.trackTa obj.Ta_C];
            obj.trackTw = [obj.trackTw mean(Tw)];
            output.Ta = obj.trackTa;
            output.Tw = obj.trackTw;
            output.ETo = obj.trackETo;
            
            obj.stepTime;
        end
    end
end

