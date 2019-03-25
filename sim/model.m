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
        RHmax       % Maximum indoor humidity
        wv          % Wind velocity (m/s)
        dt          % Model time step (s)
        no_people   % Number of people BioSim supports
        
        trackTa     % Air temperature monitoring array
        trackTw     % Average water body temperature array
        trackETo    % Reference evapotranspiration (mm/h)
        
        % Tracking of water cycling
        trackVcrop
        trackVwetland
        trackVvapor
        trackVstore
        trackVsed
        trackVreject
        
        % Tracking of crop vars
        trackLeafB  % Array of leafB vals
        trackStorB  % Array of storB vals
        trackStemB  % Array of stemB vals
        trackRootB  % Array of rootB vals
        trackDVS    % Array of DVS vals
        
        % Tracking of lighting energy vars
        tracklLight
        trackgLight
        
        % Tracking of treatment energy vars
        trackRO
        trackUV
        
        hour        % Tracks the hour of the day
        day         % Tracks the day since starting
        date        % Tracks the actual date
        hours       % Array of hours used for plotting
        
        H           % Heat equation handler
        heatperm2   % Heat delivered per m2
        QQsn        % Total heat delivered for this timestep (kW)
        lLight      % lightObj for living space
        gLight      % lightObj for growing space
        
        W           % Water handler
        
        waterBlock  % waterBody handle
        cropBlock   % crop handle
        filterBlock % filterUnit handle
        wetlandBlock% wetland handle
        sedBlock    % sedimentBasin handle
        
        % Water parameters that change each step
        dMv         % delta Mv i.e. change in vapor mass (kg)
        Vi          % irrigation needed for crops (m3)
        
        % Crop parameters for tracking calorie supply/demand
        currentSupply       % Total amount of biomass available (kg)
        currentHourlyDemand % Demand for biomass every hour (kg/hr)
    end
    
    properties (Constant)
        epsilon = 0.622;% epsilon in gas law equations (=Ra/Rv)
        ca = 1.000;     % air specific heat (kJ/kg-deg C)
        cw = 4.186;     % water specific heat (kJ/kg-deg C)
        cv = 1.996;     % water vapor specific heat (kJ/kg-deg C)
        kWall = 0.001;  % thermal conductivity of building wall material (kW/m-K)
        sWall = 0.7;    % wall thickness (m)
        
        hw = 3;         % Height of water bodies (m)
        
        % Water quality parameters
        WW_dp=[12,37.5,468.5,1000]*1e-6         % Raw wastewater particle diameters (m)
        WW_dist=[0.26,0.27,0.33,0.14]           % Raw wastewater particle distribution (-)
        WW_RHOp=1050                            % Raw wastewater particle density (kg/m3)
        WW_TSS_0=1000                           % Raw wastewater TSS concentration (mg/L)
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
        RHOvs       % Saturation vapor density (kg/m3)
        Ta_K        % Air temperature (K)
        To_K        % Outside air temperature (K)
        RH          % Current relative humidity
        ETo         % Reference evapotranspiration (mm/hr)
        Rs          % Crop "solar" radiation (MJ)
        
        WW_Q        % Wastewater flow (m3/d)
    end
    
    methods
        function obj = model(no_people, P, Ta, To, RH, dt)
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
                [At, Ac, Aw] = model.determineArea(no_people);
                Area = At + Ac + Aw;
                obj.Length = sqrt(Area);
                obj.Width = Area / obj.Length;
                obj.Height = 10;
                obj.Ta_C = Ta;
                obj.To_C = To; 
                RHOa = P / ((Ta+273.15) * refEQ.Ra);
                obj.Ma = RHOa * obj.Volume;
                obj.RHmax = RH;
                obj.wv = 0.1;
                obj.no_people=no_people;
                
                obj.dt = dt;
                
                obj.hour = 1; obj.day = 1;
                obj.date = datetime('today');
                obj.hours = [];

                % Create empty (soil) array of blocks to initialize area
                % Blocks are 1x1 soil blocks
%                 for j = 1:1:obj.Width
%                     for i = 1:1:obj.Length
%                         blocks(j, i) = block(j, i, 1, 1, soil(1, 'soil'));
%                     end
%                 end
%                 obj.blocks = blocks;
                
                obj.setup(Ac, Aw);
            end
        end
        
        function A = get.Area(obj)
            A = obj.Length * obj.Width;
        end
        
        function V = get.Volume(obj)
            V = obj.Area * obj.Height;
        end
        
        function M = get.M(obj)
            M = obj.W.Mvapor + obj.Ma; 
        end
        
        function RHO = get.RHO(obj)
            RHO = obj.RHOa + obj.RHOv;
        end
        
        function RHOa = get.RHOa(obj)
            RHOa = obj.Ma / obj.Volume;
        end
        
        function RHOv = get.RHOv(obj)
            RHOv = obj.W.Mvapor / obj.Volume;
        end
        
        function RHOvs = get.RHOvs(obj)
            RHOvs = refEQ.Pvs(obj.Ta_C) / (refEQ.Rv * obj.Ta_K);
        end
        
        function P = get.P(obj)
            P = obj.Pa + obj.Pv;
        end
        
        function Pa = get.Pa(obj)
            Pa = obj.RHOa * refEQ.Ra * obj.Ta_K;
        end
        
        function Pv = get.Pv(obj)
            Pv = obj.RHOv * refEQ.Rv * obj.Ta_K;
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
        
        function ETo = get.ETo(obj)
            ETo = refEQ.ETo(obj);
        end
        
        function Rs = get.Rs(obj)
            Rs = obj.gLight.Q(time.HOUR)/1000;
        end
        
        function WW_Q = get.WW_Q(obj)
            WW_Q = 0.15*obj.no_people;
        end
        
        function stepCropConsumption(obj)
            if isempty(obj.currentSupply)
                obj.currentSupply = crop.determineYearlyCropDemand(obj.no_people, obj.cropBlock.type);
            else
                if obj.cropBlock.DVS > 2
                    % Add current harvest to available supply            
                    obj.currentSupply = [obj.currentSupply ...
                        obj.currentSupply(end) + obj.cropBlock.StorB/1000 - obj.currentHourlyDemand];
                else
                    obj.currentSupply = [obj.currentSupply ...
                        obj.currentSupply(end) - obj.currentHourlyDemand];
                end
            end
        end
        
        % These functions take too long to run!
        
%         function waterBlocks = get.waterBlocks(obj)
%             c = arrayfun(@(b) class(b.data), obj.blocks, 'UniformOutput', false);
%             mask = strcmp(c, 'water');
%             w = arrayfun(@(b) b.data, obj.blocks(mask));
%             waterBlocks = unique(w);
%         end
%         
%         function cropBlocks = get.cropBlocks(obj)
%             f = arrayfun(@(b) class(b.data), obj.blocks, 'UniformOutput', false);
%             mask = strcmp(f, 'crop');
%             c = arrayfun(@(b) b.data, obj.blocks(mask));
%             cropBlocks = unique(c);
%         end
        
        function setup(obj, Ac, Aw)
            % Create a lake object
%             w.x=1;
%             w.y=1;
%             w.w=3;
%             w.h=10;
%             w.d=1;
            
%             for i = w.x:1:(w.x+w.w-1)
%                 for j = w.y:1:(w.y+w.h-1)
%                 obj.blocks(i, j).data = wb;
%                 end
%             end
%             obj.waterBlocks = [obj.waterBlocks wb];
            
            % Create a crop object
%             c.x=4;
%             c.w=7;
%             c.y=1;
%             c.h=10;
%             
%             for i = c.x:1:(c.x+c.w-1)
%                 for j = c.y:1:(c.y+c.h-1)
%                     obj.blocks(i, j).data = wheatCrop;
%                 end
%             end
%             obj.cropBlocks = [obj.cropBlocks wheatCrop];
            
            % Create crop parameters
            obj.cropBlock = crop(Ac,cropType.Wheat,obj);
            obj.currentHourlyDemand = crop.determineYearlyCropDemand(obj.no_people,...
                obj.cropBlock.type) / (365.25*24);
            
            % Create a treatment train
            obj.waterBlock = waterBody(Aw, obj.hw, 20);
            obj.sedBlock = sedimentBasin(63e-6,obj.WW_Q,obj);
            obj.wetlandBlock = wetland1(obj.WW_Q,obj.sedBlock.TSS,1,obj);
            obj.filterBlock = filterUnit(obj.WW_Q,obj.wetlandBlock.TSS,obj);
            
            % Create lighting objects
            obj.gLight = lightObj(Ac, 640, 2, 1, obj, 6, 23, 'Growing light');
            obj.lLight = lightObj(obj.Area-obj.gLight.A, 12, 1000, 1000, obj, 6, 23, 'Living area light');
        
            % Create necessary handlers & variables
            obj.H = heat(obj);
            obj.W = water(obj);
        end
        
        function stepTime(obj)
            %% TODO: Update to change based off dt
            if obj.hour == 24
                obj.hour = 1;
                obj.day = obj.day + 1;
            else
                obj.hour = obj.hour + 1;
            end
            if isempty(obj.hours)
                obj.hours = 1;
            else
                obj.hours = [obj.hours max(obj.hours)+1];
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
            %% Calculate the amount of heat exchange
            obj.QQsn = 0;
            for j = 1:1:(3600/obj.dt)
                % Calculate amount of water released into air from
                % all water blocks
                obj.dMv = 0;
                w = obj.waterBlock;
                % Calculate the heat absorbed by water due to lighting
                % Valid for "daylight" hours (5AM - 8PM, for example)
                Qsn = obj.lLight.Q(obj.dt) * w.Aw;

                % Solve nonlinear equation for equilibrium temp. of water
                fun = @(T) obj.H.water(T, Qsn);
                [Ts_C, ~] = fzero(fun, w.Tw_C);

                % Calculate change in water temp based off surface
                % temp change
                ds = 0.01; % Surface depth = 1cm
                RHOs = refEQ.RHOw(Ts_C);
                Vs = w.Aw * ds;
                w.Tw_C = (w.RHOw*(w.Vw-Vs)*w.Tw_C + RHOs*Vs*Ts_C)/(w.RHOw*(w.Vw-Vs) + RHOs*Vs);
                dQsn = Qsn/1000;
                dQan = obj.H.Qan/1000;
                dQbr = obj.H.Qbr/1000;
                dQe = obj.H.Qe/1000;
                dQc = obj.H.Qc/1000;
                obj.dMv = obj.H.dMv;
                
                % Add short-wave radiation to accumulator
                obj.QQsn = obj.QQsn + dQsn;

                % Calculate heat exchange due to walls of BioSim
                qWall = (obj.kWall / obj.sWall) * (2*obj.Length*obj.Height + 2*obj.Width*obj.Height) * (obj.To_C - obj.Ta_C) * obj.dt;

                % Update temperature of air
                obj.Ta_C = (qWall+(dQbr+dQc+dQe-dQan))/(obj.W.Mvapor*obj.cv + obj.Ma*obj.ca) + obj.Ta_C;
                
                c = obj.cropBlock;
                if obj.dt == time.MINUTE
                    if mod(j,60)==0
                        c.step;
                        % Calculate necessary irrigation based off ET
                        obj.Vi = c.ET*1e-3; % convert mm/hr to m3/hr
                    end
                else
                    c.step;
                    % Calculate necessary irrigation based off ET
                    obj.Vi = c.ET*1e-3; % convert mm/hr to m3/hr
                end
                
                % Step water cycling
                obj.W.step;
            end
            
            % TODO: Calculate soil water conduct and EC values
            
            obj.trackTa         = [obj.trackTa obj.Ta_C];
            obj.trackTw         = [obj.trackTw w.Tw_C];
            obj.trackVcrop      = [obj.trackVcrop obj.W.Vcrop];
            obj.trackVstore     = [obj.trackVstore obj.W.Vstore];
            obj.trackVwetland   = [obj.trackVwetland obj.W.Vwetland];
            obj.trackVsed       = [obj.trackVsed obj.W.Vsed];
            obj.trackVreject    = [obj.trackVreject obj.W.Vreject];
            obj.trackVvapor     = [obj.trackVvapor obj.W.Vvapor];
            obj.trackDVS        = [obj.trackDVS obj.cropBlock.DVS];
            obj.trackLeafB      = [obj.trackLeafB obj.cropBlock.LeafB];
            obj.trackStemB      = [obj.trackStemB obj.cropBlock.StemB];
            obj.trackRootB      = [obj.trackRootB obj.cropBlock.RootB];
            obj.trackStorB      = [obj.trackStorB obj.cropBlock.StorB];
            obj.trackETo        = [obj.trackETo obj.ETo];
            obj.tracklLight     = [obj.tracklLight obj.lLight.EnergyUsage];
            obj.trackgLight     = [obj.trackgLight obj.gLight.EnergyUsage];
            obj.trackRO         = [obj.trackRO filterUnit.RO_P];
            obj.trackUV         = [obj.trackUV filterUnit.UV_P];
            output.Ta = obj.trackTa;
            output.Tw = obj.trackTw;
            
            obj.stepCropConsumption;
            obj.stepTime;
        end
    end
    
    methods (Static)
        function [At, Ac, Aw] = determineArea(no_people)
            At = model.determineTreatmentArea(no_people);
            Ac = model.determineCropArea(no_people);
            Aw = model.determineWaterStorageArea(no_people);
        end
        
        function At = determineTreatmentArea(no_people)
            % Determine wastewater flux
            WW_flux = 0.150 * no_people / 24; % m3/hr
            % Determine number of units
            no_UF = ceil(mod(WW_flux, filterUnit.UF_max_flux));
            no_RO = ceil(mod(WW_flux, filterUnit.RO_max_flux));
            no_UV = ceil(mod(WW_flux, filterUnit.UV_flux));
            A_UF = no_UF * filterUnit.UF_footprint;
            A_RO = no_RO * filterUnit.RO_footprint;
            A_UV = no_UV * filterUnit.UV_footprint;
            A_wetland = wetland1.determineArea(no_people);
            
            At = A_UF + A_RO + A_UV + A_wetland;
        end
        
        function Ac = determineCropArea(no_people)
            % Determine number of calories necessary
            calories = 2000 * 365.25 * no_people;
            % Determine area for number of calories
            caloriesPerM2 = 1000;
            Ac = calories / caloriesPerM2;
        end
        
        function Aw = determineWaterStorageArea(no_people)
            Aw = no_people * 20 / model.hw;
        end
    end
end

