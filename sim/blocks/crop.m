classdef crop < handle
    %CROP Summary of this class goes here
    %   https://aces.nmsu.edu/aes/irrigation/wheat.html
    
    properties
        type        % Type of crop being grown (i.e. wheat)
        Ac          % Area of the crop (m2)
        tPlanted    % The time the crop was planted (used for tracking ET)
        Tbase       % Base temperature (C)
        Tcutoff     % Cutoff temperature (C)
        
        GDDemer     % GDD required to reach emergence
        GDDflow     % GDD required to reach flowering stage (C-day)
        GDDmat      % GDD required to reach maturity
        
        model       % A handle for the model
        
        % GENECROP Model variables
        CPPL        % Partitioning coefficient for leaves
        CPPP        % Partitioning coefficient for storage organs
        CPR         % Partitioning coefficient for roots
        rrsen       % Relative rate of leaf senescence
        
        POOL        % Total energy pool available from photosynthesis (g)
        LeafB       % Biomass accumulated in leaf (g)
        StorB       % Biomass accumulated in storage (g)
        StemB       % Biomass accumulated in stems (g)
        RootB       % Biomass accumulated in roots (g)
        MaxStemB    % Maximum stem biomass (g)
        
        trackLeafB  % Array of leafB vals
        trackStorB  % Array of storB vals
        trackStemB  % Array of stemB vals
        trackRootB  % Array of rootB vals
    end
    
    properties (Dependent)
        Kc          % Crop coefficient
        yield       % Yield if harvested now (kg)
        LAI         % Leaf area index
        Y           % Crop yield (kg)
        
        DVS         % Development stage 0 (emergence) - 1 (flowering) - 2 (maturity)
        GDD         % Growing degree days (oC-d)
        ET          % Evapotranspiration at time t (now)
        RAD         % Radiation received for current timestep (MJ)
        RUE         % Radiation use efficiency (randomized)
    end
    
    properties (Constant)
        % GENECROP Model parameters
        k = 0.6         % coefficient of light extinction
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
                obj.tPlanted = model.date;
                
                switch(obj.type)
                    case 'Wheat'
                        obj.Tbase = 0;
                        obj.Tcutoff = 30;
                        obj.GDDemer = 120;
                        obj.GDDflow = 1075;
                        obj.GDDmat = 1825;
                    case 'Corn'
                        obj.Tbase = 10;
                        obj.Tcutoff = 30;
                        %obj.GDDemer = 10 * rand + 55;
                        obj.GDDemer = 60;
                        obj.GDDflow = 1400;
                        obj.GDDmat = 2700;
                    case 'Soy'
                        obj.Tbase = 10;
                        obj.Tcutoff = 30;
                        % emer -http://www.coolbean.info/pdf/soybean_research/early_season/Predicting_soy_emergence.pdf
                        obj.GDDemer = 72;
                        % flow,mat =https://edisciplinas.usp.br/pluginfile.php/4111255/mod_resource/content/1/EC%233_2017_Souza2013.pdf
                        obj.GDDflow = 642;
                        obj.GDDmat = 1741;
                end
                
                % Load GENECROP model parameters
                S = load('data/cropParams');
                obj.CPPL = S.CPPL;
                obj.CPPP = S.CPPP;
                obj.CPR = S.CPR;
                obj.rrsen = S.rrsen;
                
                % Initialize remaining model parameters
                obj.POOL = 0;
                obj.LeafB = 0;
                obj.StorB = 0;
                obj.StemB = 0;
                obj.RootB = 0;
                obj.MaxStemB = 6;
            else
                error('Invalid crop type')
            end
        end
        
        function RAD = get.RAD(obj)
            %RAD Returns the amount of radiation in MJ for this timestep
            RAD = obj.model.gLight.Qperm2 * 60 * 1e-6 * obj.Ac;
        end
        
        function RUE = get.RUE(obj)
            RUE = refEQ.rand(0.7,1.2);
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
        
        function LAI = get.LAI(obj)
           if obj.DVS == -1
               LAI = 0;
               return
           end
           SLA = [0 0.037; 1 0.018; 2 0.017];
           LAI = obj.LeafB / obj.Ac * interp1(SLA(:,1),SLA(:,2),obj.DVS); 
        end
        
        function ET = get.ET(obj)
            ET = obj.Kc * obj.model.ETo;
        end
        
        function Y = get.Y(obj)
            Y = obj.StorB / 1000;
        end
        
        function DVS = get.DVS(obj)
            % part of GENECROP model
            % https://www.apsnet.org/edcenter/advanced/topics/BotanicalEpidemiology/Pages/CropGrowthModeling.aspx
            if obj.GDD > obj.GDDemer
                if obj.GDD < obj.GDDflow
                    DVS = obj.GDD / obj.GDDflow;
                else
                    DVS = 1 + (obj.GDD - obj.GDDflow)/(obj.GDDmat - obj.GDDflow);
                end
            else
                % Set to -1 if crop has not reached emergence
                DVS = -1;
            end
        end
        
        function step(obj)
            % Only continue if plant has emerged and not reached maturity
            if obj.DVS == -1 || obj.DVS > 2
                obj.trackLeafB = [obj.trackLeafB obj.LeafB];
                obj.trackStemB = [obj.trackStemB obj.StemB];
                obj.trackRootB = [obj.trackRootB obj.RootB];
                obj.trackStorB = [obj.trackStorB obj.StorB];
                return
            end
            
            % Some biomass needs to be established upon emergence
            % before GENECROP will work
            if obj.LeafB == 0
                obj.LeafB = 2*obj.Ac;
                obj.StemB = 1*obj.Ac;
                obj.RootB = 1*obj.Ac;
                obj.StorB = 0;
            end
            
            % Simulation run hourly whereas variables (except RG)
            % are based off daily parameters
            dt = 1/24;
            
            % Interpolate necessary variables from model parameters
            CPPLt = interp1(obj.CPPL(:,1),obj.CPPL(:,2),obj.DVS);
            CPPPt = interp1(obj.CPPP(:,1),obj.CPPP(:,2),obj.DVS);
            CPRt = interp1(obj.CPR(:,1),obj.CPR(:,2),obj.DVS);
            rrsent = interp1(obj.rrsen(:,1),obj.rrsen(:,2),obj.DVS);
            
            % Calculate coefficients of partitioning
            CPL = CPPLt * (1-CPRt);
            CPP = CPPPt * (1-CPRt);
            CPS = (1-CPL-CPP)*(1-CPRt);
            
            % Calculate partitioning of energy
            PartL   = CPL*obj.POOL;     % Rate of partitioning of assimilate towards leaves (g hr^-1)
            PartS   = CPS*obj.POOL;     % Rate of partitioning of assimilate towards stems (g hr^-1)
            PartSO  = CPP*obj.POOL;     % Rate of partitioning of assimilate towards storage organs (g hr^-1)
            PartR   = CPRt*obj.POOL;    % Rate of partitioning of assimilate towards roots (g hr^-1)
            
            % Calculate the rate of growth based off ET (g)
            RG = obj.RAD*obj.RUE*(1-exp(-obj.k*obj.LAI));
            
            % Caculate the amount of assimilates available (g)
            obj.POOL = obj.POOL + (RG-PartS-PartL-PartSO-PartR);
            
            % Calculate accumulation of biomass in different parts
            % Stems
            rmaxstemb = PartS;
            obj.MaxStemB = obj.MaxStemB + rmaxstemb;
            ddist = 0.005*obj.MaxStemB;
            % Translocation from stems to storage organs (g/day)
            if obj.DVS > 1
                RTransloc = ddist;
            else
                RTransloc = 0;
            end
            obj.StemB = obj.StemB + (PartS-RTransloc*dt);
            % Storage
            obj.StorB = obj.StorB + (PartSO+RTransloc*dt);
            % Roots
            obj.RootB = obj.RootB + PartR;
            % Leaves
            RSenL = rrsent*obj.LeafB;
            obj.LeafB = obj.LeafB + (PartL-RSenL*dt);
            
            obj.trackLeafB = [obj.trackLeafB obj.LeafB];
            obj.trackStemB = [obj.trackStemB obj.StemB];
            obj.trackRootB = [obj.trackRootB obj.RootB];
            obj.trackStorB = [obj.trackStorB obj.StorB];
        end
    end
end

