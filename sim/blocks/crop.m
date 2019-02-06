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
        CPPL
        CPPP
        CPR
        rrmort
        rrsen
        DVE
        
        POOL        % 
        LeafB       % Biomass accumulated in leaf (kg)
        StorB       % Biomass accumulated in storage (kg)
        StemB       % Biomass accumulated in stems (kg)
        RootB       % Biomass accumulated in roots (kg)
        MaxStemB
        VTIL        % Vegetative tiller shoots
        REPTIL      % Reproductive tiller shoots
    end
    
    properties (Dependent)
        Kc          % Crop coefficient
        yield       % Yield if harvested now (kg)
        LAI         % Leaf area index
        Y           % Crop yield (kg)
        
        DVS         % Development stage 0 (emergence) - 1 (flowering) - 2 (maturity)
        GDD         % Growing degree days (oC-d)
        ETTot       % Total evapotranspiration until now (cm)
        ET          % Evapotranspiration at time t (now)
        RAD         % Radiation received for current timestep (MJ)
    end
    
    properties (Constant)
        % GENECROP Model parameters
        RUE = 1.2       % Radiation use efficiency (g MJ^-1)
        k = 0.6         % coefficient of light extinction
        STW = 20        % Dry biomass of one new tiller (g)
        Maxtil = 900    % Maximum number of tillers
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
                obj.rrmort = S.rrmort;
                obj.rrsen = S.rrsen;
                obj.DVE = S.DVE;
                
                % Initialize remaining model parameters
                obj.POOL = 0;
                obj.LeafB = 0;
                obj.StorB = 0;
                obj.StemB = 0;
                obj.RootB = 0;
                obj.MaxStemB = 6;
                obj.VTIL = 0;
                obj.REPTIL = 0;
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
        
        function LAI = get.LAI(obj)
           if obj.DVS == -1
               LAI = 0;
               return
           end
           SLA = [0 0.037; 1 0.018; 2 0.017];
           LAI = obj.LeafB * interp1(SLA(:,1),SLA(:,2),obj.DVS); 
        end
        
        function ET = get.ET(obj)
            ET = obj.Kc * obj.model.ETo;
        end
        
        function DVS = get.DVS(obj)
            % part of GENECROP model
            % https://www.apsnet.org/edcenter/advanced/topics/BotanicalEpidemiology/Pages/CropGrowthModeling.aspx
            if obj.GDD > obj.GDDemer
                if obj.GDD < obj.GDDflow
                    DVS = obj.GDD / obj.GDDflow;
                else
                    DVS = 1 + ((obj.GDD - obj.GDDflow)/(obj.GDDmat - obj.GDDflow));
                end
            else
                % Set to -1 if crop has not reached emergence
                DVS = -1;
            end
        end
        
        function step(obj)
            % Only continue if plant has emerged and not reached maturity
            if obj.DVS == -1 || obj.DVS > 2
                return
            end
            
            dt = 1/obj.model.dt; % Define separate dt var for calculations
            
            % Interpolate necessary variables from model parameters
            CPPLt = interp1(obj.CPPL(:,1),obj.CPPL(:,2),obj.DVS);
            CPPPt = interp1(obj.CPPP(:,1),obj.CPPP(:,2),obj.DVS);
            CPRt = interp1(obj.CPR(:,1),obj.CPR(:,2),obj.DVS);
            rrmortt = interp1(obj.rrmort(:,1),obj.rrmort(:,2),obj.DVS);
            rrsent = interp1(obj.rrsen(:,1),obj.rrsen(:,2),obj.DVS);
            DVEt = interp1(obj.DVE(:,1),obj.DVE(:,2),obj.DVS);
            
            % Calculate coefficients
            CPL = CPPLt * (1-CPRt);
            CPP = CPPPt * (1-CPRt);
            CPS = (1-CPL-CPP)*(1-CPRt);
            
            % Calculate partitioning of energy
            PartL   = CPL*obj.POOL;
            PartS   = CPS*obj.POOL;
            PartSO  = CPP*obj.POOL;
            PartR   = CPRt*obj.POOL;
            PartLS  = PartL+PartS;
            
            % Calculate the rate of growth based off ET
            RG = obj.RAD*obj.RUE*(1-exp(-obj.k*obj.LAI));
            
            % Caculate the amount of assimilates available for plant growth
            obj.POOL = obj.POOL + (RG-PartS-PartL-PartSO-PartR)*dt;
            
            % Caclulate production and death of tillers
            Rmrtv = rrmortt*obj.VTIL;       % Rate of mortality of vegetative tillers
            Rmortr = rrmortt*obj.REPTIL;    % Rate of mortality of reproductive tillers
            Rmat = rrmat*obj.VTIL;          % Number of tillers shifting from veg -> rep
            Rtil = PartLS*obj.STW*(1-(obj.VTIL/obj.Maxtil))*DVEt; % Rate of tiller production
            
            obj.REPTIL = obj.REPTIL + (Rmat-Rmortr)*dt;
            obj.VTIL = obj.VTIL + (Rtil-Rmat-Rmrtv)*dt;
            
            % Calculate accumulation of biomass in different parts
            % Stems
            rmaxstemb = PartLS;
            obj.MaxStemB = obj.MaxStemB + rmaxstemb*dt;
            ddist = 0.005*obj.MaxStemB;
            if obj.DVS > 1
                RTransloc = ddist;
            else
                RTransloc = 0;
            end
            obj.StemB = obj.StemB + (PartS-RTransloc)*dt;
            % Storage
            obj.StorB = obj.StorB + (PartSO+RTransloc)*dt;
            % Roots
            obj.RootB = obj.RootB + PartR*dt;
            % Leaves
            RSenL = rrsent*obj.LeafB;
            obj.LeafB = obj.LeafB + (PartL-RSenL)*dt;
        end
    end
end

