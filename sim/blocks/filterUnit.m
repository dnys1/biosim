classdef filterUnit
    %FILTER Filter object which houses UF, RO and UV systems
    
    properties
        Q       % Influent flow rate (m3/d)
        TSS_0   % Influent TSS conc. = wetlands effluent (mg/L array)
        
        model   % Reference to main model
    end
    
    properties (Dependent)
        P       % Total power consumption of treatment train (kW)
        TDS     % Effluent TDS levels (mg/L)
        UF_J    % Temperature corrected UF flux (L/m2/hr)
    end
    
    properties (Constant)
        % UF unit: https://www.waterfilters.net/pentek-u440-ultrafiltration-system-160380.html
        % Assumption: A 20 micron cartridge filter is used before this.
        UF_pore_size=0.065e-6       % UF pore size (m)
        UF_J0=21.41                 % UF measured flux at 25 deg C (L/m2/h)
        UF_footprint=0.093;         % UF footprint (m2)
        UF_min_flux=0;              % Min UF flux for single unit (m3/hr)
        UF_max_flux=1.416;          % Max UF flux for single unit (m3/hr)
        
        % RO unit: https://www.espwaterproducts.com/content/FLEXEON_AT-500_1000_Spec%20Sheet.pdf
        RO_recovery=0.41            % RO Recovery rate (%)
        RO_salt_rejection=0.985     % RO salt rejection (%)
        RO_perm_flow=2.61           % RO permeate flow (lpm)
        RO_min_feed_flow=6.4        % RO minimum feed flow (lpm)
        RO_max_feed_flow=15.14      % RO maximum feed flow (lpm)
        RO_min_conc_flow=3.78       % RO minimum conc flow (lpm)
        RO_min_flux=0.384;          % RO min feed flux (m3/hr)
        RO_max_flux=0.9084;         % RO max feed flux (m3/hr)
        RO_P=0.746                  % RO power consumption (kW)
        RO_footprint=0.18;          % RO footprint (m2)
        
        % UV unit: https://viqua.com/wp-content/uploads/LIT-520326_TapFamily_SpecSheet-17-HR.pdf
        UV_dose=40                  % UV dose (mJ/cm2)
        UV_flux=0.2                 % UV flux (m3/hr)
        UV_P=0.013                  % UV power consumption (kW)
        UV_footprint=0.093;         % UV footprint (m2)
    end
    
    methods
        function obj = filterUnit(Q, TSS_0, model)
            %FILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Q = Q;
            obj.TSS_0 = TSS_0;
            obj.model = model;
        end
        
        function P = get.P(obj)
            
        end
        
        function UF_J = get.UF_J(obj)
            J_corr=[25;1];
            UF_J = obj.UF_J0;
        end
    end
end

