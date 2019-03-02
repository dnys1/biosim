classdef wetland1 < handle
    %WETLAND 1-layer wetland system
    %   Detailed explanation goes here
    
    properties
        A           % Wetland area (m2)
        z           % Wetland height (m)
        q           % Hydraulic loading rate (m/s)
        model       % Main model reference
    end
    
    properties (Dependent)
        t           % Estimated Green-Ampt HRT (s)
        psi_f       % Wetting front suction head (m)
    end
    
    properties (Constant)
        phi=0.41        % Porosity
        theta_s=0.41    % Saturated water content
        theta_r=0       % Residual water content
        psi_ae=-1/44    % Air-entry suction head (m)
        lambda=0.18     % Fitted soil-water parameter
        Ks=0.2          % Saturated hydraulic conductivity (m/s)
        kV_COD=2.64/24  % Volumetric rate constant for COD (hr^-1)
        kV_BOD=3.68/24  % Volumetric rate constant for BOD (hr^-1)
        kV_TSS=2.59/24  % Volumetric rate constant for TSS (hr^-1)
        kV_NH4=0.66/24  % Volumetric rate constant for NH4 (hr^-1)
        kV_TP =0.40/24  % Volumetric rate constant for TP  (hr^-1)
        kA_COD=0.20/24  % Areal rate constant for COD (m/hr)
        kA_BOD=0.27/24  % Areal rate constant for BOD (m/hr)
        kA_TSS=0.19/24  % Areal rate constant for TSS (m/hr)
        kA_NH4=0.05/24  % Areal rate constant for NH4 (m/hr)
        kA_TP =0.03/24  % Areal rate constant for TP  (m/hr)
        d10=6e-3        % d10 diameter (m)
        d60=8.5e-3      % d60 diameter (m)
    end
    
    methods
        function obj = wetland1(A, z, q, model)
            %WETLAND Construct an instance of this class
            %   Detailed explanation goes here
            obj.A = A;
            obj.z = z;
            obj.q = q;
            obj.model = model;
        end
        
        function [COD, BOD, TSS, NH4, TP] = kinetics(obj, C)
            %KINETICS Outputs the concentrations of COD,BOD,TSS,NH4,TP
            %given input concentrations. C must be a vector of length 5
            if len(C) ~= 5
                error('C must be vector of length 5 [COD, BOD, TSS, NH4, TP]');
            end
            COD = C(1)*exp(-obj.kV_COD*obj.t);
            BOD = C(2)*exp(-obj.kV_BOD*obj.t);
            TSS = C(3)*exp(-obj.kV_TSS*obj.t);
            NH4 = C(4)*exp(-obj.kV_NH4*obj.t);
            TP  = C(5)*exp(-obj.kV_TP*obj.t);
        end
        
        function psi_f = get.psi_f(obj)
            psi_f = obj.psi_ae / 2;
        end
        
        function t = get.t(obj)
            if obj.q > obj.Ks
                h0 = obj.q - obj.Ks;
            else
                h0 = 0;
            end
            h = h0 - obj.psi_f;
            t = (obj.theta_s - obj.theta_r)/obj.Ks * (obj.z - h*log(1 + obj.z/h));
        end
    end
end

