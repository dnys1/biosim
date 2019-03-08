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
        f           % Estimated Green-Ampt infiltration rate (m/s)
        F           % Estimated Green-Ampt wetted front (m)
        psi_f       % Wetting front suction head (m)
        NH3_0       % Initial un-ionized ammonia (mg/L)
        NH4_0       % Initial ionized ammonia (mg/L)
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
        COD_0=600       % Initial COD in wastewater (mg/L)
        BOD_0=300       % Initial BOD5 in wastewater (mg/L)
        NH3T_0=25       % Initial NH3T (=[NH3]+[NH4+]) in wastewater (mg/L)
        TP_0=7          % Initial TP in wastewater (mg/L)
        pH_0=7.0        % Initial pH in wastewater
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
        
        function [COD, BOD, TSS, NH4, TP] = kinetics(obj)
            %KINETICS Outputs the concentrations of COD,BOD,TSS,NH4,TP
            %given initial concentrations which are assumed unchanging
            COD = obj.COD_0*exp(-obj.kV_COD*obj.t);
            BOD = obj.BOD_0*exp(-obj.kV_BOD*obj.t);
            TSS = obj.TSS_0*exp(-obj.kV_TSS*obj.t);
            NH4 = obj.NH4_0*exp(-obj.kV_NH4*obj.t);
            TP  = obj.TP_0*exp(-obj.kV_TP*obj.t);
        end
        
        function NH3_0 = get.NH3_0(obj)
        % Calculate the ammount of un-ionized ammonia using
        % Emmerson relationship
            pK = 0.09018 + 2729.92/obj.model.Ta_K;
            Fu = 1/(1+(10^-obj.pH_0/10^-pK));
            NH3_0 = obj.NH3T_0 * Fu;
        end
        
        function NH4_0 = get.NH4_0(obj)
        % Calculate the ammount of ionized ammonia using
        % Emmerson relationship
            NH4_0 = obj.NH3T_0 - obj.NH3_0;
        end
        
        function psi_f = get.psi_f(obj)
        % Estimate the wetting front suction head using
        % air-entry suction head
            psi_f = obj.psi_ae / 2;
        end
        
        function t = get.t(obj)
        % Calculation of hydraulic residence time t using
        % Green-Ampt modeling for a single-layer wetland
            if obj.q > obj.Ks
                h0 = obj.q - obj.Ks;
            else
                h0 = 0;
            end
            h = h0 - obj.psi_f;
            t = (obj.theta_s - obj.theta_r)/obj.Ks * (obj.z - h*log(1 + obj.z/h));
        end
        
        function f = get.f(obj)
        % Calculation of infiltration rate
            if obj.q < obj.Ks
                f = @(t)obj.q;
            else
                tp = obj.Ks*abs(obj.psi_f)*(obj.phi-obj.theta_r)/(obj.q*(obj.q-obj.Ks));
                T_ = abs(obj.psi_f)*(obj.phi-obj.theta_r)/obj.Ks;
                Fp = obj.q * tp;
                tc = Fp/obj.Ks - (abs(obj.psi_f)*(obj.phi-obj.theta_r))/obj.Ks * ...
                        log(1+Fp/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
                f = @(t)(obj.Ks*(0.707*(((t-tp+tc)+T_)/(t-tp+tc))^0.5 + 0.667 - 0.236*((t-tp+tc)/((t-tp+tc)+T_))^0.5 - ...
                        0.138 * ((t-tp+tc)/(t-tp+tc)+T_)));
            end
        end
        
        function Fstep = Fstep(obj, t)
            Fguess = obj.Ks*t;
            Fstep = obj.Ks*t + abs(obj.psi_f)*(obj.phi-obj.theta_r)*log(1+Fguess/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
            while abs(Fguess-Fstep) > 1e-5
                Fguess = obj.Ks*t + abs(obj.psi_f)*(obj.phi-obj.theta_r)*log(1+Fstep/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
                Fstep = obj.Ks*t + abs(obj.psi_f)*(obj.phi-obj.theta_r)*log(1+Fguess/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
            end
        end
        
        function F = get.F(obj)
        % Calculation of wetting front depth
            if obj.q < obj.Ks
                F = @(t)(obj.Fstep(t));
            else
                tp = obj.Ks*abs(obj.psi_f)*(obj.phi-obj.theta_r)/(obj.q*(obj.q-obj.Ks));
                T_ = abs(obj.psi_f)*(obj.phi-obj.theta_r)/obj.Ks;
                Fp = obj.q * tp;
                tc = Fp/obj.Ks - (abs(obj.psi_f)*(obj.phi-obj.theta_r))/obj.Ks * ...
                        log(1+Fp/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
                F = @(t)(obj.Ks*(0.529*(t-tp+tc) + 0.471*(T_*(t-tp+tc)+(t-tp+tc)^2)^0.5 + 0.0128*T_* ...
                        (log((t-tp+tc)+T_)-log(T_)) + 0.471*T_*(log((t-tp+tc)+0.5*T_+(T_*(t-tp+tc)+(t-tp+tc)^2)^0.5) - ...
                        log(0.5*T_))));
            end
        end
    end
end