classdef wetland1 < handle
    %WETLAND 1-layer wetland system
    %   Detailed explanation goes here
    
    properties
        Q           % Volumetric loading rate (m3/d)
        no          % Wetland number in series of wetlands
        TSS_0       % Influent TSS = sed basin effluent (mg/L array)

        model       % Main model reference
    end
    
    properties (Dependent)
        % Design parameters
        A           % Wetland area (m2)
        Vw          % Volume of water in the wetlands (m3)
        
        Ks          % Saturated hydraulic conductivity (m/s)
        t           % Hydraulic residence time (hr)
        tp          % Estimated Green-Ampt time to ponding (s)
        f           % Estimated Green-Ampt infiltration rate (m/s)
        F           % Estimated Green-Ampt wetted front (m)
        hL          % Clean head loss from TSS loading (m)
        psi_f       % Wetting front suction head (m)
        NH3_0       % Initial un-ionized ammonia (mg/L)
        NH4_0       % Initial ionized ammonia (mg/L)
        TSS         % Effluent TSS concentration (mg/L)
    end
    
    properties (Constant)
        z=0.8           % Wetland height (m)
        q=0.07/(3600*24)% Hydraulic loading rate (m/s)
        qd=0.07         % Hydraulic loading rate (m/d)
        phi=0.41        % Porosity
        theta_s=0.41    % Saturated water content
        theta_r=0       % Residual water content
        psi_ae=-1/44    % Air-entry suction head (m)
        lambda=0.18     % Fitted soil-water parameter
        kA_BOD=0.43     % Areal rate constant for BOD (m/day)
        kA_TSS=0.28     % Areal rate constant for TSS (m/day)
        kA_NH4=0.56     % Areal rate constant for NH4 (m/day)
        d10=5.5e-4      % d10 diameter (m)
        d60=3.1e-3      % d60 diameter (m)
        BOD_0=300       % Initial BOD5 in wastewater (mg/L)
        NH3T_0=25       % Initial NH3T (=[NH3]+[NH4+]) in wastewater (mg/L)
        pH_0=7.0        % Initial pH in wastewater
        kappaV=150      % Ergun coefficient for viscous loss
        kappaI=1.75     % Ergun coefficient for inertial loss
    
        NTIS=1          % Number of tanks in series (1=plug flow)
    end
    
    methods
        function obj = wetland1(Q, TSS_0, no, model)
            %WETLAND Construct an instance of this class
            %   Detailed explanation goes here
            obj.Q = Q;
            obj.TSS_0 = TSS_0;
            obj.no = no;
            obj.model = model;
        end
        
        function A = get.A(obj)
            A = obj.Q / obj.qd;
        end
        
        function Vw = get.Vw(obj)
            Vw = obj.A * obj.z * obj.phi;
        end
        
        function [BOD, TSS, NH4] = kinetics(obj)
            %KINETICS Outputs the concentrations of COD,BOD,TSS,NH4,TP
            %given initial concentrations which are assumed unchanging
            BOD = obj.BOD_0*(1+obj.kA_BOD*obj.z/(obj.NTIS*obj.qd))^-obj.no;
            TSS = obj.TSS_0*(1+obj.kA_TSS*obj.z/(obj.NTIS*obj.qd))^-obj.no;
            NH4 = obj.NH4_0*(1+obj.kA_NH4*obj.z/(obj.NTIS*obj.qd))^-obj.no;
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
        
        function Ks = get.Ks(obj)
        % Estimated hydraulic conductivity based on the void fraction and 
        % particle size
            Ks = 2.4622e-2*((obj.d10*10^3)^2*obj.phi^3/(1+obj.phi))^0.7825;
        end
        
        function t = get.t(obj)
        % Calculation of hydraulic residence time t using
        % Green-Ampt modeling for a single-layer wetland
            if obj.q > obj.Ks
                h0 = obj.q - obj.Ks - obj.hL;
            else
                h0 = -obj.hL;
            end
            h = h0 - obj.psi_f;
            t = (obj.theta_s - obj.theta_r)/obj.Ks * (obj.z - h*log(1 + obj.z/h));
        end

        function tp = get.tp(obj)
        % Calculation of the time to ponding using Green-Ampt
            tp = obj.Ks*abs(obj.psi_f)*(obj.phi-obj.theta_r)/(obj.q*(obj.q-obj.Ks));
        end
        
        function f = get.f(obj)
        % Calculation of infiltration rate
            % When q <= Ks, then the infiltration rate = loading rate
            if obj.q <= obj.Ks
                f = @(t)obj.q;
            % When q > Ks, then the infiltration rate can be estimated
            else
                T_ = abs(obj.psi_f)*(obj.phi-obj.theta_r)/obj.Ks;
                Fp = obj.q * obj.tp;
                tc = Fp/obj.Ks - (abs(obj.psi_f)*(obj.phi-obj.theta_r))/obj.Ks * ...
                        log(1+Fp/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
                f = @(t)(obj.Ks*(0.707*(((t-obj.tp+tc)+T_)/(t-obj.tp+tc))^0.5 + 0.667 - 0.236*((t-obj.tp+tc)/((t-obj.tp+tc)+T_))^0.5 - ...
                        0.138 * ((t-obj.tp+tc)/(t-obj.tp+tc)+T_)));
            end
        end
        
        function Fstep = Fstep(obj, t)
        % Guess and check method for Green-Ampt F(t)
            Fguess = obj.Ks*t;
            Fstep = obj.Ks*t + abs(obj.psi_f)*(obj.phi-obj.theta_r)*log(1+Fguess/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
            while abs(Fguess-Fstep) > 1e-5
                Fguess = obj.Ks*t + abs(obj.psi_f)*(obj.phi-obj.theta_r)*log(1+Fstep/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
                Fstep = obj.Ks*t + abs(obj.psi_f)*(obj.phi-obj.theta_r)*log(1+Fguess/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
            end
        end
        
        function F = get.F(obj)
        % Calculation of wetting front depth
            % If q <= Ks, then we use the guess and check method above
            if obj.q <= obj.Ks
                F = @(t)(obj.Fstep(t));
            % If q > Ks, then the direct calculation method can be applied
            else
                T_ = abs(obj.psi_f)*(obj.phi-obj.theta_r)/obj.Ks;
                Fp = obj.q * obj.tp;
                tc = Fp/obj.Ks - (abs(obj.psi_f)*(obj.phi-obj.theta_r))/obj.Ks * ...
                        log(1+Fp/(abs(obj.psi_f)*(obj.phi-obj.theta_r)));
                F = @(t)(obj.Ks*(0.529*(t-obj.tp+tc) + 0.471*(T_*(t-obj.tp+tc)+(t-obj.tp+tc)^2)^0.5 + 0.0128*T_* ...
                        (log((t-obj.tp+tc)+T_)-log(T_)) + 0.471*T_*(log((t-obj.tp+tc)+0.5*T_+(T_*(t-obj.tp+tc)+(t-obj.tp+tc)^2)^0.5) - ...
                        log(0.5*T_))));
            end
        end
        
        function hL = get.hL(obj)
        % Clean head loss in the system (Ergun equation)
            hL = obj.kappaV*(1-obj.phi)^2/obj.phi^3*refEQ.mu(obj.model.Ta_C)*obj.z*obj.q ...
                 /(refEQ.RHOw(obj.model.Ta_C)*9.81*obj.d60^2) + obj.kappaI*(1-obj.phi)/obj.phi^3 ...
                 *obj.z*obj.q^2/(9.81*obj.d60);
        end
        
        function TSS = get.TSS(obj)
        % Estimate the removal of TSS using rapid filtration equations
            y   = (1-obj.phi)^(1/3);
            kB  = 1.381e-23;
            Ha  = 1e-20;
            mu  = refEQ.mu(obj.model.Ta_C);
            RHOw= refEQ.RHOw(obj.model.Ta_C);
            T   = obj.model.Ta_K;
            dp  = obj.model.WW_dp;
            RHOp= obj.model.WW_RHOp;
            
            
            As  = 2*(1-y^5)/(2-3*y+3*y^5-2*y^6);
            Pe  = 3*pi*mu*dp*obj.d60*obj.q/(kB*T);
            NR  = dp./obj.d60;
            NV  = Ha/(kB*T);
            NG  = 9.81*(RHOp-RHOw)*dp.^2/(18*mu*obj.q);
            NA  = Ha./(3*pi*mu*dp.^2*obj.q);
            
            nuD=2.4*As^(1/3).*NR.^(-0.081).*NV.^(0.052).*Pe.^(-0.715);
            nuG=0.22*NR.^-0.24.*NV.^0.053.*NG.^1.11;
            nuI=0.55*As.*NA.^(1/8).*NR.^1.675;
            
            nu = nuD + nuG + nuI;
            a  = 1.0;
            k  = 3*(1-obj.phi)*nu*a/(2*obj.d60);
            
            TSS = sum(obj.TSS_0.*exp(-k*obj.z));
        end
    end
    
    methods (Static)
        function A = determineArea(no_people)
            A = 0.15 * no_people / wetland1.qd;
        end
    end
end