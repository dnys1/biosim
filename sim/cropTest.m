load data/cropParams;

POOL = 0;           
LeafB = 10;         
MaxStemb = 6;
REPTIL = 0;
RootB = 5;
StemB = 6;
STEMP = 320;
StorB = 0;
VTIL = 250;
ddist = 0.005*MaxStemb;
FST = 0.05;
k = 0.6;            % coefficient of light extinction
Maxtil = 900;       % maximum number of tillers (tillers m^-2)
RAD = 17;           % MJ m^-2 d^-1
RRMAT = 0.3;
RUE = 0.7;          % Radiation use efficiency (g MJ^-1)
STW = 20;           % Dry biomass of one new tiller (g)
TBASE = 8;
TFLOW = 1500;
TMAT = 2000;
TMAX = 20;
TMIN = 19;
Totil = VTIL+REPTIL;
DVS = 0;
dt = 1;

SLA = [0 0.037; 1.0 0.018; 2.0 0.017];

i=0;
while DVS < 2
    i=i+1;
    CPPLt = interp1(CPPL(:,1),CPPL(:,2),DVS);
    CPPPt = interp1(CPPP(:,1),CPPP(:,2),DVS);
    CPRt = interp1(CPR(:,1),CPR(:,2),DVS);
    rrmortt = interp1(rrmort(:,1),rrmort(:,2),DVS);
    rrsent = interp1(rrsen(:,1),rrsen(:,2),DVS);
    DVEt = interp1(DVE(:,1),DVE(:,2),DVS);
    SLAt = interp1(SLA(:,1),SLA(:,2),DVS);
    
    CPL = CPPLt * (1-CPRt);
    CPP = CPPPt * (1-CPRt);
    CPS = (1-CPL-CPP)*(1-CPRt);
    
    PartL = CPL*POOL;
    PartS = CPS*POOL;
    PartSO= CPP*POOL;
    PartR = CPRt*POOL;
    PartLS = PartL + PartS;
    trackPartL(i)=PartL;
    trackPartS(i)=PartS;
    trackPartSO(i)=PartSO;
    trackPartR(i)=PartR;
    
    LAI = LeafB*SLAt;
    trackLAI(i) = LAI;
    
    RGrowth = RAD*RUE*(1-exp(-k*LAI));
    
    POOL = POOL + (RGrowth-PartS-PartL-PartSO-PartR)*dt;
    
    RSenL = rrsent*LeafB;
    LeafB = LeafB + (PartL-RSenL)*dt;
    
    rmaxstemb = PartS;
    MaxStemb = MaxStemb + rmaxstemb*dt;
    ddist = 0.005*MaxStemb;
    
    Rmrtv = rrmortt * VTIL;
    Rmat = RRMAT*VTIL;
    Rmortr = rrmortt*REPTIL;
    REPTIL = REPTIL + (Rmat-Rmortr)*dt;
    
    RootB = RootB + PartR*dt;
    
    if DVS > 1
        RTransloc = ddist;
    else
        RTransloc = 0;
    end
    StemB = StemB + (PartS-RTransloc)*dt;
    
    Dtemp = ((TMAX+TMIN)/2)-TBASE;
    STEMP = STEMP + Dtemp*dt;
    
    StorB = StorB + (PartSO+RTransloc)*dt;
    
    Rtil = PartLS*STW*(1-(VTIL/Maxtil))*DVEt;
    VTIL = VTIL + (Rtil-Rmat-Rmrtv)*dt;
    
    trackRootB(i)=RootB;
    trackStorB(i)=StorB;
    trackStemB(i)=StemB;
    trackLeafB(i)=LeafB;
    
    if STEMP < TFLOW
        DVS = STEMP / TFLOW;
    else
        DVS = 1+((STEMP-TFLOW)/(TMAT-TFLOW));
    end
end

i=1:1:i;
figure;
plot(i,trackRootB,i,trackStorB,i,trackStemB,i,trackLeafB);
hold on
legend('RootB','StorB','StemB','LeafB');
hold off