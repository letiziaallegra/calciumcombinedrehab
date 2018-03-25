%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find number of sub-peaks  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [numPksDer_Force numPksDer_Fluo valuePksDer_Force locPksDer_Force valuePksDer_Fluo locPksDer_Fluo ] = ParFindNumberSubPeaks(StartForcePeak, StartFluoPeak, EndFluoPeak, der_IntervalForce_Sig_Long, der_IntervalFluoROI_Sig_Long)

%%%% find positive peaks in Force derivative
Int_Peak_Force    = [StartForcePeak:EndFluoPeak]';
[pksValue pksLoc] = findpeaks(der_IntervalForce_Sig_Long(Int_Peak_Force));
%
%%%%%%%%%
%%% Force peaks check
%%%%%%%%
prePoint = 5; preDiffValue = 0.1;
i_delPk = [];
for i=1:length(pksLoc)
    pk    = (pksValue(i));
    if (pksLoc(i)-prePoint)>0
        prePoint_Ok = prePoint;
    else
        prePoint_Ok = 1;
    end
    
    pkPre = (der_IntervalForce_Sig_Long(Int_Peak_Force(pksLoc(i)-prePoint_Ok)));
    if (pk-pkPre)<=preDiffValue
        i_delPk = [i_delPk; i];
    end
    
end
%%%%%%%%

pksValue(i_delPk) = [];
pksLoc(i_delPk) = [];
%
numPksDer_Force        = length(pksLoc);
valuePksDer_Force      = pksValue;
locPksDer_Force        = (Int_Peak_Force(pksLoc)); %[points]
% locPksDer_Force_sec    = t_long_sec(Int_Peak_Force(pksLoc)); %[sec]


%%%% find positive peaks in Fluo derivative
Int_Peak_Fluo          = [StartFluoPeak:EndFluoPeak]';
[pksValue pksLoc]      = findpeaks(der_IntervalFluoROI_Sig_Long(Int_Peak_Fluo));
numPksDer_Fluo         = length(pksLoc);
valuePksDer_Fluo       = pksValue;
% locPksDer_Fluo_sec     = t_long_sec(Int_Peak_Fluo(pksLoc)); %[sec]
locPksDer_Fluo         = (Int_Peak_Fluo(pksLoc)); %[points]

clear pksValue pksLoc