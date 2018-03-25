% [ParForce ParFluo] = ParFromCorrespondingForceFluoPeaks(SigForce, SigFluo, StartForce, StartFluo, Dur, TrialsNum, Status, Fs)
% Compute info on one peak of force and its corresponding peak of fluo
% Inputs:
% Outputs:
%       DP: each column refers to a peak (the interval over threshold), with its amplitude values
%       InfoPicchi: column 1 -> trial
%                   column 2 -> status
%                   column 3 -> peak start point
%                   column 4 -> peak duration
%                   column 5 -> max value peak
%                   column 6 -> max value peak point (IstAmpMaxWSig)
%                   column 7 -> full width at half maximum (FWHM)
%                   column 8 -> area under the peak curve (AUC)
%                   column 9 -> peak-to-peak amplitude (PtPAmp)
%                   column 10 -> NumberOfRealPeak inside the peak     %non pi� il num of subpeaks inside the main peak (numPks)


%last updating 19 11 2015


function [ParForce ParFluo] = ParFromCorrespondingForceFluoPeaks_v2(SigForce, SigFluo, StartForcePeak, StartFluoPeak, EndForcePeak ,EndFluoPeak, numPksDer_Force, numPksDer_Fluo, stF_real, TimeWindow, dimOVS, trial_peak, status_peak, Fs, chTime)

for i=1:2
    
    %
    if i==1
        %Force
        IstInPicco      = StartForcePeak;                      %[points]
        IstEndPicco     = EndFluoPeak;                          %[points]
        SigLong         = SigForce;
        %num 
        numPk           = numPksDer_Force;  
    else
        %Fluo
        IstInPicco      = StartFluoPeak ;                     %[points]
        IstEndPicco     = EndFluoPeak;                          %[points] 
        SigLong         = SigFluo;
        numPk           = numPksDer_Fluo;
    end

    %%%%
    % trial when peak occurs
    TrialsNumPeak = trial_peak;
    %%%%
    
    %%%%
    % status when peak occurs
    statusPeak = status_peak;
    %%%%
    
    %%%%
    % start peak
    IstInPicco_InTheWholeTrial = IstInPicco + stF_real - TimeWindow(1);  %[points]
    %%%%
    
    %%%%
    IntervalTimePoint  = IstInPicco:IstEndPicco;    %[points]
    Sig                = SigLong(IntervalTimePoint);
    %%%%
    
    %%%%
    % peak duration &  NumberOfReaPeak 
    if i==1
        %per la forza considero la durata come la somma dei tempi dei
        %singoli picchi di forza considerati all'interno del picco di fluo
        threshForce     = SigLong(round(TimeWindow(1)+1)); %-> considero come soglia di forza il valore della forza all'istante di rilevazione del picco di fluo
        BinaryForce     = Sig>=threshForce;
        %numero picchi compresi nel picco di fluo
        NumberOfRealPeak    = sum((diff(BinaryForce)<0));
        %istante di fine picco 
        BinaryForceDiff = find(diff(BinaryForce)<0);
        
        if ~isempty(BinaryForceDiff)
            DurataPicco  = BinaryForceDiff(end)+1;         %[points]
        else
            %altrimenti uso fine picco fluo
            DurataPicco  = EndFluoPeak-IstInPicco+1;
        end
        
    else
        %numero picchi di fluo
        NumberOfRealPeak = 1;
        %peak duration
        DurataPicco  = IstEndPicco-IstInPicco+1;         %[points]
    end
    %%%%
    

    % max value peak;
    [AmpMax IstAmpMax] = max(Sig);
    % max value peak timepoint (reported for the whole signal)
    IstAmpMaxWSig = IntervalTimePoint(IstAmpMax);                              %[points]
    IstAmpMaxWSig_InTheWholeTrial = IstAmpMaxWSig + stF_real - TimeWindow(1);  %[points]
    %%%%
    
    %%%%
    % full width at half maximum (FWHM)
    %half maximum
    HM = AmpMax/2;
    %minimum difference
    %first point
    [HM_1 IstHM_1] = min(abs(Sig-HM));
    if IstAmpMax >= IstHM_1
        %second point
        [HM_2 IstHM_2] = min(abs(Sig(IstAmpMax:end)-HM));
    else
        %second point
        [HM_2 IstHM_2] = min(abs(Sig(1:IstAmpMax)-HM));
    end
    % FWHM
    FWHM = round(abs(IstHM_1-IstHM_2));
    %%%%
    
    %%%%
    % AUC
    AUC = trapz(IntervalTimePoint,Sig);
    %%%%
    
    %%%%
    % PtPAmp
    PtPAmp = abs(AmpMax-Sig(1));
    %%%%
    
    
    %store info
    PeaksInfo = [...
        TrialsNumPeak...
        statusPeak...
        IstInPicco_InTheWholeTrial...
        DurataPicco...
        AmpMax(1)...
        IstAmpMaxWSig_InTheWholeTrial...
        FWHM...
        AUC...
        PtPAmp...
        NumberOfRealPeak...
        ];
    %
    
    
    if chTime
        %ch == 1 -> time [sec]
        PeaksInfo(:,3) = PeaksInfo(:,3)/Fs;
        PeaksInfo(:,4) = PeaksInfo(:,4)/Fs;
        PeaksInfo(:,6) = PeaksInfo(:,6)/Fs;
        PeaksInfo(:,7) = PeaksInfo(:,7)/Fs;
    end
    
    
    
    if i==1
        ParForce = PeaksInfo;
    else
        ParFluo  = PeaksInfo;
    end
    
    
end




end


