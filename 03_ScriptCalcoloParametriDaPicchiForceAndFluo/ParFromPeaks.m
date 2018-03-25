% [InfoPicchi] = PeaksFinder(SIG, Soglia, Fs, ch)
%Inputs: SIG   = filtered signal
%       Fs     = Sampling Frequency
%       StartStopPeak = Start and Stopo time points of the peaks
%       TrialsNum = num of the current trial
%       Status = Task Status
%       ch     = 0->vector point 1->time; 
%      %% Cor    = Correction -> it takes the start and the stop of the peaks (minima before and after the peak)
%Outputs:
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

%last updating 20 Oct 2015


function [PeaksInfo] = ParFromPeaks(SIG, Fs, StartStopPeak, TrialsNum, Status, ch)

%condizione iniziale
PeaksInfo = [];
num_picchi = size(StartStopPeak,1);

for i=1:num_picchi    
    
    IstInPicco   = StartStopPeak(i,1);
    IstEndPicco  = StartStopPeak(i,2);
    
    %%%%
    % trial when peak occurs
    TrialsNumPeak = TrialsNum(IstInPicco);
    %%%%
    
    %%%%
    % status when peak occurs
    statusPeak = Status(IstInPicco);
    %%%%
    
    %%%%
    % peak duration
    DurataPicco = IstEndPicco-IstInPicco+1;
    %%%%

    IntervalTimePoint  = IstInPicco:IstEndPicco;
    IntervalSig        = SIG( IntervalTimePoint);
    
    %%%%
    % max value peak;
    [AmpMax IstAmpMax] = max(IntervalSig);
    % max value peak timepoint (reported for the whole signal)
    IstAmpMaxWSig = IntervalTimePoint(IstAmpMax);
    %%%%
    
    %%%%
    % full width at half maximum (FWHM)
    %half maximum
    HM = AmpMax/2;
    %minimum difference
    %first point
    [HM_1 IstHM_1] = min(abs(IntervalSig-HM));
    if IstAmpMax >= IstHM_1
        %second point
        [HM_2 IstHM_2] = min(abs(IntervalSig(IstAmpMax:end)-HM));
    else
        %second point
        [HM_2 IstHM_2] = min(abs(IntervalSig(1:IstAmpMax)-HM));
    end
    % FWHM
    FWHM = round(abs(IstHM_1-IstHM_2));
    %%%%
    
    %%%%
    % AUC
    AUC = trapz(IntervalTimePoint,IntervalSig);
    %%%%
    
    %%%%
    % PtPAmp
    PtPAmp = abs(AmpMax-SIG(IstEndPicco));    
    %%%%
    
    %store info
    PeaksInfo = [PeaksInfo;...
        [...
        TrialsNumPeak...
        statusPeak...
        IstInPicco...
        DurataPicco...
        AmpMax(1)...
        IstAmpMaxWSig...
        FWHM...
        AUC...
        PtPAmp...
        ]];
    %
    
end

if ch
    %ch == 1 -> time
    PeaksInfo(:,3) = PeaksInfo(:,3)/Fs;
    PeaksInfo(:,4) = PeaksInfo(:,4)/Fs;
    PeaksInfo(:,6) = PeaksInfo(:,6)/Fs;
    PeaksInfo(:,7) = PeaksInfo(:,7)/Fs;
end

end


