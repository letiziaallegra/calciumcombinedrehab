% [InfoPicchi] = PeaksFinder(SIG, Soglia, Fs, ch)
%Inputs: SIG   = signal
%       Soglia = Threshold
%       Fs     = Sampling Frequency
%       ch     = 0->vector point 1->time; 
%      %% Cor    = Correction -> it takes the start and the stop of the peaks (minima before and after the peak)
%Outputs:
%       DP: each column refers to a peak (the interval over threshold), with its amplitude values
%       InfoPicchi: column 1 -> peak start point
%                   column 2 -> peak duration
%                   column 3 -> max value peak
%                   column 4 -> max value peak point

%last updating 25 Feb 2014


function [InfoPicchi] = PeaksFinder_v3(SIG, Soglia,Fs,ch)

%NO
% %plot
% figure
% plot(SIG)
% hold on
% SIG_Abs = abs(SIG);
% plot(SIG_Abs,'r')
% 
% SIG_zScored = zscore(SIG);
% SIG = abs(SIG_zScored);
% Soglia=2;
%


%condizione iniziale
InfoPicchi = [];
num_picchi = 0;

%prendo solo i valori maggiori della soglia
% indici_SIG_soglia = find(SIG_Abs>Soglia);
indici_SIG_soglia = find(SIG>Soglia);

%vettore della durata dei picchi in cui (durata = tutti 1(o 2 o 3), quando incontra un altro valore (>3) è cominciato un altro picco) 
picchi = [1; diff(indici_SIG_soglia)];
  
    
if length(indici_SIG_soglia) >= 1 
    % +1 perchè devo considerare anche il primo picco altrimenti non
    %   conteggiato
    num_picchi = length(find(picchi>1))+1;
end


if(num_picchi ~= 0)
    
    f = find(picchi>1);
    f = [1;f];
    
    for i=1:num_picchi
        
       
            IstInPicco  = indici_SIG_soglia(f(i));
            if i~=num_picchi
                DurataPicco = round(f(i+1)-f(i));
            else
                DurataPicco = round(length(picchi)-f(i))+1;
            end
%             IntervalSig        = SIG_Abs( IstInPicco:IstInPicco+DurataPicco-1);
            IntervalSig        = SIG( IstInPicco:IstInPicco+DurataPicco-1);
            [AmpMax IstAmpMax] = max(IntervalSig);
        
        
            InfoPicchi = [InfoPicchi; [IstInPicco   DurataPicco  AmpMax(1)  IstInPicco+IstAmpMax(1)-1]];
        
    end
        
    if ch
        InfoPicchi(:,1) = InfoPicchi(:,1)/Fs;
        InfoPicchi(:,2) = InfoPicchi(:,2)/Fs;
        InfoPicchi(:,4) = InfoPicchi(:,4)/Fs;               
    end

   
end