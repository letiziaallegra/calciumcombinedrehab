function [ signal_filt ] = NotchFilter(signal, Fs, freq)

    %Design
    f_notch = freq;
    wo = f_notch/(Fs/2);
    bw = wo/35;
    [b,a] = iirnotch(wo,bw); 
    %Application
    signal_filt = filtfilt(b, a, signal);

end

