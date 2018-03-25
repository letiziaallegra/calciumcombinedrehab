function sigfilt = butterHPfilt(sig,f,order,Fs,ripple)
% cheb2BPfilt band-pass filters the signal 'sig' between frequencies ['f_low' 'f_high']
% It uses a chebyshev type 2 filter of order 'order'
% 'Fs' is the sampling frequency of the signal 'sig'
if(nargin<5)
    ripple = 20;
end

Wn = f/(Fs/2); % Normalized stopband edge frequency
[B, A] = butter(order,Wn,'high'); 
sigfilt = filtfilt(B, A, sig);