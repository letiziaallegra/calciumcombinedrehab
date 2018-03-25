function [ss f] = PowerSpectrumCompute_fun(signal, winDim, Fs)

%calcolo vettore potenza tramite finestre di Welch
[ss f] = pwelch(signal,winDim,winDim/2,[],Fs);

end

