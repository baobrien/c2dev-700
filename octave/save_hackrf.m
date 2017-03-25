
function save_hackrf(outfilename,samps)
    % Set up a buffer for output samples
    outsamps = zeros(1,length(samps)*2);
    plot(20*log10(abs(fft(samps(1:10000)))));
end
