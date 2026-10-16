% MUMT 307 Final Project: Phase Vocoder Implementation in MATLAB
% Owen MacDonald

function y = phaseVocoder(x, tStretch, pitchShift, winLen, overlap)

    N = length(x);
    
    % Initialize parameters
    pFactor = 2^(pitchShift/12); % Pitch factor
    stretch = tStretch * pFactor; % Final stretch

    H = round(winLen * (1-overlap)); % Sample hop length
    strH = round(H * stretch); % Stretched hop length
    win = hann(winLen, 'periodic'); % Periodic Hanning window for spectral analysis
    
    numFrames = floor((N - winLen) / H) + 1; % Number of frames for STFT
    
    y = zeros(numFrames * strH + winLen, 1); % Initialize output signal
    
    % Initialize phase accumulator to first frame
    frame1 = x(1:winLen) .* win;
    X1 = fft(frame1);
    X1 = X1(1:winLen/2+1);
    prevPhase = angle(X1);
    phaseAcc = prevPhase;
    
    expPhase = 2 * pi * (0:winLen/2)' * H / winLen; % Calculate expected phase advance
    
    % Spectral analysis loop
    for i = 0:numFrames-1
    
        frame = x(i*H+1 : i*H+winLen) .* win; % Frame to analyze
        X = fft(frame);
        X = X(1:winLen/2+1); % Keep positive frequencies only
    
        % Convert complex FFT data to magnitude and phase
        mag = abs(X);
        phase = angle(X);
    
        % Calculate phase difference
        phaseDiff = phase - prevPhase;
        prevPhase = phase; % Update previous phase
    
        % Remove expected phase advance
        phaseDiff = phaseDiff - expPhase;
    
        % Reduce to [-pi, pi] range
        phaseDiff = phaseDiff - 2 * pi * round(phaseDiff / (2*pi));
    
        % Update phase accumulator
        phaseAcc = phaseAcc + (phaseDiff + expPhase) * (strH / H);
        
        % Resynthesize time scaled output
        Y = mag .* exp(1i * phaseAcc); % Convert back to complex form for resynthesis
        Y = [Y; conj(Y(end-1:-1:2))]; % Mirror negative frequencies

        frameOutput = real(ifft(Y)) .* win;
        y(i*strH+1 : i*strH+winLen) = y(i*strH+1 : i*strH+winLen) + frameOutput;
    
    end
    
    % Resample at 1/pFactor times the original sample rate
    [p, q] = rat(1/pFactor); % Returns two integers p and q to approximate 1/pFactor
    y = resample(y, p, q);
    
    % Pad or trim to exact desired length
    len = round(N * tStretch);
    
    if length(y) > len
        y = y(1:len);
    else
        y = [y; zeros(len - length(y), 1)];
    end

    % Scale to [-1, 1] range
    y = y / max(abs(y));

end