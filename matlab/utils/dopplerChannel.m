%SUPPORT FUNCTION
function y = dopplerChannel(x,fs,chanParams)
    % Form an output vector y comprising paths of x with different
    % delays, Dopplers, and complex gains
    numPaths = length(chanParams.pathDelays);
    maxPathDelay = max(chanParams.pathDelays);
    txOutSize = length(x);
    
    y = zeros(txOutSize+maxPathDelay,1);
   
    for k = 1:numPaths
        pathOut = zeros(txOutSize+maxPathDelay,1);

        % Doppler
        pathShift = frequencyOffset(x,fs,chanParams.pathDopplerFreqs(k));
    
        % Delay and gain
        pathOut(1+chanParams.pathDelays(k):chanParams.pathDelays(k)+txOutSize) = ...
            pathShift * chanParams.pathGains(k);
            
        y = y + pathOut;
    end
end