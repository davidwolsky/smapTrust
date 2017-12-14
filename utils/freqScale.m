function [freqScaled, freqUnit, scaleFact] = freqScale(freq)

 % Scale frequency
    if freq(1) > 1e12
        scaleFact = 1e12;
        freqUnit = 'THz';
    elseif freq(1) > 1e9
        scaleFact = 1e9;
        freqUnit = 'GHz';
    elseif freq(1) > 1e6
        scaleFact = 1e6;
        freqUnit = 'MHz';
    elseif freq(1) > 1e3
        scaleFact = 1e3;
        freqUnit = 'kHz';
    else
        scaleFact = 1;
        freqUnit = 'Hz';
    end
    
    freqScaled = freq./scaleFact;
end