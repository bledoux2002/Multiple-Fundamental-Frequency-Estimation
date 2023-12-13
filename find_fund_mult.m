% find_fund_mult.m - Multiple Fundamental Frequency Estimation by Ben Ledoux
% This program returns a vector of multiple detected fundamental frequencies.
% By analyzing peaks of amplitude in the signal window, and ignoring harmonic
% partials by omitting any freauencies determined to be 2x the frequency of
% a previous element (found starting from low frequencies to high)

function frequency_estimates = find_fund_mult(spectrum, sFreq, fft_size)
    F = (0:fft_size/2-1)*sFreq/fft_size;                                    % vector of frequencies in spectrum
    [m, i] = max(spectrum);                                                 % find max amp
    threshold = 0.1 * m;                                                    % set threshold to 10% of max amp
    fr_list = [];                                                           % initialize significant frequencies vector

    for i = 2 : length(spectrum) - 1                                        % go through spectrum of signal
        if (((spectrum(i) > spectrum(i - 1)) && (spectrum(i) > spectrum(i + 1))) && (spectrum(i) > threshold)) % if amp is peak and above threshold
            fr_list = [fr_list F(i)];                                       % append associated freq
        end
    end

    if (length(fr_list) > 0)                                                % if significant freqs are present in signal
        fr_list_filt = [fr_list(1)];                                        % init filtered vector to omit partials
    else                                                                    % if no significant freqs found in signal
        frequency_estimates = [];
        return;                                                             % return empty vector
    end

    qtone = 2 ^ (1 / 24);                                                   % used to find freqs within musical quartertone of expected partial (might be slightly off due to fft variations)
    for i = 2 : length(fr_list)                                             % check every significant peak found, assume first peak is a f0
        for j = 1 : i                                                       % compare to every significant peak leading up to it
            if ((fr_list(i) >= ((fr_list(j) * 2) / qtone)) && (fr_list(i) <= ((fr_list(j) * 2) * qtone))) % if peak is a partial of a previous peak
                                                                            % by check all sig freqs, can assume partial will be 2x one of them ,dont need to check if 4x, 8x, etc of established f0
                break;                                                      % skip over it
            elseif (j == i)                                                 % else if peak is not a partial of any lower peaks
                fr_list_filt = [fr_list_filt fr_list(i)];                   % concat to filtered list as new f0
            end
        end
    end
    frequency_estimates = fr_list_filt;
end