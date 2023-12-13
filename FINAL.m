% Final Project: Multiple Fundamental Frequency Estimation by Ben Ledoux
% The purpose of this program is to find multiple fundamental frequencies
% in a sound file. Using my own find_fund_mult.m file, which is modelled
% after the partial difference model of estimating fundamental frequencies,
% along with concepts I learned in researching my paper, I have developed a
% convincing method of estimating and reproducing multiple fundamental
% frequencies in a MIDI file.

% Main function to loop through dataset and create midis for each file
function main
    % Test against collection of 3 notes in different consecutive
    % combinations (D, F, A, DF, FA, DA, DFA)
    estFundFreq("sound.wav");
%{
    % Testing against chords of 3 notes
    files = ["A#aug.wav" "A#dim.wav" "A#maj.wav" "A#min.wav" "Aaug.wav" "Adim.wav" "Amaj.wav" "Amin.wav" "Baug.wav" "Bdim.wav" "Bmaj.wav" "Bmin.wav" "C#aug.wav" "C#dim.wav" "C#maj.wav" "C#min.wav" "Caug.wav" "Cdim.wav" "Cmaj.wav" "Cmin.wav" "D#aug.wav" "D#dim.wav" "D#maj.wav" "D#min.wav" "Daug.wav" "Ddim.wav" "Dmaj.wav" "Dmin.wav" "Eaug.wav" "Edim.wav" "Emaj.wav" "Emin.wav" "F#aug.wav" "F#dim.wav" "F#maj.wav" "F#min.wav" "Faug.wav" "Fdim.wav" "Fmaj.wav" "Fmin.wav" "G#aug.wav" "G#dim.wav" "G#maj.wav" "G#min.wav" "Gaug.wav" "Gdim.wav" "Gmaj.wav" "Gmin.wav"];
    for i = 1 : length(files)
        estFundFreq(files(i)); % outputs MIDI file conversion
    end
%}
end

% Function to process audio file and return MIDI object for output
function estFundFreq(signal)

    fileName = signal;

% Read and analyze audio file
    [sound, fs] = audioread(signal);
    t = 0 : 1/fs : (length(sound) / fs);
    
    fft_size = 2^14;
    noverlap = floor(fft_size/1.1);
    [S, F, T] = spectrogram(sound, hanning(fft_size), noverlap, fft_size, fs);
    S = abs(S);
%{
    figure(1);
    imagesc( S.^0.5 );
    set(gca, 'YDir', 'normal' );
    xlabel('Time (Frames)');
    ylabel('Bins (Frequency)');
%}

% Create matrix of fundamental frequencies according, height denotes
% frames, width denotes total # of fundamental frequencies found
    FUND = [];                                                              % matrix holding organized f0's with 0's filling in empty spaces
    header = [];                                                            % vector of column hedaers identifying which column of FUND has which f0's (inherently if they exist or not)
    for i = 1 : width(S)
        est_list = [find_fund_mult(S(:, i), fs, fft_size)];                 % find fundamental freq for each frame
                                                                            % if f0 has not been reported yet, add it to header vector for later identification
        if (~isempty(header))                                               % first set of f0's don't need to be checked
            header = catHeader(header, est_list);                           % custom function for adding f0's to header
        else
            header = est_list;                                              % if first set of f0's, just set header to be same as list
        end

        if (length(header) > width(FUND))                                   % if more frequencies have been added to header
            FUND = [FUND zeros(height(FUND), (length(header) - width(FUND)))]; % add new column of zeros to account for new frequencies
        end
                                                                            % concat to FUND
        FUND = [FUND; zeros(1, length(header))];                            % new row for new frame
        for k = 1 : length(est_list)                                        % for every f0 estimated
            for l = 1 : length(header)                                      % check against every item in header
                if (est_list(k) == header(l))                               % locate column of f0
                    FUND(height(FUND), l) = est_list(k);                    % insert f0 to matrix at proper column
                end
            end
        end
    end

% Convert f0's to MIDI nums
    MIDI = zeros(size(FUND));                                               % matrix of MIDI nums corresponding to f0's of FUND matrix
    score = [];
    for i = 1 : width(FUND)                                                 % go thru FUND matrix
        for j = 1 : height(FUND)
            num = (round(12 * log2(FUND(j, i) / 27.5)) + 21);               % convert f0's to MIDI nums
            MIDI(j, i) = num;                                               % assign to corresponding spot in MIDI matrix
        end
    end

% MIDI Filter
    for h = 1 : width(MIDI)                                                 % go thru midi
        for i = 1 : height(MIDI)
            if MIDI(i, h) > 80                                              % if A5 or higher
                MIDI(i, h) = -Inf;                                          % remove
            end
        end
    end

% Combine consecutive notes into one note for efficiency
    window = ((length(sound) / fs) / height(FUND));                         % length of time in seconds for single MIDI note

    for i = 1 : width(MIDI)                                                 % go thru columns of MIDI nums
        MIDInote = [MIDI(1, i)];                                            % vector of MIDI note values without consecutive repeats
        MIDIstart = [0.0000];                                               % start time of each MIDInote
        MIDIdur = [window];                                                 % duration of each MIDInote
    
        for p = 2 : height(MIDI)                                            % go thru MIDI nums
            if (MIDI(p, i) == MIDI((p - 1), i))                             % if current MIDI num same as prev MIDI num
                MIDIdur(length(MIDIdur)) = MIDIdur(length(MIDIdur)) + window; % add index of time
            else                                                            % else current MIDI num NOT same as prev MIDI num
                MIDInote = [MIDInote MIDI(p, i)];                           % append next midi note
                MIDIstart = [MIDIstart (MIDIstart(length(MIDIstart)) + MIDIdur(length(MIDIdur)))]; % append prev midi start time + dur
                MIDIdur = [MIDIdur window];                                 % append new dur
            end
        end

        MIDIend = MIDIstart + MIDIdur;                                      % combine start times and durations for end times
        score = [score; MIDInote' MIDIstart' MIDIend'];                     % concat into matrix for MIDI, order doesn't matter due to definition of start/end times
    end

% Display MIDI data for debugging
%{
    for q = 1 : height(score)
        output = [num2str(score(q, 1)), ' ', char(9), num2str(score(q, 2)), ' s ', char(9), num2str(score(q, 3)), ' s']; %format output
        disp(output);
    end
%}

% Create MIDI file
    N = height(score);                                                      % length of file
    Notes = zeros(N,8);                                                     % vector of zeros to be replaced with MIDI data
    Notes(:, 1) = 1;                                                        % track 1
    Notes(:, 2) = 0;                                                        % channel 1
    Notes(:, 4) = 127;                                                      % FULL VOLUME
    for r = 1 : N                                                           % go thru all notes in matrix
        Notes(r, 3) = score(r, 1);                                          % note MIDI numbers
        Notes(r, 5) = score(r, 2);                                          % note on: start time
        Notes(r, 6) = score(r, 3);                                          % note off: start time + duration
    end

    NotesFilt = zeros(0, width(Notes));                                     % matrix to hold filtered MIDI matrix
    for i = 1 : height(Notes)                                               % go thru data of each note in MIDI
        if ((Notes(i, 3) ~= -Inf) && ((Notes(i, 6) - Notes(i, 5)) > 0.0708)) % if note is valid note and longer than one frame
            NotesFilt = [NotesFilt; Notes(i, :)];                           % concat to filtered matrix
        end
    end

    scoreMidi = matrix2midi(NotesFilt);                                         % convert to MIDI
%    scoreMidiInfo = midiInfo(scoreMidi,[], 1)                               %test to make sure matrix2midi() worked
    outName = convertStringsToChars(fileName);                              % create custom file names for each sound file, could use loop to omit .wav but not necessary
    outFile = [outName '.mid'];
    writemidi(scoreMidi, outFile);
    %https://midiplayer.ehubsoft.net/ to test midi file, or Windows Media Player

% unfinished attempt to combine consecutive notes into one longer note to
% prevent splicing issues (spikes in white noise between notes)
%{
    consecutive = [];
    for p = 1 : width(FUND)
        consecCol = []
        j = 1;
        while j <= height(FUND)
            k = 1;
            while ((FUND(j, p) == FUND(j + k, p)) && ((j + k)))
                k = k + 1;
            end
            consecCol = [consecCol; k];
            j = j + k;
        end

    end

    for p = 2 : height(FUND) %for each frame
        for m = 1 : width(FUND) %for each f0
            if (FUND(p, m) == FUND(p - 1, m))
                consecutive(p, m) = consecutive(p, m) + window; %$increment duration of consecutive note
            end
        end
    end
%}

% Using generate_note() to create the synthesized rendition of sound.wav.
% Works but has splicing issues cause spikes of white noise between notes,
% replaced with MIDI conversion.
%{
    window = ((length(sound) / fs) / height(FUND));                         % length of time in seconds to generate note
    EST = [];                                                               % initiate vector for final sound
    for k = 1 : height(FUND)                                                % for every row of FUND
        note = zeros(1, ceil(window * fs));                                 % prepare vector to hold combined generated notes
        noteCount = 0;                                                      % initialize counter for how many overlapping notes are present
        for m = 1 : width(FUND)                                             % for each element of a frame (row of FUND)
            if (FUND(k, m) == 0)                                            % if no f0 present
                continue;                                                   % skip
            else                                                            % if f0 present
                estNote = generate_note(FUND(k, m), window, fs, 1);         % generate note based on f0
                note = note + estNote;                                      % stack concurrent notes
                noteCount = noteCount + 1;                                  % increment stacked notes counter
            end
            if (noteCount > 0)                                              % if notes present 
                note = note / noteCount;                                    % average out amplitude
            end
        end
        EST = [EST note];                                                   % concat frame to final sound
    end


    t = 0 : 1/fs : (length(EST) / fs);
    fft_size = 2^14;
    noverlap = floor(fft_size/1.1);
    [S, F, T] = spectrogram(EST, hanning(fft_size), noverlap, fft_size, fs);
    S = abs(S);

    figure(2);
    imagesc( S.^0.5 );
    set(gca, 'YDir', 'normal' );  % make the y axis run bottom up
    xlabel('Time (Frames)');
    ylabel('Bins (Frequency)');

    audiowrite("sound_out.wav", EST, fs, "BitsPerSample", 16);
    figure(3);
    t = 0 : 1/fs : (length(EST) / fs);
    plot(t(2:end), EST);
%}
end

% Function to concatenate new frequencies into header. Not saving much on
% code complexity, just looks neater
function vec = catHeader(head, lis)
    vec = head;
    for q = 1 : length(lis)                                                 % for each item in list of reported f0's
        headBool = false;                                                   % bool if f0 exists in header vector already
        for r = 1 : length(vec)                                             % check against header vector
            if (lis(q) == vec(r))                                           % if it already exists
                headBool = true;                                            % set bool to true
                break;                                                      % skip rest of checks for that f0
            end
        end
        if (~headBool)                                                      % if f0 not in list
            vec = [vec lis(q)];                                             % concat f0 onto list
        end
    end
end