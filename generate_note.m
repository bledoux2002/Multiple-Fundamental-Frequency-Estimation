%MATLAB - Ben Ledoux
%generate_note function

function s = generate_note(freq, dur, sampFreq, N)
    t = 0 : 1/sampFreq : dur;
    s = sin(2 * pi * freq * t);
    for i = 2 : N
        %freq = freq * (freq + 1); %makes a really cool sound
        s = s + ((1 / i) * (sin(2 * pi * freq * t * i)));
    end
    s = s / (1.01 * max(max(s), -min(s)));
    %{
    fade = 0 : 8/sampFreq : 1;
    eighth = floor(sampFreq / 8);
    s(:, 1 : eighth + 1) = s(:, 1 : eighth + 1) .* fade;
    s(:, length(s) - eighth : length(s)) = s(:, length(s) - eighth : length(s)) .* flip(fade);
    %}
end