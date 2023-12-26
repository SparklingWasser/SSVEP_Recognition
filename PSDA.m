function output = PSDA(Stim_freq, Data, fs, Nh)
%
% ==== Input Arguments ====
%  - Stim_freq : stimulation frequencies
%  - Data : channel x time serise
%  - fs : sampling rate
%  - Nh : number of harmonic frequencies to be employed for classification (default : 5)
%
% ==== Output Arguments ====
%  - output : classification result as an index for 'Stim_freq'
%
% Authors: Hodam Kim (rhg910907@hanyang.ac.kr), Seonghun Park (s.park7532@gmail.com), Jinuk Kwon (kowm2000@naver.com)
% Computational Neuroengineering (CoNE) Laboratory, Hanyang University, South Korea
%


if nargin == 3
    Nh = 5;
end

win_length = length(Data)/fs;
t = 1/fs:1/fs:win_length;
NFFT = length(t);

sig_length = length(Data);

Data = (hann(size(Data, 2))'.*Data);
Data = [Data zeros(size(Data, 1), win_length*fs-size(Data, 2))];

FFT_full = (2*abs(fft(Data, NFFT, 2))/sig_length).^2;

n_harmonic = 1:Nh;

for i = 1:length(Stim_freq)
    FFT(i, :, :) = FFT_full(:, round(win_length*Stim_freq(i)*n_harmonic+1));
end

P_reshape = reshape(FFT, [], 1);

[P_sort, P_idx] = sort(P_reshape, 'descend');

for i = 1:length(P_idx)
    idx = P_idx(i);
    if length(find(P_sort == P_sort(idx))) == 1
        break;
    end
end

[output, ~, ~] = ind2sub(size(FFT), idx);
