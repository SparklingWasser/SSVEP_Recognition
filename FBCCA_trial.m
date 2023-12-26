function output = FBCCA_trial(Stim_freq, Data, fs, a, b, subband)
%
% ==== Input Arguments ====
%  - Stim_freq : stimulation frequencies
%  - Data : channel x time serise
%  - fs : sampling rate
%  - Nh : number of harmonic frequencies to be employed for classification (default : 5)
%  - a, b, subband: hyper-parameters of FBCCA (default: a = 1.25, b = 0.5, subband = 10)
%
% ==== Output Arguments ====
%  - output : classification result as an index for 'Stim_freq'
%
% Authors: Hodam Kim (rhg910907@hanyang.ac.kr), Seonghun Park (s.park7532@gmail.com), Jinuk Kwon (kowm2000@naver.com)
% Computational Neuroengineering (CoNE) Laboratory, Hanyang University, South Korea
%
% Reference : Chen, Xiaogang, et al. "Filter bank canonical correlation analysis for implementing a high-speed SSVEP-based brain-computer interface." Journal of neural engineering 12.4 (2015): 046008.
%

if nargin == 3
    a = 1.25;
    b = 0.5;
    subband = 10;
end

min_freq = min(Stim_freq);
max_freq = min(Stim_freq) + (max(Stim_freq)-min(Stim_freq))*10;
Subband.freq_step = max(Stim_freq)-min(Stim_freq);

Data = Data';

for subband_idx = 1:subband
    
    Subband.min_freq = min_freq + Subband.freq_step*(subband_idx-1);
    
    [B, A] = butter(3, [Subband.min_freq-2 max_freq+2]/(fs/2));
    subband_sig = filtfilt(B, A, Data);
    
    CCA_sig = subband_sig(:, :)';
    [~, ~, raw_rho(subband_idx, :)] = CCA(Stim_freq, CCA_sig, fs);
    
    weight(subband_idx) = subband_idx^(-a)+b; 
end

for trial_all_freq = 1:length(Stim_freq)
    weighted_sum_rho(trial_all_freq) = squeeze(weight(1:subband_idx))*squeeze(raw_rho(1:subband, trial_all_freq)); % weight sum
end

[~, output] = max(weighted_sum_rho);















