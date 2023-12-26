function output = EMSI(Stim_freq, Data, fs, Nh, sub_hf, L)
%
% ==== Input Arguments ====
%  - Stim_freq : stimulation frequencies
%  - Data : channel x time serise
%  - fs : sampling rate
%  - Nh : number of harmonic frequencies to be employed for classification (default : 5)
%  - sub_hf : whether to use sub-harmonic frequency component (default: 0)
%  - L : data length in seconds
%
% ==== Output Arguments ====
%  - output : classification result as an index for 'Stim_freq'
%
% Authors: Hodam Kim (rhg910907@hanyang.ac.kr), Seonghun Park (s.park7532@gmail.com), Jinuk Kwon (kowm2000@naver.com)
% Computational Neuroengineering (CoNE) Laboratory, Hanyang University, South Korea
%
% Reference : Zhang, Yangsong, et al. "Multivariate synchronization index for frequency recognition of SSVEP-based brain-computer interface." Journal of neuroscience methods 221 (2014): 32-40.
%

if nargin == 3 || nargin == 6
elseif nargin > 6
    error('At most 6 input arguments are acceptable');
elseif nargin < 3
    error('At least 3 input arguments are needed.');
end

if nargin == 3
    Nh = 5;
    sub_hf = false;
    L = size(Data, 2)/fs;
elseif nargin == 4
    sub_hf = false;
    L = size(Data, 2)/fs;
elseif nargin == 5
    L = size(Data, 2)/fs;
end

Ref_time = pi * Stim_freq' * linspace(0 , L, fs*L);

Ref_signal = zeros(fs*L, 2*length(Stim_freq)*Nh);
for i = 1 : Nh
    Ref_signal(:, 2*length(Stim_freq)*(i-1)+1:2*length(Stim_freq)*i) = [sin(2*i*Ref_time); cos(2*i*Ref_time)]';
end

Reference = permute(reshape(Ref_signal, fs*L, length(Stim_freq), 2*Nh), [3, 1, 2]);    

if sub_hf
    Reference(end+1, :, :) = sin(Ref_time)';
    Reference(end+1, :, :) = cos(Ref_time)';
end

S = zeros(1,length(Stim_freq));
X = Data(:,1:fs*L);
X = [X;circshift(X,1,2)];
X = (X-mean2(X))/std2(X);
for this_freq = 1:length(Stim_freq)
    Y = squeeze(Reference(:,:,this_freq));
    Y = (Y-mean2(Y))/std2(Y);    
    C = cov(X',Y');
    U = C.^(-1/2);
    U(1,2) = 0; U(2,1) = 0;
    R = U*C*U';
    e = eig(R);
    e = e/trace(R);
    S(this_freq) = 1+e'*log(e)/log(trace(R));
end

[~, output]  = max(S);


