function output = MSI(Stim_freq, Data, fs, Nh, L)
% 
% ==== Input Arguments ====
%  - Stim_freq : stimulation frequencies
%  - Data : channel x time serise
%  - fs : sampling rate
%  - Nh : number of harmonic frequencies to be employed for classification (default : 5)
%  - L : data length in seconds
%
% ==== Output Arguments ====
%  - output : classification result as an index for 'Stim_freq'
%
% Authors: Hodam Kim (rhg910907@hanyang.ac.kr), Seonghun Park (s.park7532@gmail.com), Jinuk Kwon (kowm2000@naver.com)
% Computational Neuroengineering (CoNE) Laboratory, Hanyang University, South Korea
%
% Reference : Zhang, Yangsong, et al. "Multivariate synchronization index for frequency recognition of SSVEP-based brain-computer interface." Journal of neuroscience methods 221 (2014): 32-40.


if nargin == 3 || nargin == 5
elseif nargin > 5
    error('At most 5 input arguments are acceptable');
elseif nargin < 3
    error('At least 3 input arguments are needed.');
end

if nargin == 3
    Nh = 5;
    L = size(Data, 2)/fs;
elseif nargin == 4
    L = size(Data, 2)/fs;
end

Ref_time = pi * Stim_freq' * linspace(0 , L, fs*L);

Ref_signal = zeros(fs*L,2*length(Stim_freq)*Nh);
for i = 1 : Nh
    Ref_signal(:,2*length(Stim_freq)*(i-1)+1:2*length(Stim_freq)*i) = [sin(2*i*Ref_time); cos(2*i*Ref_time)]';
end

Reference = permute(reshape(Ref_signal, fs*L, length(Stim_freq), 2*Nh), [3, 1, 2]);

S = zeros(1,length(Stim_freq));
X = Data(:,1:fs*L)';
X = detrend(X)';
X = 2*(X-min(X,[],2))./(max(X,[],2)-min(X,[],2))-1; % normalization
Nc = size(X,1); % Number of channels of X
for this_freq = 1:length(Stim_freq)
    Y = squeeze(Reference(:,:,this_freq));
    C = cov([X',Y']);
    U = zeros(size(C));
    U(1:Nc,1:Nc) = C(1:Nc,1:Nc)^(1/2);
    U(Nc+1:end,Nc+1:end) = C(Nc+1:end,Nc+1:end)^(1/2);
    R = U\C/U';
    e = eig(R);
    e = e/trace(R);
    S(this_freq) = 1+e'*log(e)/log(trace(R));
end

[~, output]  = max(S);

