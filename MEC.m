function output = MEC(Stim_freq, Data, fs, Nh)
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
% Reference : Friman, O., Volosyak, I., & Graser, A. (2007). Multiple channel detection of steady-state visual evoked potentials for brain-computer interfaces. IEEE transactions on biomedical engineering, 54(4), 742-750.

if nargin == 3
    Nh = 5;
end

%% Preprocessing

wo = 60/(fs/2);  bw = wo/35;
[b,a] = iirnotch(wo,bw);
Data = filtfilt(b,a,Data');

E = sum(Data.^2);
Data = (Data./sqrt(E))';

Ny = size(Data, 1);               % # of electrodes used

%% Calculating Needed Parameters 

Nf = length(Stim_freq);

Nt = size(Data, 2);

Y = Data';

%% Calculating X Matrix
t = linspace(0, Nt/fs, Nt);
X = zeros(Nf, Nt, 2*Nh);
for j = 1:Nf
    for i = 1:Nh
        X(j, :, 2*i-1) = sin(2*pi*i*Stim_freq(j)*t);
        X(j, :, 2*i) = cos(2*pi*i*Stim_freq(j)*t);
    end
end

% Calculating Ytild
Ytild = zeros(Nf, Nt, Ny);
for i = 1:Nf
    X_this_freq = squeeze(X(i,:,:));  
    Ytild(i,:,:) = Y - X_this_freq*inv(X_this_freq'*X_this_freq)*X_this_freq'*Y;
end

%% Calculating Ytild_trans*Ytild=Yeig and eigs and Weights

Yeig = zeros(Nf, Ny, Ny);
Ns = Ny;              % # of Artificial Channals

for i = 1:Nf
   Yeig(i,:,:) = (squeeze(Ytild(i, :, :)))'*squeeze(Ytild(i, :, :));
end

W = zeros(Nf, Ny, Ns);
for i = 1:Nf
    [V, D] = eig(squeeze(Yeig(i, :, :)));
    eigValues = abs(diag(D)');
    eigValues = sort(eigValues);
    
    
    denom = sum(eigValues);
    num = 0;
    for j = 1:length(eigValues)
        num = num + eigValues(j);
        if num/denom > 0.1 && j == length(eigValues)
            Ch_num(i) = j-1;
            break;
        elseif num/denom > 0.1
            Ch_num(i) = j;
            break;
        end
    end
    
    sqrt_eigValues = sqrt(eigValues);
    sorted_vectors = V(:, n);
    W(i, :, :) = sorted_vectors./repmat(sqrt_eigValues, Ny, 1) ;   
end

%% Calculating S 

S = zeros(Nf, Nt, Ns);
for i = 1:Nf
   W_this = squeeze(W(i, :, :));
   S(i,:,:) = Y*W_this;
end




%% Calculating P_hat & P & P_prim & Classifying

P_hat=zeros(1,Nf);
for i=1:Nf
    % for each freq
    S_this=squeeze(S(i,:,:));
    X_this=squeeze(X(i,:,:));
    for l=1:Ch_num(i)
       for k=1:Nh
          P_hat(i)=P_hat(i)+ norm(X_this(:,k)'*S_this(:,l))^2;
       end
    end
    
end

P_hat=P_hat/(Ns*Nh);
%P_hat Done 

P=P_hat/sum(P_hat);

[~, output] = max(P);




