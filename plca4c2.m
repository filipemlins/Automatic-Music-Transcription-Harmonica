function [ww,ff,pp, bb, ppb, ppb20, xa, mask] = plca4c2(ww, x,  iter)
%
% Guitar transcription based on PLCA model
%
% Inputs:
%  ww    initial value of ww (basis 4-D tensor = dictionary template)
%  x     input distribution
%  iter  number of EM iterations [default = 35]
%
% Outputs:
%  ww   spectral bases
%  ff   pitch shifting
%  pp   pitch detection
%  pp5  pitch detection (20 cents resolution)
%  xa  approximation of input
%
% Original script from Emmanouil Benetos 2015
% (https://code.soundsoftware.ac.uk/projects/amt_plca_5d)
% Modificated by
% by Rodrigo Schramm 2017/2018

rng(666);

% // TODO // reducing matrix size to reduce computation (usualy, there is no
% relevant info below the 150th bin.
x = x(150:end,:);

% Get sizes
[M,N] = size(x);
sumx = sum(x);

pp=[];


K = size(ww,2); % pitch range
F = size(ww,3); % pitch shift
H = size(ww,4); % num spectrum templates per pitch and shift,

maskC=[];
maskBb=[];
maskG=[];

maskC = build_harmonica_constraints(96, N, 0); % K N B   Harmonica in C
maskBb = build_harmonica_constraints(96, N, -2); % K N B   Harmonica in Bb
maskG = build_harmonica_constraints(96, N, -5); % K N B   Harmonica in G

maskC = maskC(54:54+K-1,:,:); %%
maskBb = maskBb(54:54+K-1,:,:); %%
maskG = maskG(54:54+K-1,:,:); %%

mask = cat(3,maskC,maskBb,maskG);
%mask = maskG;
%mask = ones([size(maskBb,1), size(maskBb,2), 1]);

B = size(mask,3);

ww =  repmat(ww(150:end,:,:,:,1), [1 1 1 1 B]);

%sumx(sumx<0.005)=0;
%ss = diag(sumx>.005);

%ss = diag((sumx>5)); %% TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
%x = x*double(ss);


% Default training iterations
if ~exist( 'iter')
    iter = 30;
end

% Initialize
%% ----------------------------------------
if ~exist( 'ww') || isempty(ww)
    error('Please, define the input dictionary ww.');
end


%normalize spectral templates
% // TODO -- replace for loops below (this can be done previously offline)
% for k=1:K % num of pitches
%     for f=1:F  % num of freq bins (shift)
%         ww(:,k,f) = ww(:,k,f) ./ (sum(ww(:,k,f))+eps);
%     end
% end


%% ----------------------------------------
% pp = pitch activation
if ~exist( 'pp') || isempty(pp)
    pp = rand(K, N);
end
n=1:N;
pp(:,n) = repmat(sumx(n),K,1) .* (pp(:,n) ./ (repmat( sum( pp(:,n), 1), K, 1)+eps));

%% ----------------------------------------
% ff = 20 cent pitch shift
if ~exist( 'ff') || isempty( ff)
    ff = rand(F, K, N);
end
% normalize
f_resh = reshape(ff,[F K*N]);
f_resh = (f_resh ./ (repmat(sum(f_resh,1),[F 1])+eps));
ff = reshape(f_resh,[F K N]);

%% ----------------------------------------
% hh = template type (for many harmonica kinds)
if ~exist( 'hh') || isempty(hh)
    hh = rand(H, K, N);
end
% normalize
h_resh = reshape(hh,[H K*N]);
h_resh = (h_resh ./ (repmat(sum(h_resh,1),[H 1])+eps));
hh = reshape(h_resh,[H K N]);

%% ----------------------------------------
% bb = mask activation (harmonic contraints)
if ~exist( 'bb') || isempty(bb)
    bb = rand(B, N);
end
% normalize
b_resh = reshape(bb,[B N]);
b_resh = (b_resh ./ (repmat(sum(b_resh,1),[B 1])+eps));
bb = reshape(b_resh,[B N]);


%% ----------------------------------------
% Initialize update parameters
w_reshaped = reshape(ww,[M K*F*H*B]);  % freqBins x pitch_20cents
sumx = diag(sumx);


%% fixed masks (harmonic constraints)
mask_big = permute(repmat(mask, [1 1 1 F H]), [2 1 4 5 3]);% (t p f h b) <= K N B

fig1 = figure;
% Iterate
fprintf('iteration:   ');
for it = 1:iter
    %fprintf(' %d' , it);
    %% E-step
    
    ff_big = permute(repmat(ff, [1 1 1 H B]), [3 2 1 4 5]); % pitch shift (t p f h b)
    hh_big = permute(repmat(hh, [1 1 1 F B]), [3 2 4 1 5]); % template/harmonica type (t p f h b)
    bb_big = permute(repmat(bb, [1 1 F H K]), [2 5 3 4 1]); % harmonic contraint activ (t p f h b)
    pp_big = permute(repmat(pp, [1 1 F H B]), [2 1 3 4 5]); % pitch activ (t p f h b)    
    ffhhbbpp = ff_big.*hh_big.*(mask_big.*bb_big).*pp_big;    
    ffhhbbpp_reshaped = reshape(ffhhbbpp, [N K*F*H*B]); %% TODO
    %% posterior
    xa = w_reshaped * ffhhbbpp_reshaped'; % dot product estima o 'spectrogram input' como
    % combinacoes lineares dos "spectrum atoms"do dicionario
    %% scale probs based on the fitting (responsabilities)
    %--------------
    D = x ./ (xa+eps);
    
    %M-step
    WD = D' * w_reshaped;
    FP_big = ffhhbbpp_reshaped.* WD;
    %--------------
    
    %% marginalise variables
    FP = reshape(FP_big, [N K F H B]);
    ff = permute(squeeze(sum(sum(FP, 4),5)), [3 2 1]);
    hh = permute((sum(sum(FP, 3),5)), [4 2 1, 3]);
    %%d
    bb = permute(squeeze(sum(sum(sum(FP,2),3),4)), [2 1]);
    %%
    pp = permute(squeeze(sum(sum(sum(FP,3),4),5)), [2 1]);
    
    
    %% temporal constraints
    if 1
        bb = bb.^1.1;
    end
    
    if it>15
        pp = pp.^1.1;
        bb = bb.^5; % enforce sparsity across mask activations (force only 1 mask)
    end
    
    if 1
        bb = medfilt1(bb,20,[],2);
        pp = medfilt1(pp,10,[],2);
    end
    
    
    %% normalizations
    f_resh = reshape(ff,[F K*N]);
    f_resh = (f_resh ./ (repmat(sum(f_resh,1),[F 1])+eps));
    ff = reshape(f_resh,[F K N]);
    
    h_resh = reshape(hh,[H K*N]);
    h_resh = (h_resh ./ (repmat(sum(h_resh,1),[H 1])+eps));
    hh = reshape(h_resh,[H K N]);
    
    bb_resh = bb;
    bb_resh = (bb_resh ./ (repmat(sum(bb_resh,1),[B 1])+eps));
    bb = bb_resh;
    
    p_resh = pp;
    p_resh = (p_resh ./ (repmat(sum(p_resh,1),[K 1])+eps));
    pp = p_resh * sumx;
    
   
 
    %% ploting
    figure(fig1);
    dm1 = zeros([96, size(x,2)]);
    %dm2 = zeros([96, size(x,2)]);
    dm1(60:60+K-1,:) = pp;
    %dm2(60:60+K-1,:) = pp2*sumx;
    %subplot(2,1,1); im
    imagesc(dm1); axis xy; axis([1 size(dm1,2) 50 96]);
    %subplot(2,1,2); imagesc(dm2); axis xy;
    %subplot(2,1,2); imagesc(pp_aux); axis xy;
    title(it);shg; pause(.05);
    
     
    fprintf('\b\b%2d', it);  
end


ppb=zeros(size(pp));
[~,I]=max((bb),[],1);
for i=1:length(I)
    ppb(:,i)= pp(:,i).*mask(:,i,I(i));
end
ppb = medfilt1(ppb, 10,[],2);
%ppb = (ppb ./ (repmat(sum(ppb,1),[K 1])+eps));



%% generates the 20 cent resolution pitch activations
ffe = reshape(ff,[F*K size(pp,2)]);
ppb20 = zeros(size(ffe));
for j=1:K; ppb20(j*5-4:j*5,:) = repmat(pp(j,:),5,1); end
if F==5
    ppb20 = ffe.*ppb20;
end

close(fig1);

fprintf(' done.\n');

%ff = circshift(ff, [-1, 0]);
