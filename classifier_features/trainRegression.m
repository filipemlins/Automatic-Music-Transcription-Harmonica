
addpath('utils');
addpath('cqt');
addpath('metrics');

load  dict_harm1f5d5_r.mat;

root_dir = [pwd, '/'];

%% configure dataset path 
%path_dir = 'datasets/chords/';
%path_dir = 'datasets/teste01/';
path_dir = 'datasets/Annotaded/*/*/';
path_dir = 'datasets/Annotaded/Harmonica_em_C/*/';

output_dir = 'output/';
% 
% %filesw = rdir([path_dir,'cancer*.wav']); 
% filesw_l{1} = rdir([path_dir,'heyjudeG05*.wav']); 
% %filesw = rdir([path_dir,'*C2posGpart1*.wav']); 
% filesw_l{2} = rdir([path_dir,'*rout*.wav']); 
% 
% %filesw = dir([path_dir,'untit*.wav']); 
% filesw_l{3} = rdir([path_dir,'chor*.wav']); 
% 
% %filesw = dir([path_dir,'*04.wav']); 
% %filesw = dir([path_dir,'*76.wav']); 
% %filesw = dir([path_dir,'*20.wav']); 
% 
% filesw_l{4} = rdir([path_dir,'heyjudeG03*.wav']); 
% filesw_l{5} = rdir([path_dir,'heyjudeG01*.wav']); 
% filesw_l{6} = rdir([path_dir,'byebyebirdpart1.wav']); 
% 
% filesw_l{7} = rdir([path_dir,'cancerianosemlarsolo1.wav']); 
% filesw_l{8} = rdir([path_dir,'OhSusanapart2GaitaC.wav']); 
% filesw_l{9} = rdir([path_dir,'intervalos.wav']); 
% 

filesw = rdir([path_dir,'*.wav']); 

XX=[];
YY=[];
XX2=[];

% for vv = 1:length(filesw_l);
%     filesw = filesw_l{vv};

nfiles = length(filesw);
for k=1:nfiles;
filename = filesw(k).name;

[fpath, filename ,fext] = fileparts(filename);
input_wav = [fpath,'/', filename ,fext];
input_gt = [fpath,'/', filename, '.txt'];
 
fprintf('Processing %s: ', filename);


%% load input wav file
[y,Fs] = audioread(input_wav);
y=y(:,1);

%% compute CQT
% maxFreq = Fs/2;
% bins = 12*bins_div;  % bins per octave
% sparKernel= sparseKernel(minFreq, maxFreq, bins, Fs);
% intCQT = schramm_cqt(y,sparKernel);
% intCQT = circshift(intCQT,[-2,0]);

fmin = 27.5;
%[intCQT, Xcq] = computeCQT(input_wav);
[intCQT, Xcq] = computeVQT(input_wav, 30);
%intCQT(1:150,:)=0;    
%intCQT = circshift(intCQT, [1,0]);

%% spectrogram factorisation
iter = 15;
[ww,ff,pp,pp5,pp2, bb, ppb, xa, mask] = plca4c(sd5d(:,:,:,:,:), intCQT, iter);

%pp5 = circshift(pp5, [-2,0]);

% padd pitch range at begining (because I removed it from the dictionary)
ppe = zeros(size(intCQT));
ppe(54:54+size(pp5,1)-1,:) = pp5; %first pitch in the dictionary is 55 (54 is for noise)

notes_gt = dlmread(input_gt);
[gt,gt2,gts] = create_gt_matrix(notes_gt, ppe, Fs, length(y));
gts = gts(1:96,1:size(intCQT,2));
%

if 1
[a1,a2, fY,fX, fX2, fYmin, fYmax] = find_best_threshold(ppb, intCQT, gts);
%[a1,a2, fY, fX, ]=find_best_threshold(ppb, intCQT, gts);

XX = [XX; fX];
XX2 = [XX2; fX2];
YY = [YY; fY];
end

if 0
     plca_mat = zeros(size(gts));
    plca_mat(54:54+size(ppb,1)-1,:)= ppb;

[cp,ci,cy]=comparePitchActivations(gts,plca_mat, mask,bb);

XX = [XX; cp];
YY = [YY; cy];

end
end



% end

%treeB = TreeBagger(100,XX,YY,'oobpred','on', 'Method','regression');
treeB = TreeBagger(100,XX,YY,'oobpred','on', 'Method','classification');

save TreeB_r1 treeB XX YY XX2;
