
addpath('utils');
addpath('cqt');
addpath('constant_q_toolbox');
addpath('constant_q_toolbox/VQT');
addpath('constant_q_toolbox/CQT');
addpath('metrics');
addpath('ACA-Code-master');
addpath('other_amt_methods');
addpath('other_amt_methods/Vincent_multipitch');
addpath('other_amt_methods/Pertusa_multipich');
addpath('other_amt_methods/Essentia_multipitch');

% load plca 
load  dictionary/dict_harm1f5d5_r.mat;

root_dir = [pwd, '/'];

%% configure dataset path
%path_dir = 'datasets/chords/';
%path_dir = 'datasets/teste001/';
%path_dir = 'datasets/Annotaded/*/*/';
path_dir = 'datasets/Annotaded/Harmonica_em_C/Intervalos/';
%path_dir = 'datasets/Annotaded/Harmonica_em_Bb/Hucklebuckintro2/';
% path_dir = 'datasets/Annotaded/Harmonica_em_G/*/';
output_dir = 'output/';
filesw = rdir([path_dir,'*.wav']);

%% avaliar isso
pitch_act_threshold = 0.01;

measures=[];

nfiles = length(filesw);
for k=1:nfiles
    
    filename = filesw(k).name;
    
    [fpath, filename ,fext] = fileparts(filename);
    input_wav = [fpath,'/', filename ,fext];
    input_gt = [fpath,'', filename, 'Hz.txt'];
    
    fprintf('Processing %s: ', filename);
    
    %% load input wav file
    [y,Fs] = audioread(input_wav);
    y=y(:,1);
    
    fmin = 27.5;
    [intCQT, Xcq] = computeVQT(input_wav);
    intCQT(1:190,:)=0;
    intCQT= intCQT(:,1:5:end); %% reducing the time resolution (this accelerates the PLCA)
    
    %% spectrogram factorisation    
    iter = 20;
    [ww,ff,pp, bb, ppb,ppb20, xa, mask] = plca4c2(sd5d(:,:,:,:,:), intCQT, iter);
    
    % //TODO
    %ppb20 = circshift(ppb20, [-2,0]);
    
    % padd pitch range at begining (because I removed it from the dictionary)
    pact_s = zeros([96,size(intCQT,2)]);
    pact_s(54:54+size(ppb,1)-1,:) = ppb; % first pitch in the dictionary is 61
    
    pact_20 = zeros(size(intCQT));
    pact_20(54:54+size(ppb20,1)-1,:) = ppb20; % first pitch in the dictionary is 61
    
    
    %% load/creates the ground-truth data
    if ~exist(input_gt,'file'); return; end
    
    notes_gt = dlmread(input_gt);
    [gt,gt2,gts] = create_gt_matrix(notes_gt, pact_20, Fs, length(y));
    gts = gts(1:96,1:size(intCQT,2));
    gt = gt(1:size(intCQT,1),1:size(intCQT,2));
    
    pitch_act_threshold = 0.04;
    pact_s_bin = pact_s > pitch_act_threshold;
    pact_20_bin = pact_20 > pitch_act_threshold;
    
    %% pitch refinement
    %[pact_20_bin,pact_s_bin] = pitch_refinement(intCQT, pact_20, pitch_act_threshold);
    
    % %% extracts metrics
    [Pre_s, Rec_s, F_s] = compute_fmeasure(pact_s_bin, gts, pitch_act_threshold);
    [Pre, Rec, F] = compute_fmeasure(pact_20_bin, gt, pitch_act_threshold);
    fprintf('PLCA(s)    -> Pre:%2.2f Rec:%2.2f F:%2.2f    th=%2.3f\n', [Pre_s, Rec_s, F_s, pitch_act_threshold]);
    fprintf('PLCA(20)   -> Pre:%2.2f Rec:%2.2f F:%2.2f    th=%2.3f\n', [Pre, Rec, F, pitch_act_threshold]);
    measures(k).plca = [Pre_s, Rec_s, F_s, pitch_act_threshold];
    measures(k).plca20 = [Pre, Rec, F, pitch_act_threshold];
    
    
    figure; imagesc(gts(1:end,:) + 5*double(pact_s>pitch_act_threshold)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar; title('PLCA');
    
    %save([output_dir, filename, '_G'], 'filename', 'pact_s', 'pact_20', 'ppb', 'bb' , 'mask','pp', 'xa', 'intCQT');
    
    
%     
    if 0
        
        [h, notes] = mf0_pertusa([root_dir, input_wav]);
        [~, ~, pitch_pertusa] = create_gt_matrix(notes, gts, Fs, length(y));
        pitch_pertusa = pitch_pertusa(1:size(gts,1),1:size(gts,2));
        [Pre_p, Rec_p, F_p] = compute_fmeasure(pitch_pertusa, gts, 0.01);
        fprintf('Pertusa -> Pre:%2.2f Rec:%2.2f F:%2.2f \n', [Pre_p, Rec_p, F_p]);
                figure; imagesc(gts(1:end,:) + 5*double(pitch_pertusa(1:end,:)>.0)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar;title('pertusa');
        measures(k).pertusa = [Pre_p, Rec_p, F_p];
        

        % bock
        bock_dir = fullfile(root_dir, 'BockExperimentFiles/BockResults/');
        bock_filename = fullfile(bock_dir, [filename,'.out']);
        notes = dlmread(bock_filename);  
        [~, ~, pitch_bock] = create_gt_matrix(notes, gts, Fs, length(y));
        pitch_bock = pitch_bock(1:size(gts,1),1:size(gts,2));
        [Pre_b, Rec_b, F_b] = compute_fmeasure(pitch_bock, gts, 0.01);
        fprintf('Bock -> Pre:%2.2f Rec:%2.2f F:%2.2f \n', [Pre_b, Rec_b, F_b]);
                figure; imagesc(gts(1:end,:) + 5*double(pitch_bock(1:end,:)>.0)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar;title('bock');
        measures(k).bock = [Pre_b, Rec_b, F_b];   
        
        %magenta
        magenta_dir = fullfile(root_dir, 'MagentaExperiment/MagentaResults/');
        magenta_filename = fullfile(bock_dir, [filename,'.out']);
        notes_magenta = dlmread(magenta_filename);  
        [~, ~, pitch_magenta] = create_gt_matrix(notes_magenta, gts, Fs, length(y));
        pitch_magenta = pitch_magenta(1:size(gts,1),1:size(gts,2));
        [Pre_b, Rec_b, F_b] = compute_fmeasure(pitch_magenta, gts, 0.01);
        fprintf('Magenta -> Pre:%2.2f Rec:%2.2f F:%2.2f \n', [Pre_b, Rec_b, F_b]);
                figure; imagesc(gts(1:end,:) + 5*double(pitch_magenta(1:end,:)>.0)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar;title('magenta');
        measures(k).bock = [Pre_b, Rec_b, F_b];         
        
        
        % run vincent method
        [~,filename,~] = fileparts([root_dir, input_wav]);
        mf0file = [root_dir,'other_amt_methods/Vincent_multipitch/vincent_mf0_out/' , filename, '.txt'];
        multipitch_tracking(input_wav,mf0file);
        notes = dlmread(mf0file);
        [~, ~, pitch_vincent] = create_gt_matrix(notes, gts, Fs, length(y));
        pitch_vincent=pitch_vincent(1:96,:);
        [Pre_s, Rec_s, F_s] = compute_fmeasure(pitch_vincent, gts, 0.01);
        fprintf('Vincent -> Pre:%2.2f Rec:%2.2f F:%2.2f \n', [Pre_s, Rec_s, F_s]);
        figure; imagesc(gts(1:end,:) + 5*double(pitch_vincent(1:end,:)>.0)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar;title('vincent');
        measures(k).vincent = [Pre_s, Rec_s, F_s];
        
        % run klapuri and melodia (essentia)
        mf0file = [root_dir,'other_amt_methods/Essentia_multipitch/essentia_mf0_out/' , filename];
        [h, notes_klapuri, notes_melodia] = mf0_essentia([root_dir, input_wav],mf0file);
        
        m_klapuri = convert_essentia2matrix(gts,length(y)/Fs, notes_klapuri);
        m_melodia = convert_essentia2matrix(gts,length(y)/Fs, notes_melodia);
        
        [Pre_s, Rec_s, F_s] = compute_fmeasure(m_klapuri, gts, 0.01);
        fprintf('Klapuri -> Pre:%2.2f Rec:%2.2f F:%2.2f \n', [Pre_s, Rec_s, F_s]);
        figure; imagesc(gts(1:end,:) + 5*double(m_klapuri(1:end,:)>.0)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar;title('klapuri');
        measures(k).klapuri = [Pre_s, Rec_s, F_s];
        
        [Pre_s, Rec_s, F_s] = compute_fmeasure(m_melodia, gts, 0.01);
        fprintf('Melodia -> Pre:%2.2f Rec:%2.2f F:%2.2f \n', [Pre_s, Rec_s, F_s]);
        figure; imagesc(gts(1:end,:) + 5*double(m_melodia(1:end,:)>.0)); axis xy ; axis([1,size(gts,2), 55, 96]); colorbar;title('melodia');
        measures(k).melodia = [Pre_s, Rec_s, F_s];
        
       close all;
     end
end

save([output_dir, 'all_files.mat'], 'measures');
