

output_dir = 'output/';
d = dir([output_dir, '*.mat']); 
measures =[];

for j=1:length(d)

% XX = [];
% YY = [];
% fprintf('\nloading data (%d examples) for training: ', length(d));
% for i=1:length(d)
% %    if (i==j), continue; end; %% one leave, on out
%     fprintf(' %d', i)
%    fname = d(i).name;
%    load([output_dir, fname]);
%      
%    [a1,a2, fY, fX, fX2, fY_min, fY_max]=find_best_threshold(ppb, intCQT, gts);
%    %XX = [XX; [fX, fX2]];
%    %YY = [YY; fY];
%    XX{i} =  [fX, fX2];
%    YY{i} =  fY;
% end
%save featuresXXYY XX YY;

XXt = [];
YYt = [];
fprintf('\nloading data (%d examples) for training: ', length(d));
for i=1:length(d)
    if (i==j), continue; end; %% one leave, on out  
    XXt = [XXt; XX{i}];
    YYt = [YYt; YY{i}];    
end

fname = d(j).name;
load([output_dir, fname], 'plca_mat', 'gts', 'f', 'f2');
%
fprintf('\n[%d] training randomForest... against %s',j, fname);
tic;
%help templateTree
%randomForest = TreeBagger(30, XXt, YYt, 'oobpred','on', 'Method','regression','MinLeafSize', 20);
randomForest = TreeBagger(40, XXt, YYt, 'oobpred','on', 'Method','regression','NumVariablesToSample','all');
%Mdl = fitrtree([XXt,YYt],'MPG','PredictorSelection','curvature','Surrogate','on');

fprintf(' %.2fs  done.\n', toc);
plot(oobError(randomForest));
xlabel('number of grown trees');
ylabel('out-of-bag classification error');

lim = predict(randomForest, [f,f2]);
fprintf('File: %s\n' , fname);    
[Pre_s, Rec_s, F_s, pianoRoll] = fmeasure_dynamic(plca_mat, gts, lim);
fprintf('PLCA  -> Pre:%2.2f Rec:%2.2f F:%2.2f  (regression)\n', [Pre_s, Rec_s, F_s]);

measures(j,:,1) = [Pre_s, Rec_s, F_s];

th = 0.035; % best (time constant) option
[Pre_s, Rec_s, F_s] = compute_fmeasure(plca_mat, gts, th);
fprintf('PLCA  -> Pre:%2.2f Rec:%2.2f F:%2.2f   (th=%2.2f)\n', [Pre_s, Rec_s, F_s, th]);

measures(j,:,2) = [Pre_s, Rec_s, F_s];

end



