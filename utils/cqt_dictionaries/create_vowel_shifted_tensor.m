
%function [shiftedW,singers] = create_vowel_shifted_tensor()

% Loads sound state templates and creates a 5-D tensor:
% pitch shifting x pitch x sound state x frequency x instrument

% Load note templates
%W = zeros(5,size(W,1),size(W,2),size(W,3),size(W,4));

%path_wav = '/Users/schramm/Desktop/vowel_templates/';
path_wav = './';


voice_types = {'soprano', 'alto' , 'tenor', 'baritone', 'bass'};
%voice_types = {'soprano', 'alto' , 'tenor'};
vowels = {'a', 'e' , 'i', 'o', 'u'};
%vowels = {'a', 'e'};

W=[];
singers = {};
for vv = 1:length(voice_types)
    for oo = 1:length(vowels)
        
        file_name = [path_wav, 'nmf*_dic_', voice_types{vv}, '_', vowels{oo} ,'.mat'],
        
        d=rdir(file_name);        
        min_samples = Inf;
        
        
        for k=1:length(d)
            fname = d(k).name;
            fprintf('[%d] loading dic: [%s].\n', k,fname);
            load(fname);
            [~,singer] = fileparts(fname);
            sg = eval(singer(6:end)); %% TODO
            l = length(sg);
            min_samples = min(l,min_samples);
            singers{k,1} = singer;
            singers{k,2} = l;
            
            for i=1:l;
                %[k,i,],
                X = sg{i}.tpl;
                [Y,~] = fill_gaps(X); %figure; imagesc(Y)
                %[p,s,w]
                %imagesc(Y'); shg; pause(.01);
                W(:,:,oo,vv,i) = Y';
            end
        end
    end
end

        % [M R K F O V]
        P = size(W,1);
        M = size(W,2);
        O = size(W,3);
        V = size(W,4);        
        S = size(W,5);
        F = 5;
        
        %shiftedW = zeros([M,S,P,F,O,V]);
        shiftedW = zeros([F,P,M,O,V,S]);
                
        for i=1:S
        for o=1:O
        for v=1:V        
            shiftedW(1,:,:,o,v,i) =  cat(2, W(:,3:545,o,v,i), zeros(88,2));
            shiftedW(2,:,:,o,v,i) =  cat(2, W(:,2:545,o,v,i), zeros(88,1));
            shiftedW(3,:,:,o,v,i) =  W(:,:,o,v,i);
            shiftedW(4,:,:,o,v,i) =  cat(2, zeros(88,1), W(:,1:544,o,v,i));
            shiftedW(5,:,:,o,v,i) =  cat(2, zeros(88,2), W(:,1:543,o,v,i));
        end
        end
        end;
        
        
        save('shiftedW','shiftedW', 'singers');
        fprintf('5-D tensor [shiftedW] saved.\n');
        if 0
            for i=1:25;
                subplot(5,5,i);
                imagesc(squeeze(shiftedW(3,:,2,:,i))');
                voz = strrep(singers{i,1}, '_', '-');
                xlabel([sprintf('%d - ', i), voz(10:end)]);
            end;
        end
        
