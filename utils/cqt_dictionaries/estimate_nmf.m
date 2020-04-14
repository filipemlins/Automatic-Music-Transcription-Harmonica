

d=rdir('dic_*.mat');
%d=rdir('dic_tenor_a.mat');

for i=1:length(d)
    %
    [p,n,e]=fileparts(d(i).name);
    fprintf('Processing %s.: ' , d(i).name);
    load(d(i).name);
    dics_all = eval(n);
    
    
    %
    parfor j=1:length(dics_all)
        dic = dics_all{j};        
        if isempty(dic) 
            continue;
        end        
        tpl =  zeros([size(dic.m,1), 88]);
        tpl3 = zeros([size(dic.m,1), 88, 3]);
        tpl5 = zeros([size(dic.m,1), 88, 5]);
        
        %% robust fitt of pitch evolution (used to avoid octave jumps)
        yy = dic.pitch;
        xx = 1:length(yy);
        brob = robustfit(xx,yy);
                    
        for k=1:length(dic.dict_cell);
            %k,
            
            if (isnan(dic.pitch(k))) || (abs(brob(1)+brob(2)*k - yy(k))> 2) % 2 semitones
                continue;
            end
           a =  dic.dict_cell{k};
           % why elements of a can be below zero? Maybe because of the
           %    interpolation during the subresolution shifting.
           a(a<0)=0;
           % skip small dicts
           if size(a,2)<5
               continue;
           end
           if min(sum(a,2)) == 0
               a = a+0.00001;
           end
           
           %% TODO verify the settings of VQT/CQT ?
           % pitch location
           idx = round((dic.pitch(k)-4)/5);
           
           %[W,H,errs,vout] = nmf_kl_sparse_v(a(:,2:end-1), 5, 'alpha', .01);
           [W,H,errs,vout] = nmf_kl_sparse_v(a, 1, 'alpha', .01);
           dic.nmf{k} = W;
           tpl(:,idx) = W;
           
           [W,H,errs,vout] = nmf_kl_sparse_v(a, 3, 'alpha', .01);
           dic.nmf3{k} = W;
           tpl3(:,idx,:) = reshape(W, [size(dic.m,1), 1, 3]);
           
           [W,H,errs,vout] = nmf_kl_sparse_v(a, 5, 'alpha', .01);
           dic.nmf5{k} = W;
           tpl5(:,idx,:) = reshape(W, [size(dic.m,1), 1, 5]);
        end     
        dic.tpl = tpl;
        dic.tpl3 = tpl3;
        dic.tpl5 = tpl5;
        %
        dics_all{j} = dic;
                
    end
    assignin('caller', n, dics_all);
    save(['nmf_',n], n);
    fprintf('.  [%s]  saved.\n' , ['nmf_',n]);
    clear(n);
    
end    



