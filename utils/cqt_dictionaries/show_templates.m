%d=rdir('nmf_dic_*.mat');
d=rdir('nmf_dic_tenor_e.mat');

fig1 = figure;
for i=1:length(d)
    %
    [p,n,e]=fileparts(d(i).name);
    fprintf('Loading %s.: ' , d(i).name);
    load(d(i).name);
    fprintf('done. \n');
    dics_all = eval(n(5:end));
    
    for j=1:length(dics_all)
        dic = dics_all{j};        
        if ~isempty(dic) & sum(sum(dic.tpl))>0
            j,
            figure(fig1);         
            imagesc(dic.tpl);
            pause;
        end
    end
end
