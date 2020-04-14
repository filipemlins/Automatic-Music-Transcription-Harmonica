function [h, notes] = mf0_pertusa(input_wav)
h = 0;
notes=[];

[path_wav,wav_name,ext] = fileparts(input_wav);

%mf0_path = '/Users/schramm/Desktop/Pertusa_multipich/';
%mf0_path = '/Users/schramm/GoogleDrive/CARREIRA/UFRGS/RESEARCH/PROJECTS_UFRGS/AMT/automaticharmonicatranscription/other_amt_methods/Pertusa_multipich/';
mf0_path = 'other_amt_methods/Pertusa_multipich/';
%mf0_tool = [mf0_path, 'pertusa_exe.sh']; 
mf0_tool = ['./pertusa.sh']; 


curdir = pwd; 
cd(mf0_path);

dtemp = tempname;
[ss,mm,mi] = mkdir(dtemp);
%dir(input_wav),
%[mf0_tool, ' ', input_wav, ' ', [dtemp, '/'], ' ',  wav_name ],
[status,result] = system([mf0_tool, ' ', input_wav, ' ', [dtemp, '/'], ' ',  wav_name ]);

%status,
%fprintf(['Moving: ', [dtemp, '/*.csv'], '   to  ',  [path_wav, '/'], '\n']);
%dir([dtemp]),
d=dir([dtemp, '/*.csv']);
if length(d)>0
    %movefile([dtemp, '/*.csv'],  [path_wav, '/pertusa_mf0_mirex/']);    
    movefile([dtemp, '/*.csv'],  ['./pertusa_mf0_mirex/']);    
    notes = dlmread(['./pertusa_mf0_mirex/' , wav_name, '_vamp_ua-vamp-plugins_mf0ua_mf0ua.csv']);
    % convert to  onset, offset, hz  
    notes(:,2) = notes(:,1)+notes(:,2);
    notes(:,3) = midi2freq(notes(:,3));

    rmdir(dtemp);
    cd(curdir);
else
    fprintf(['File error: ', input_wav, '\n']);
    h=1;    
end

end