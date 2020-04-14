function m = build_harmonica_constraints(max_pitch, num_frames, pitch_ref)

if nargin<3
    pitch_ref = 0; %% this means = 60
end

% pitch_ref =  0 % harmonica in C
% pitch_ref = -5 % harmonica in G
% pitch_ref = -2 % harmonica in Bb

ss = [max_pitch, num_frames];
m = zeros(ss);

ii=1;

b=[];
b(ii,:) = [60 64 67 72 0  0  0  0  0  0 ];ii=ii+1;  
b(ii,:) = [0  64 67 72 76 0  0  0  0  0 ];ii=ii+1;  
b(ii,:) = [0  0  67 72 76 79 0  0  0  0 ];ii=ii+1; 
b(ii,:) = [0  0  0  72 76 79 84 0  0  0 ];ii=ii+1; 
b(ii,:) = [0  0  0  0  76 79 84 88 0  0 ];ii=ii+1; 
b(ii,:) = [0  0  0  0  0  79 84 88 91 0 ];ii=ii+1; 
b(ii,:) = [0  0  0  0  0  0  84 88 91 96];ii=ii+1; 

b(ii,:) = [62 67 71 74 0  0  0  0  0  0 ];ii=ii+1; 
b(ii,:) = [0  67 71 74 77 0  0  0  0  0 ];ii=ii+1; 
b(ii,:) = [0  0  71 74 77 81 0  0  0  0 ];ii=ii+1; 
b(ii,:) = [0  0  0  74 77 81 83 0  0  0 ];ii=ii+1; 
b(ii,:) = [0  0  0  0  77 81 83 86 0  0 ];ii=ii+1; 
b(ii,:) = [0  0  0  0  0  81 83 86 89 0 ];ii=ii+1; 
b(ii,:) = [0  0  0  0  0  0  83 86 89 93];ii=ii+1; 

b(ii,:) = 60;ii=ii+1; 
b(ii,:) = 64;ii=ii+1; 
b(ii,:) = 67;ii=ii+1; 
b(ii,:) = 72;ii=ii+1; 
b(ii,:) = 76;ii=ii+1; 
b(ii,:) = 79;ii=ii+1; 
b(ii,:) = 84;ii=ii+1; 
b(ii,:) = 88;ii=ii+1; 
b(ii,:) = 91;ii=ii+1; 

b(ii,:) = [62];ii=ii+1; 
b(ii,:) = [67];ii=ii+1; 
b(ii,:) = [71];ii=ii+1; 
b(ii,:) = [74];ii=ii+1; 
b(ii,:) = [77];ii=ii+1; 
b(ii,:) = [81];ii=ii+1; 
b(ii,:) = [83];ii=ii+1; 
b(ii,:) = [86];ii=ii+1; 
b(ii,:) = [89];ii=ii+1; 


b(ii,:) = [61];ii=ii+1; 
b(ii,:) = [65];ii=ii+1; 
b(ii,:) = [66];ii=ii+1; 
b(ii,:) = [68];ii=ii+1; 
b(ii,:) = [69];ii=ii+1; 
b(ii,:) = [70];ii=ii+1; 
b(ii,:) = [73];ii=ii+1; 
b(ii,:) = 78;ii=ii+1;  %% THIS IS NOT PRESENT IN THE HARMONICA!!! BUT IT WAS IN THE ANNOTATIONS!!!!!
b(ii,:) = [80];ii=ii+1; 

b(ii,:) = [87];ii=ii+1; 
b(ii,:) = [90];ii=ii+1; 
b(ii,:) = [95];ii=ii+1; 


%%



%b = b-5; % G
%b = b-2; % bb
b = b+pitch_ref;
b(b<0)=0;



a = [1:max_pitch];
for j=1:size(b,1)
    bb = b(j,:);  bb(bb==0)=[];
    idx = ismember(a,bb);
    m(idx,:,j) = 1;
    
    if j>15 && j<34
        m(idx,:,j) = 1.5;
    end
end




