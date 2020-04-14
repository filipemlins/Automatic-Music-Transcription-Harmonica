function frames = convert_notes2frame(notes,int_cqt, num_samples, Fs)
%% TEMP
s = size(notes,1);
v = zeros([size(int_cqt,2),s]);
tstep = num_samples/Fs/size(int_cqt,2);

frames = zeros(size(int_cqt));

for i=1:s;
    if round(notes(i,1)/tstep)==0
        notes(i,1) = tstep;
    end
    if round(notes(i,2)/tstep)>=size(int_cqt,2)
        notes(i,2) = size(int_cqt,2)*tstep;
    end
        
    
     v(round(notes(i,1)/tstep):round(notes(i,2)/tstep),i) = notes(i,3);
        
end
    

for i=1:size(int_cqt,2);
    f = v(i,:); 
    f(f==0)=[];
    
    frames(freq2midi(f)*5,i) = f;
    
end




