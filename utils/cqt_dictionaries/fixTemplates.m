function [noteTemplates] = fixTemplates(noteTemplates)


% Find existing templates
existingNotes = find(sum(noteTemplates')>0);


% Create auxiliary note template
noteTemplates_aux = noteTemplates + [noteTemplates(:,61:545) zeros(88,60)];


% For each existing note: shift to ideal tuning
for i=1:length(existingNotes)
    
    ind = existingNotes(i);
    
    %if (ind>24)
    if (ind>1)
    
    idealTuning = (ind*5)-4;    
    tuningNeighborhood = [idealTuning-7:idealTuning+9];
    [B,I] = max(noteTemplates_aux(ind,tuningNeighborhood)); % Find peak in tuning neighborhood
    
    % Micro-shift each existing template
    switch I
        case 1
            noteTemplates(ind,:) = [0 0 0 0 0 0 0 noteTemplates(ind,1:end-7)];
        case 2
            noteTemplates(ind,:) = [0 0 0 0 0 0 noteTemplates(ind,1:end-6)];
        case 3
            noteTemplates(ind,:) = [0 0 0 0 0 noteTemplates(ind,1:end-5)];
        case 4
            noteTemplates(ind,:) = [0 0 0 0 noteTemplates(ind,1:end-4)];
        case 5
            noteTemplates(ind,:) = [0 0 0 noteTemplates(ind,1:end-3)];
        case 6
            noteTemplates(ind,:) = [0 0 noteTemplates(ind,1:end-2)];
        case 7
            noteTemplates(ind,:) = [0 noteTemplates(ind,1:end-1)];
        case 8
            noteTemplates(ind,:) = noteTemplates(ind,:);
        case 9
            noteTemplates(ind,:) = [noteTemplates(ind,2:end) 0];
        case 10
            noteTemplates(ind,:) = [noteTemplates(ind,3:end) 0 0];
        case 11
            noteTemplates(ind,:) = [noteTemplates(ind,4:end) 0 0 0];
        case 12
            noteTemplates(ind,:) = [noteTemplates(ind,5:end) 0 0 0 0];
        case 13
            noteTemplates(ind,:) = [noteTemplates(ind,6:end) 0 0 0 0 0];
        case 14
            noteTemplates(ind,:) = [noteTemplates(ind,7:end) 0 0 0 0 0 0];
        case 15
            noteTemplates(ind,:) = [noteTemplates(ind,8:end) 0 0 0 0 0 0 0];
        case 16
            noteTemplates(ind,:) = [noteTemplates(ind,9:end) 0 0 0 0 0 0 0 0];
        case 17
            noteTemplates(ind,:) = [noteTemplates(ind,10:end) 0 0 0 0 0 0 0 0 0];            
    end
    
    end;
    
end;



if 1
    %Apply harmonic comb to remove noise
    harmonicComb = createHarmonicComb();
    activePitches = existingNotes(1):existingNotes(end);
    for i=1:length(activePitches) 
        noteTemplates(activePitches(i),:) = noteTemplates(activePitches(i),:).*harmonicComb(:,activePitches(i))';
    end;
end

% Normalise
for i=1:88 
    noteTemplates(i,:) = noteTemplates(i,:)/(sum(noteTemplates(i,:))+eps);
end


