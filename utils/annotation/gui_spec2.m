function varargout = gui_spec2(varargin)
% GUI_SPEC MATLAB code for gui_spec.fig
%      GUI_SPEC, by itself, creates a new GUI_SPEC or raises the existing
%      singleton*.
%
%      H = GUI_SPEC returns the handle to a new GUI_SPEC or the handle to
%      the existing singleton*.
%
%      GUI_SPEC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPEC.M with the given input arguments.
%
%      GUI_SPEC('Property','Value',...) creates a new GUI_SPEC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_spec_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_spec_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_spec2

% Last Modified by GUIDE v2.5 15-Dec-2017 01:07:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_spec2_OpeningFcn, ...
    'gui_OutputFcn',  @gui_spec2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_spec2 is made visible.
function gui_spec2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_spec2 (see VARARGIN)

% Choose default command line output for gui_spec2
handles.output = hObject;

handles.startNote = [];                
handles.endNote = [];
handles.note_list = {};

axes(handles.log_spec);
hold off; imagesc([]);


if 1
set(hObject, ...
    'WindowButtonDownFcn',   @mouseDownCallback, ...
    'WindowButtonUpFcn',     @mouseUpCallback,   ...
    'WindowButtonMotionFcn', @mouseMotionCallback);

[FileName,PathName,FilterIndex] = uigetfile('*.wav','Load wav file');
wavfile = [PathName,FileName];
handles.wavfile = wavfile;

[y,Fs, mix, x, window, noverlap, nfft] = load_wav(wavfile);
axes(handles.log_spec);
handles.sec = length(y)/Fs;

wav_player = audioplayer(y,Fs);
objdata=[];


minFreq = 130.8128;  % midi = 48
min_midi = freq2midi(minFreq) ;
bins_div = 5;  % bins per semitone

 maxFreq = Fs/2;
 bins = 12*bins_div;  % bins per octave
 sparKernel= sparseKernel(minFreq, maxFreq, bins, Fs);
 intCQT = schramm_cqt(mix,sparKernel);

% [yy,f,t,p] = spectrogram(mix, window, noverlap, nfft, Fs, 'yaxis');
% surf(t,f,(abs(p)),'EdgeColor','none'); 
%  colorbar
%  axis tight;
%  view(0, 90);
%  colormap(jet);
%  set(handles.log_spec,'Yscale','log')
 
intCQT = circshift(intCQT,[-2,0]);

imagesc([0 length(mix)/Fs], [1 size(intCQT,1)/bins_div]+min_midi -0.4,10*log( intCQT )); axis xy;
tstep = (length(mix)/Fs) / size(intCQT,2);
handles.time_step =  tstep;
handles.bins_div = bins_div;
handles.min_midi = min_midi;

colorbar;
colormap('jet');
xlabel('time (s)');
ylabel('midi scale (Hz)');
axis xy;
hold on;

handles.intSpec = intCQT;
handles = load_annotation(handles); 


wav_player.UserData = handles;
handles.wav_player = wav_player;
wav_player.TimerPeriod = .05;
wav_player.TimerFcn = {@fcn_show_time_pointer, wav_player, handles.log_spec, objdata};
wav_player.StopFcn = {@fcn_player_stop, wav_player, handles};
wav_player.StartFcn = {@fcn_player_start, wav_player, handles};

handles.Fs = Fs;
annotation_wave(hObject, eventdata, handles);
handles.wav_player2 = load_wavplayer2(handles);
handles.wav_player3 = load_wavplayer3(handles);



end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_spec2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function wav_player2 = load_wavplayer2(handles)
objdata2=[];

[w,Fs, mix, x, window, noverlap, nfft] = load_wav('annotation.wav');
wav_player2 = audioplayer(w,Fs);
wav_player2.UserData = handles;

wav_player2.TimerPeriod = .05;
wav_player2.TimerFcn = {@fcn_show_time_pointer, wav_player2, handles.log_spec, objdata2};
wav_player2.StopFcn = {@fcn_player_stop, wav_player2, handles};
wav_player2.StartFcn = {@fcn_player_start, wav_player2, handles};

function wav_player3 = load_wavplayer3(handles)
objdata3=[];
[z,Fs, mix, x, window, noverlap, nfft] = load_wav('mix.wav');
wav_player3 = audioplayer(z,Fs);
wav_player3.UserData = handles;


wav_player3.TimerPeriod = .05;
wav_player3.TimerFcn = {@fcn_show_time_pointer, wav_player3, handles.log_spec, objdata3};
wav_player3.StopFcn = {@fcn_player_stop, wav_player3, handles};
wav_player3.StartFcn = {@fcn_player_start, wav_player3, handles};






function [y,Fs, mix, x, window, noverlap, nfft] = load_wav(wfile)
[y,Fs] = audioread(wfile);
[l,c] = size(y);

left = y(:,1);
if(l == 2)
right = y(:,2);
mix = (left + right)/2;
end
mix = left;
x=round(0.10*Fs); %samples taken over 96 msec (20 ms)?
window=hamming(x); %window with a size of 512 points (2205) ?
%noverlap=512; %the noverlaps its the number of points for repeating the window
noverlap=x-round(0.05*Fs);
nfft=4096; %size of the fft
% --- Outputs from this function are returned to the command line.
function varargout = gui_spec2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function play_mix_bt_Callback(hObject, eventdata, handles)

handles.wav_player3 = load_wavplayer3(handles);
guidata(hObject, handles);
wav_player3 = handles.wav_player3;

if wav_player3.isplaying
    stop(wav_player3);
    %handles.play_bt.String = 'Play';
else
    %play(wav_player3);
    startStopPos = xlim(handles.log_spec) * wav_player3.SampleRate;  
    if (startStopPos(1)<1) startStopPos(1)=1; end;
    if (startStopPos(2)>wav_player3.TotalSamples) startStopPos(2)=wav_player3.TotalSamples; end;    
    play(wav_player3, startStopPos);
    %handles.play_bt.String = 'Stop';
end
     

function play_annotation_bt_Callback(hObject, eventdata, handles)

handles.wav_player2 = load_wavplayer2(handles);
guidata(hObject, handles);

wav_player2 = handles.wav_player2;

if wav_player2.isplaying
    stop(wav_player2);
    %handles.play_bt.String = 'Play';
else
    startStopPos = xlim(handles.log_spec) * wav_player2.SampleRate;      
    if (startStopPos(1)<1) startStopPos(1)=1; end;
    if (startStopPos(2)>wav_player2.TotalSamples) startStopPos(2)=wav_player2.TotalSamples; end;    
    play(wav_player2, startStopPos);
    %handles.play_bt.String = 'Stop';
end


% --- Executes on button press in play_bt.
function play_bt_Callback(hObject, eventdata, handles)
% hObject    handle to play_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wav_player = handles.wav_player;

if wav_player.isplaying
    stop(wav_player);
    
else
    startStopPos = xlim(handles.log_spec) * wav_player.SampleRate;  
    if (startStopPos(1)<1) startStopPos(1)=1; end;
    if (startStopPos(2)>wav_player.TotalSamples) startStopPos(2)=wav_player.TotalSamples; end;    
    play(wav_player, startStopPos);
    
end



% --- Executes on key press with focus on play_bt and none of its controls.
function play_bt_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to play_bt (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%% do nothing...

function fcn_player_stop(obj, ...            % refers to the object that called this function (necessary parameter for all callback functions)
    eventdata, ...      % this parameter is not used but is necessary for all callback functions)
    player, ...
    handles)

hMarker = findobj(handles.log_spec, 'Color', 'b');
delete(hMarker);
handles.play_bt.String = 'Stop';


function fcn_player_start(obj, ...            % refers to the object that called this function (necessary parameter for all callback functions)
    eventdata, ...      % this parameter is not used but is necessary for all callback functions)
    player, ...
    handles)
handles.play_bt.String = 'Play';


function fcn_show_time_pointer(...
    obj, ...            % refers to the object that called this function (necessary parameter for all callback functions)
    eventdata, ...      % this parameter is not used but is necessary for all callback functions
    player, ...         % we pass the audioplayer object to the callback function
    figHandle, ...      % pass the figure handle also to the callback function
    objdata)           % finally, we pass the data necessary to draw the new marker

%t = (player.CurrentSample) / player.TotalSamples / 0.1;
%if (t<.3) return; end;

%idx = (player.CurrentSample) / player.TotalSamples/0.1 + .2;
%idx = (player.CurrentSample) / player.TotalSamples;
idx = (player.CurrentSample) / 44100;
%[player.CurrentSample, player.TotalSamples, idx],


if strcmp(player.Running, 'on')
    % get the handle of current marker and delete the marker
    hMarker = findobj(figHandle, 'Color', 'b');
    delete(hMarker);
    stem(figHandle,idx, 100, 'b');
end

function mouseMotionCallback (object, eventdata)
C = get (gca, 'CurrentPoint');
magic_yoffset = 0; %//TODO!
title(gca, ['(X,Y) = (', num2str(C(1,1)), 's, ',num2str(C(1,2)-magic_yoffset), ' midi    Hz = ', num2str(midi2freq(C(1,2)-magic_yoffset)), ')']);


function mouseDownCallback (object, eventdata)
handles = guidata(object);

if ~strcmp(get(object,'SelectionType'),'normal')
    return;
end

C = get (gca, 'CurrentPoint');

%handles.startNote,
if isempty(handles.startNote)
    handles.startNote = [C(1,1), C(1,2)];
    plot(handles.startNote(1), handles.startNote(2), 'g+');    
else
    handles.endNote = [C(1,1), C(1,2)];
   
   size(handles.intSpec),
    %% mapping of input frequency to midi scale (and VQT bins)    
    % gets the middle point (between onset and offset times)
    time_frame = round((handles.startNote(1) + (handles.endNote(1) - handles.startNote(1))/2) /handles.time_step),    
    freq_pos_m = (C(1,2)-handles.min_midi)*handles.bins_div;
    w = ceil(handles.bins_div/2); %duvida nesse valor de w numero de picos por vez  %//TODO           
    
    pitch_pos = [];
    cc=0;
    
    if handles.startNote(1)>handles.endNote(1)
        aux = handles.startNote;
        handles.startNote = handles.endNote;
        handles.endNote = aux;
    end
    pitch_pos(1) = round(freq_pos_m);
    
    %track_peak(handles.intSpec, round(handles.startNote(1)/handles.time_step), freq_pos_m, round(handles.endNote(1)/handles.time_step),freq_pos_m,15);
    
    for i=round(handles.startNote(1)/handles.time_step) : round(handles.endNote(1)/handles.time_step); 
        if i>0 || i<=size(handles.intSpec,2)
            cc=cc+1;
            
            %pitch_pos(cc) = adjust_freq_to_peak(handles.intSpec(:,i), round(freq_pos_m ), w, handles.bins_div);
            pitch_pos(cc) = adjust_freq_to_peak(handles.intSpec(:,i), round(freq_pos_m ), w);
        end
    end
    
    new_freq_pos = median(pitch_pos); %% //TODO +2 is important! why? 5 bins? = 2+1+2 
    %new_freq_pos = adjust_freq_to_peak(handles.intSpec(:,time_frame), round(freq_pos_m ), w),
    pitch_midi = new_freq_pos/handles.bins_div + handles.min_midi;
    %%    
    x  = handles.startNote(1);   
    x2 = handles.endNote(1);    
    
    handles.startNote=[];    
    hMarker = findobj(object, 'Color', 'g');
    delete(hMarker);
    
    %[C(1,2), pitch_midi],
    ht = text(x-0.00,pitch_midi+0.4,sprintf('%s / %2.2f / %2.2f',getNoteChar(round(pitch_midi)), pitch_midi,midi2freq(pitch_midi)),'Color',[0.0 0.0 0.0]);
    hl = line([x,x2],[pitch_midi, pitch_midi],  'LineWidth',2,'Color',[0.0 .0 .0]);    
    %hl2 = line([x,x2],[pitch_midi+24, pitch_midi+24],  'LineWidth',2,'Color',[0.0 .9 .0]);    
    set(hl,'ButtonDownFcn',@fcn_remove_note);  % Now click on a line...
    
    idx = length(handles.note_list);
    handles.note_list{idx+1} = {hl, ht};
end
guidata(object, handles);


function mouseUpCallback (object, eventdata)
%% do nothing...

function fcn_remove_note(object, eventdata)
handles = guidata(object);
if ~strcmp(class(object), 'matlab.graphics.primitive.Line')
    return;
end

if strcmp(get(handles.figure1,'SelectionType'),'alt')
    for i=1:length(handles.note_list)
        obj = handles.note_list{i};
        if isempty(obj)
            continue;
        end
        if (obj{1}.XData(1)==object.XData(1) && obj{1}.XData(2)==object.XData(2))
            delete(obj{1});
            delete(obj{2});
            handles.note_list{i}={};
            guidata(handles.figure1, handles);
            return;
        end
    end
end

function annotation_wave(hObject, eventdata, handles)
% GERAR ONDA SENOIDAL BASEADO NAS ANOTACOES E TOCA-LA

% define common parameters
sec = handles.sec;
fs = handles.Fs;                              % sampling frequency (Hz)
amplitudes = [0.5 0.25 0.25];
phases = [ pi/3 0 pi/7 ];               % pure tone phase (rad/sec)
fade_durations = [ 30 10 ];             % fade-in and fade-out durations (ms)

f_input = [handles.wavfile(1:end-4), '.txt'];

if exist(f_input, 'file')


B = importdata(f_input,',');
A = sortrows(B,1);
num = size(A);
numLines = num(1);  


durationsilence1 = handles.sec*1000;
frequencysilence = 0;
     
frequencies = [frequencysilence frequencysilence*2 frequencysilence*3];
Numcols = floor(durationsilence1*1E-3*fs);            % signal length (samples)
     
            % fade-in and fade-out window function handles
fade_windows = { @(N)(hanning(N).^2) @(N)(chebwin(N,100)) };
Total = tone_generator( fs, durationsilence1, amplitudes, frequencies, phases, fade_durations, fade_windows );


for i=1:numLines
    

        frequencysilence = 0;
        frequencies = [frequencysilence frequencysilence*2 frequencysilence*3];
        durationsilence =  A(i,1)*1000;
        N = floor(durationsilence*1E-3*fs);            % signal length (samples)
        % fade-in and fade-out window function handles
        fade_windows = { @(N)(hanning(N).^2) @(N)(chebwin(N,100)) };
        silence = tone_generator( fs, durationsilence, amplitudes, frequencies, phases, fade_durations, fade_windows );
        
  
        duration = (A(i,2) - A(i,1))*1000;    
        frequency = A(i,3);
        N = floor(duration*1E-3*fs);            % signal length (samples)
        % fade-in and fade-out window function handles
        fade_windows = { @(N)(hanning(N).^2) @(N)(chebwin(N,100)) };

        frequencies = [frequency frequency*2 frequency*3];
        tones = tone_generator( fs, duration, amplitudes, frequencies, phases, fade_durations, fade_windows );
        Signal = horzcat(silence, tones); 
        
        Signal(:,end+1:Numcols+1)=0;   
        
        Total = (Signal + Total);
                
   
   
end

Total = Total/numLines;

audiowrite('annotation.wav', Total, fs);

[y1,Fs, mix1, x, window, noverlap, nfft] = load_wav([handles.wavfile(1:end-4), '.wav']);
[y2,Fs, mix2, x, window, noverlap, nfft] = load_wav('annotation.wav');

%y1=y1*2;
%y2=y2*.5;

mix1(numel(y2)) = 0;
y3 = horzcat(mix1,y2);

audiowrite('mix.wav', y3, Fs);

else
    audiowrite('annotation.wav', zeros(1), fs);
    audiowrite('mix.wav', zeros(1), fs);
end



% --- Executes on button press in save_bt.
function save_bt_Callback(hObject, eventdata, handles)
% hObject    handle to save_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
%handles = sortrows(handles,1)
%[handles.wavfile(1:end-4), '.txt'],

magic_yoffset = 0.0; %//TODO!

fid = fopen([handles.wavfile(1:end-4), '.txt'], 'w');
mirex_notes = [];
for i=1:length(handles.note_list)
    obj = handles.note_list{i};
    if isempty(obj) || isnan(obj{1}.YData(1)-magic_yoffset)
        continue
    end
    
    mirex_notes = [mirex_notes; obj{1}.XData(1),  obj{1}.XData(2),  obj{1}.YData(1)-magic_yoffset];
end

if size(mirex_notes,1)>0
    mirex_notes = sortrows(mirex_notes,1);
end

for i=1:size(mirex_notes,1)
    
    %% save to file
    % onset, offset, pitch
    fprintf(fid,'%2.2f, %2.2f, %f\n', ...
        mirex_notes(i,1),  mirex_notes(i,2),  midi2freq(mirex_notes(i,3))); % saving in Hz
    
    %% display on console
    % onset, offset, pitch
    fprintf('%2.2fs, %2.2fs, %fHz\n', ...
        mirex_notes(i,1),  mirex_notes(i,2),  midi2freq(mirex_notes(i,3))); % saving in Hz
    
    
    
end

fclose(fid);
B = dlmread([handles.wavfile(1:end-4), '.txt']);
%B = load([handles.wavfile(1:end-4), '.txt']);
A = sortrows(B,1);
dlmwrite([handles.wavfile(1:end-4), '.txt'],A);
annotation_wave(hObject, eventdata, handles);

% reload annotations (wav audio) 
[w,Fs, mix, x, window, noverlap, nfft] = load_wav('annotation.wav');
handles.wav_player2 = audioplayer(w,Fs);
[z,Fs, mix, x, window, noverlap, nfft] = load_wav('mix.wav');
handles.wav_player3 = audioplayer(z,Fs);
guidata(hObject, handles);



function handles = load_annotation(handles)
%(hObject, eventdata, handles)
f_input = [handles.wavfile(1:end-4), '.txt'];

magic_yoffset = 0; %//TODO!

if exist(f_input, 'file')
    fid = fopen(f_input,'r');
    A = textscan(fid, '%f, %f, %f\n');
    A = [A{1}, A{2}, A{3}];
    fclose(fid);
    
    for i=1:size(A,1);                        
        hz = A(i,3);        
        pitch_midi = freq2midi(hz) + magic_yoffset; %//TODO!;        
        
        x = A(i,1);
        x2 = A(i,2);        
        ht = text(x-0.00,pitch_midi+0.4,sprintf('%s / %2.2f / %2.2f',getNoteChar(round(pitch_midi)), pitch_midi,midi2freq(pitch_midi)),'Color',[0.0 0.0 0.0]);
        hl = line([x,x2],[pitch_midi,pitch_midi],  'LineWidth',2,'Color',[0.0 0.0 0.0]);
        %hl2 = line([x,x2],[pitch_midi+36, pitch_midi+36],  'LineWidth',2,'Color',[0.0 .9 .0]);    
        set(hl,'ButtonDownFcn',@fcn_remove_note);  % Now click on a line...
        
        idx = length(handles.note_list);
        handles.note_list{idx+1} = {hl, ht};
        
    end
    
end

%% auxiliary functions
function midi = freq2midi(hz)
midi = 12*log2(hz/440)+69;

function hz = midi2freq(midi)
hz = 440*2^((midi-69)/12);
