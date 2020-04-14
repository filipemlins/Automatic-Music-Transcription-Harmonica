
function [p, expression] = rwc_path_filter(path_wav, voice_type, vowel, dynamic, var_type) 
% [p, expression] = rwc_path_filter(path_wav, voice_type, vowel, dynamic, var_type) 
%
%  This  function lists the path of files which are in agreement with the input options (filters)
%  input:: 
%   path_wav: starting point for path search. If an empty matrix is passed, then the function uses the default (hard coded) starting point.
%   voice_type: one of the follow possibilities: {'soprano', 'alto', 'tenor', 'baritone', 'bass'}.
%   vowel: one of the follow possibilities: {'a', 'e', 'i', 'o', 'u'}.
%   dynamic: one of the follow possibilities: {'p', 'f', 'm', '*'}, (stantding for piano, mezzo-forte and forte).
%   var_type: one of the follow possibilities: {'n', 's', 'f', '*'}, (stantding for normal, staccato and vibrato).
%  output:: 
%   p: path list
%   expression: the regular expression create by using the filter options.
%
% Created by Rodrigo Schramm on 27/09/2016.

if nargin < 5
    error('Not enough input arguments.');
end

if isempty(path_wav)
    path_wav  = '/Users/schramm/Desktop/RWA_dataset/';
end
%dynamic = {F,M,P}
switch lower(voice_type)
    case {'soprano', 'alto', 'tenor', 'baritone', 'bass'}
        voice_type = lower(voice_type);
    otherwise
        warning('unrecognized voice_type. Options are: soprano, alto, tenor, baritono, bass.');
end

switch lower(vowel)
    case {'a', 'e', 'i', 'o', 'u'}
        vowel = upper(vowel);
    otherwise
        warning('unrecognized vowel. Options are: A, E, I, O, U.');
end

switch lower(dynamic)
    case {'p', 'f', 'm', '*'}
        dynamic = upper(dynamic);
    otherwise
        warning('unrecognized dynamic. Options are: P, M, F (stantding for piano, mezzo-forte and forte).');
end

switch lower(var_type)
    case {'n', 1}
        var_type = '1';
    case {'s',2}
        var_type = '2';
    case {'v',3}
        var_type = '3';   
    case '*'
        var_type = '*';
    otherwise
        warning('unrecognized style. Options are: n, s, v (stantding for normal, staccato and vibrato).');
end

expression = [path_wav, '/', voice_type, '/*/*', vowel, var_type, dynamic, '*.WAV'];
p = rdir(expression);

