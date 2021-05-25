function print_legend()

%   This function prints the code for the different trials aong with
%   an explanation

%   Copyright (C) May 2021
%   D. Pedrosa, Urs Kleinholdermann, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

toi = {% Different trial number available; for further details cf: ./paradigm/WCST-Tremorparadigma.sce
    {[2,3], 'Early trials (wo/ alcohol)'}, ...
    {[6,7], 'Late trials (wo/ alcohol)'}, ...
    {[102,103], 'Early trials (w/ alcohol)'}, ...
    {[106,107], 'Late trials (w/ alcohol)'}, ...
    {21:25, 'Right trials (wo/ alcohol)'}, ...
    {121:125, 'Right trials (w/ alcohol)'}, ...
    {[10, 20], 'Shift trials(wo/ alcohol)'}, ...
    {[110, 120], 'Shift trials(w/ alcohol)'}};

fprintf('Legend for all trials used during WCST\n')
fprintf('======================================\n')
for k = 1:numel(toi)
   fprintf('%g,', toi{k}{1}(1:end-1))  
   fprintf('%g\t->\t%s\n', toi{k}{1}(end), toi{k}{2})  
end
fprintf('\n')

end