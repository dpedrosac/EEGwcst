function idx_bad = remove_badblocks(data_merged, idx_fstrl, code_participant)

%   This function removes blocks in which the participant anticipates the
%   correct answer (code = 40) or is unable to make the right decision (code = 50)

%   Copyright (C) February 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


idx_bad = [];                                           % initialise variable
for m = 1:numel(idx_fstrl) % loop through vector
    if m ~= numel(idx_fstrl)                           % special case of the last condition
        if any(data_merged.trialinfo(idx_fstrl(m):idx_fstrl(m+1))==40) || ...
                any(data_merged.trialinfo(idx_fstrl(m):idx_fstrl(m+1))==140)
            idx_temp = idx_fstrl(m):idx_fstrl(m+1)-1;
            idx_bad = [idx_bad; idx_temp(:)]; clear idx_temp
        end
    else
        if any(data_merged.trialinfo(idx_fstrl(m):...
                numel(data_merged.trialinfo))==40) || ...
                any(data_merged.trialinfo(idx_fstrl(m):...
                numel(data_merged.trialinfo))==140)
            idx_temp = idx_fstrl(m):numel(data_merged.trialinfo);
            idx_bad = [idx_bad; idx_temp(:)]; clear idx_temp
        end
    end
end

idx50 = find(data_merged.trialinfo == 50 | ...
    data_merged.trialinfo == 150);
for k = 1:numel(idx50)
    idx_temp = idx50(k)-2:idx50(k);
    idx_bad = [idx_bad; idx_temp.'];
end
idx_bad = sort(idx_bad, 'ascend');

removed = numel(find(data_merged.trialinfo(idx_bad)==1)) + numel(find(data_merged.trialinfo(idx_bad)==101));
fprintf('\nIn total, %s blocks of data had to be removed for subject %s due to anticipation \nand %s trials because of errors\n', num2str(removed), code_participant, num2str(numel(idx50)))
