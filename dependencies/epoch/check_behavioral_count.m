function check_behavioral_count(behavFiles, cellCount)
behavDataCheck = behavFiles{cellCount};
[checkRunStart, checkRunEnd] = regexp(behavDataCheck,'run[0-9]{2,}');
behavDataCheck = behavDataCheck(checkRunStart:checkRunEnd);
behavDataCheck = regexp(behavDataCheck,'\d*','Match');
behavDataCheck = str2num(behavDataCheck{1});
if behavDataCheck ~= cellCount
    error('behavioral data (trial) run does not match current count, fix this');
end

end