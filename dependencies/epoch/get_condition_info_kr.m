function timestamps = get_condition_info_kr(subjid)

if strcmpi(subjid, 'UC004') || strcmpi(subjid, 'UC007') 
    timestamps.retinotopy.idx = 1:11;
else
    timestamps.retinotopy.idx = 1:12;
end

end