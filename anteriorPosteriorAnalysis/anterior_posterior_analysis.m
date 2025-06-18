function anterior_posterior_analysis()
outputData = 'statsOut/onOff_LR_significance.mat';
outputData = load(outputData).sig_onOff_LR;

data_dir = '/project/joelvoss/data';

onOff_response = [];
LR_response = [];
LR_response_dir = [];
mni_coords = [];

for j = 1:length(outputData.subjects)
    curr_onOff = outputData.subjects{j}.onOff;
    curr_LR = outputData.subjects{j}.LR;
    curr_LR_dir = outputData.subjects{j}.LRDirection(:,2);

    curr_contacts = outputData.subjects{j}.electrodeLabel;
    curr_MNI = readtable([data_dir '/' outputData.subjID{j} '/anat/' outputData.subjID{j} '_mni.csv' ]);
    [~, matchIdx] = ismember(curr_contacts, curr_MNI.channel);
    curr_MNI = curr_MNI(matchIdx, :);

    sig_on_off = sum(~isnan(curr_onOff), 2) > 0;
    sig_LR = sum(~isnan(curr_LR), 2) > 0;
    sig_LR_dir = curr_LR_dir(sig_LR);

    onOff_response = [onOff_response; sig_on_off];
    LR_response = [LR_response; sig_LR];

    mni_coords = [mni_coords; curr_MNI];
end

[~, ~, stats_onoff] = glmfit([mni_coords.x, mni_coords.y, mni_coords.z], onOff_response, 'binomial','link','logit');

onOff_response_sig = onOff_response(onOff_response == 1);
mni_coords_onOff_sig = mni_coords(onOff_response == 1,:);

LR_sig = LR_response(onOff_response == 1);
[~, ~, stats_lr] = glmfit([mni_coords_onOff_sig.x, mni_coords_onOff_sig.y, mni_coords_onOff_sig.z], LR_sig, 'binomial','link','logit');

end

