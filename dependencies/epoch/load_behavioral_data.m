function behavFiles = load_behavioral_data(indir, subject)
% indir: root subject directory
% subject: subject id in indir with "behav" directory

behavFiles = struct2cell(dir([indir '/' subject '/behav/retinotopy/*.mat']));
behavFiles = behavFiles(1,:);
end