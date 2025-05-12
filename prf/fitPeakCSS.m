clc; clear;
addpath(genpath('/scratch/midway3/rostowsky/badsCode/bads-master'));
addpath(genpath('/project/joelvoss/tmp-rostowsky/bosc/'));

%%
origImages = load('/project/joelvoss/code/retinotopy/workspace_retinotopyCaltsmash.mat');
origCircle = origImages.basiccircle(1:4:end, 1:4:end);

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);
inFile = 'tfAnalysis.mat';

allPeaks = [];
subjectData = struct();
subjectData.('subject') = cell(length(subjects),1);

makeTS = 1;
fs = 500;

%%
for j = 1:length(subjects)
    load([indir '/' subjects{j} '/' inFile], 'res', 'hippocampalDataReref');
    load([indir '/' subjects{j} '/events.mat'],'onEventIndices', 'onEventTimes', 'offEventIndices', 'offEventTimes', 'onEventImages','out');
    allPeaks = [allPeaks, [res.peakFreq{:,1}]];
    onEventIndices = cat(2, onEventIndices{:});
    offEventIndices = cat(2, offEventIndices{:});

    for k = 1:size(res.peakFreq,1) % for every contact, find the peak frequencies
        if isempty(res.peakFreq{k,1})
            continue;
        end
        
        if makeTS 
            for kk = 1:length(res.peakFreq{k,1}) % for every peak for the current electrode find the event timeseries for on and its ratio to the average off 
                for i = 1:size(onEventOut, 1) % for each event, calculate average power at peak freq
                   
                    %% for each event  
                    currEventData = res.B(:, onEventIndices(onEventIndices(:,i) ~= 0, i), k);
                    ps = nanmean(currEventData,2);
                    ps = log10(ps');

                    ap_guess = [nan, ps(1)];
                    [ap_params, ap_ps] = robust_ap_fit(res.freqs(12:57), ps(12:57), ap_guess);
                    figure; plot(res.freqs, ps); hold on; plot(res.freqs(12:57), ap_ps);

                    ps = log10(nanmean(res.B(:,onEventIndices((onEventIndices(:,i) ~= 0),k)),2));
                    ap_guess = [nan, ps(1), 0];

                    onEventOut(i) = log10(nanmean(peakData(onEventIndices((onEventIndices(:,i) ~= 0),i)))) - ap_ps(currPeakInd);
                    onEventImagesSeries = cat(3, onEventImages{:});

                end    
                
                %onEventOutRatio = onEventOut ./ avgOff;
            end
        end
    end
    subjectData.('subject'){j} = currBandCell;
end

%%
