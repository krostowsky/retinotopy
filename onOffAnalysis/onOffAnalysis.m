clc; clear;
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy/dependencies/tfAnalysis/bosc/'));

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

figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';

computeAPPerEvent = 1;

% I extend the bounds a bit to include peaks that would normally be
% excluded
lowTheta = [1 5];
highTheta = [5 11];
gamma = [80, 200];
%gammaFreqs = res.freqs - 80;
gammaFreqs = [67:80]; % these are indices from the TFA

fs = 500;

%% 
removeSpikes = 1;
freqThreshold = 11; % anything over 12Hz is outside of interest or just noise

for j = 1:length(subjects)
    mkdir([figDir '/' subjects{j} '/fitEvents']);
    mkdir([figDir '/' subjects{j} '/fitEventsRatio']);

    mkdir([figDir '/' subjects{j} '/fitPeaks']);
    mkdir([figDir '/' subjects{j} '/fitPeaksRatio']);

    % load on and off data
    load([indir '/' subjects{j} '/' inFile], 'res', 'hippocampalDataReref', 'resultLabel');
    load([indir '/' subjects{j} '/events.mat'], 'descriptorMat', 'onEventIndices', 'onEventTimes', 'offEventIndices', 'offEventTimes', 'onEventImages','out');
    onEventIndices = cat(2, onEventIndices{:});
    offEventIndices = cat(2, offEventIndices{:});
    if removeSpikes
        modMat = descriptorMat((descriptorMat(1:size(onEventIndices,2),5) == 0), :);
        onEventIndices = onEventIndices(:,(descriptorMat(1:size(onEventIndices,2),5) == 0));
        offEventIndices = offEventIndices(:,(descriptorMat(1:size(offEventIndices,2),5) == 0));
    end

    electrodePeakData = cell(size(hippocampalDataReref,1),1);
    electrodePeakDataPower = cell(size(hippocampalDataReref,1),1);
    electrodePeakDataPeriodicPower = cell(size(hippocampalDataReref,1),1);
    electrodePeakDataOff = cell(size(hippocampalDataReref,1),1);
    electrodePeakDataPowerOff = cell(size(hippocampalDataReref,1),1);
    electrodePeakDataPeriodicPowerOff = cell(size(hippocampalDataReref,1),1);
    electrodePeakDataFreqs = cell(size(hippocampalDataReref,1),1);
    for k = 1:size(hippocampalDataReref,1)
        currPowerData = res.B(:,:,k);
        currOffAvg = currPowerData(:, offEventIndices(offEventIndices ~= 0)); % no need to split by event here, assuming all off data is the same
        currOffAvg = log10(nanmean(currOffAvg,2));

        % for every electrode, compare on power to average off power for each event
        if computeAPPerEvent
            mkdir([figDir '/' subjects{j} '/fitEvents/electrode' num2str(k) ]);
            mkdir([figDir '/' subjects{j} '/fitEventsRatio/electrode' num2str(k) ]);

            mkdir([figDir '/' subjects{j} '/fitPeaks/electrode' num2str(k) ]);
            mkdir([figDir '/' subjects{j} '/fitPeaksRatio/electrode' num2str(k) ]);

            %currTotPeaks = []; % used to look at all peaks within the current electrode
            %currTotPeaksThresh = []; % used to look at all peaks below threhsold within the current electrode
            currPeaks = res.peakFreq{k,1};
            if isempty(currPeaks) % means no oscillatory activity
                fprintf(['No peaks for contact: ' num2str(k) ' subject: ' subjects{j} '\n']);
                %continue;
            end

            % avg off data, and mean power for each off event at each frequency
            % you cannot mean across events to recover the total average, the
            % events are log10 values, meaning averaging over them will not
            % give you currOffAvg
            ap_guess = [nan, currOffAvg(1)];
            [ap_params, ap_ps] = robust_ap_fit(res.freqs(17:57), currOffAvg(17:57)', ap_guess);
            currOffAvgPeriodic = currOffAvg(17:57) - ap_ps';

            currOffDataEvents = zeros(size(currOffAvgPeriodic,1), size(offEventIndices,2));
            currOffDataEventsPower = zeros(size(currOffAvgPeriodic,1), size(offEventIndices,2));
            for i = 1:size(offEventIndices, 2)
                currOffDataEvents(:, i) = log10(nanmean(currPowerData(17:57, offEventIndices(offEventIndices(:,i) ~= 0, i)), 2)) - ap_ps';
                currOffDataEventsPower(:, i) = log10(nanmean(currPowerData(17:57, offEventIndices(offEventIndices(:,i) ~= 0, i)), 2));
            end
            offEventOffAvgPeriodicRatios = currOffDataEvents ./ currOffAvgPeriodic;

            % for each on event, find peak freqs
            currOnDataEvents = zeros(length(res.freqs), size(onEventIndices,2));
            currContactPeaks = cell(size(onEventIndices,2), 2); %
            currContactPeaksPower = cell(size(onEventIndices,2), 2); %
            currContactPeaksPeriodicPower = cell(size(onEventIndices,2), 2); %
            currContactPeaksOff = cell(size(onEventIndices,2), 2); %
            currContactPeaksPowerOff = cell(size(onEventIndices,2), 2); %
            currContactPeaksFreq = cell(size(onEventIndices,2), 2); %
            for i = 1:size(onEventIndices,2)
                currEventData = currPowerData(:, onEventIndices(onEventIndices(:,i) ~= 0, i)); % this just gets the raw (non log-transformed) power data for the current event
                currOnDataEvents(:,i) = log10(nanmean(currEventData,2)); % this is the average power for each frequency for the event
                currPs = currOnDataEvents(:,i);

                apFreqs = res.freqs(17:57);
                ap_guess = [nan, currPs(1)];
                [ap_params, ap_ps] = robust_ap_fit(res.freqs(17:57), currPs(17:57)', ap_guess); % fit 1/f for each event

                % visualize fit
                % eventFit = figure('visible','off'); plot(res.freqs, currPs); hold on; plot(res.freqs(17:57), ap_ps);
                % axes('position', [.6 .6 .3 .3]);
                % box on; plot(res.freqs(17:46), currPs(17:46)); hold on; plot(res.freqs(17:46), ap_ps(1:30)); saveas(eventFit, [figDir '/' subjects{j} '/fitEvents/electrode' num2str(k) '/event-' num2str(i) '-fit.png']); close();

                % find peaks between 3 - 42Hz
                [currEventValsOrig, currEventPeaks] = findpeaks(currPs(17:57)' - ap_ps, res.freqs(17:57));

                currEventVals = currEventPeaks(currEventValsOrig > 0 & currEventPeaks <= freqThreshold);
                currEventPeaks = currEventPeaks(currEventValsOrig > 0 & currEventPeaks <= freqThreshold); % this is every peak below 10Hz
                if isempty(currEventPeaks)
                    continue;
                end
                currLtPeak = [];
                currHtPeak = [];
                findLt = find(round(currEventPeaks) <= lowTheta(2) & round(currEventPeaks) > lowTheta(1));
                findHt = find(round(currEventPeaks) <= highTheta(2) & round(currEventPeaks) > highTheta(1));
                if length(findLt) > 1
                    [~, maxInd] = max(currEventVals(findLt));
                    findLt = findLt(maxInd);
                    currLtPeak = currEventPeaks(findLt);
                else
                    currLtPeak = currEventPeaks(findLt);
                end
                if length(findHt) > 1
                    [~, maxInd] = max(currEventVals(findHt));
                    findHt = findHt(maxInd);
                    currHtPeak = currEventPeaks(findHt);
                else
                    currHtPeak = currEventPeaks(findHt);
                end
                currContactPeaksFreq{i,1} = currLtPeak;
                currContactPeaksFreq{i,2} = currHtPeak;

                [~, currEventPeaksIndsLt] = ismember(currLtPeak, res.freqs);
                [~, currEventPeaksIndsHt] = ismember(currHtPeak, res.freqs);

                [~, currEventPeaksIndsAPLt] = ismember(currLtPeak, apFreqs);
                [~, currEventPeaksIndsAPHt] = ismember(currHtPeak, apFreqs);

                if ~isempty(findLt)
                    currContactPeaks{i,1} = (currPs(currEventPeaksIndsLt) - ap_ps(currEventPeaksIndsAPLt)) / currOffAvgPeriodic(currEventPeaksIndsAPLt);
                    currContactPeaksPower{i,1} = currPs(currEventPeaksIndsLt);
                    currContactPeaksPeriodicPower{i,1} = (currPs(currEventPeaksIndsLt) - ap_ps(currEventPeaksIndsAPLt));
                end
                if ~isempty(findHt)
                    currContactPeaks{i,2} = (currPs(currEventPeaksIndsHt) - ap_ps(currEventPeaksIndsAPHt)) / currOffAvgPeriodic(currEventPeaksIndsAPHt);
                    currContactPeaksPower{i,2} = currPs(currEventPeaksIndsHt);
                    currContactPeaksPeriodicPower{i,2} = (currPs(currEventPeaksIndsHt) - ap_ps(currEventPeaksIndsAPHt));
                end
            end
            electrodePeakData{k,1} = currContactPeaks;
            electrodePeakDataPower{k,1} = currContactPeaksPower;
            electrodePeakDataPeriodicPower{k,1} = currContactPeaksPeriodicPower;
            electrodePeakDataOff{k,1} = {nanmean(offEventOffAvgPeriodicRatios(1:9, :),1), nanmean(offEventOffAvgPeriodicRatios(10:20, :),1)};
            electrodePeakDataPowerOff{k,1} = {nanmean(currOffDataEventsPower(1:9, :),1), nanmean(currOffDataEventsPower(10:20, :),1)};
            electrodePeakDataPeriodicPowerOff{k,1} = {nanmean(currOffDataEvents(1:9, :),1), nanmean(currOffDataEvents(10:20, :),1)};
            electrodePeakDataFreqs{k,1} = currContactPeaksFreq;

        else % compare on vs off but look at peak within entire timeseries
            offFreqs = res.freqs(1:57);
            currOffAvgPeriodic = currOffAvg(1:57) - res.ap_ps(:,k);
            currOffDataEvents = zeros(size(currOffAvgPeriodic,1), size(offEventIndices,2));
            %currOffDataEventsGamma = zeros(1, size(offEventIndices,2));

            for i = 1:size(offEventIndices, 2)
                currOffDataEvents(:, i) = log10(nanmean(currPowerData(1:57, offEventIndices(offEventIndices(:,i) ~= 0, i)), 2));
                currOffDataEventsPower(:, i) = log10(nanmean(currPowerData(1:57, offEventIndices(offEventIndices(:,i) ~= 0, i)), 2));
                %averageGammaPower = squeeze(nanmean(zscore(log10(currPowerData(gammaFreqs, :)), 0, 2), 1));
            end
            %offEventOffAvgPeriodicRatios = abs(currOffDataEvents ./ currOffAvgPeriodic);

            currPeaks = res.peakFreq(k,:);
            currLowThetaPeak = find(currPeaks{1} < lowTheta(2) & currPeaks{1} >= lowTheta(1));
            currHighThetaPeak = find(currPeaks{1} < highTheta(2) & currPeaks{1} >= highTheta(1));

            % if there are > 1 peaks for a range, take the one with the
            % highest "peak" (periodic power)
            if length(currLowThetaPeak) > 1
                checkPowers = currPeaks{2}(currLowThetaPeak);
                [~, maxInd] = max(checkPowers);
                currLowThetaPeak = currLowThetaPeak(maxInd);
            end
            if length(currHighThetaPeak) > 1
                checkPowers = currPeaks{2}(currHighThetaPeak);
                [~, maxInd] = max(checkPowers);
                currHighThetaPeak = currHighThetaPeak(maxInd);
            end
            currLowThetaPeak = currPeaks{1}(currLowThetaPeak);
            currHighThetaPeak = currPeaks{1}(currHighThetaPeak);

            eventRatios = cell(size(onEventIndices,2), 3);
            for i = 1:size(onEventIndices,2)
                if ~isempty(currLowThetaPeak)
                    %ltPeriodicRatio = abs((log10(nanmean(currPowerData(find(res.freqs == currLowThetaPeak), onEventIndices(onEventIndices(:,i) ~= 0 ,i))))));% - res.ap_ps(find(offFreqs == currLowThetaPeak),k))  / currOffAvgPeriodic(find(offFreqs == currLowThetaPeak)));
                    ltPeriodicRatio = log10(nanmean(currPowerData(find(res.freqs == currLowThetaPeak), onEventIndices(onEventIndices(:,i) ~= 0 ,i))));% - res.ap_ps(find(offFreqs == currLowThetaPeak),k))  / currOffAvgPeriodic(find(offFreqs == currLowThetaPeak)));
                    eventRatios{i,1} = ltPeriodicRatio;
                end
                if ~isempty(currHighThetaPeak)
                    %htPeriodicRatio = abs((log10(nanmean(currPowerData(find(res.freqs == currHighThetaPeak), onEventIndices(onEventIndices(:,i) ~= 0 ,i))))));% - res.ap_ps(find(offFreqs == currHighThetaPeak),k))  / currOffAvgPeriodic(find(offFreqs == currHighThetaPeak)));
                    htPeriodicRatio = log10(nanmean(currPowerData(find(res.freqs == currHighThetaPeak), onEventIndices(onEventIndices(:,i) ~= 0 ,i))));% - res.ap_ps(find(offFreqs == currHighThetaPeak),k))  / currOffAvgPeriodic(find(offFreqs == currHighThetaPeak)));
                    eventRatios{i,2} = htPeriodicRatio;
                end
            end

            setBinWidth = 0.5;
            lowHigh = figure; %('visible','off');
            subplot(2,1,1);
            if ~isempty(currLowThetaPeak)
                %histogram([eventRatios{:,1}], 'BinWidth', setBinWidth); hold on; histogram(offEventOffAvgPeriodicRatios(find(offFreqs == currLowThetaPeak),:), 'BinWidth', setBinWidth); legend('On:Off periodic power ratios', 'Off:Off periodic power ratios'); title(['Ratio of On and Off events to average off event power: low theta - ' num2str(currLowThetaPeak) 'Hz']);
                histogram([eventRatios{:,1}], 'BinWidth', setBinWidth); hold on; histogram(currOffDataEvents(find(offFreqs == currLowThetaPeak),:), 'BinWidth', setBinWidth); legend('on theta power', 'off theta power'); title(['on vs off theta power - ' num2str(currLowThetaPeak) 'Hz']);
            end
            subplot(2,1,2);
            if ~isempty(currHighThetaPeak)
                %histogram([eventRatios{:,2}], 'BinWidth', setBinWidth); hold on; histogram(offEventOffAvgPeriodicRatios(find(offFreqs == currHighThetaPeak),:), 'BinWidth', setBinWidth); legend('On:Off periodic power ratios', 'Off:Off periodic power ratios'); title(['Ratio of On and Off events to average off event power: high theta - ' num2str(currHighThetaPeak) 'Hz']);
                histogram([eventRatios{:,2}], 'BinWidth', setBinWidth); hold on; histogram(currOffDataEvents(find(offFreqs == currHighThetaPeak),:), 'BinWidth', setBinWidth); legend('on theta power', 'off theta power'); title(['on vs off theta power - ' num2str(currHighThetaPeak) 'Hz']);
            end
            saveas(lowHigh, [figDir '/' subjects{j} '/onOffWholeTSPeak/electrode' num2str(k) 'theta.png']); close();
        end
    end
    mkdir([figDir '/' subjects{j} '/outFigs']);
    if computeAPPerEvent
        % electrodePeakData : electrodePeakDataOff
        ratioCompare = figure; countGraph = 0;
        for jj = 1:length(electrodePeakData)
            [~, edges] = histcounts([abs([electrodePeakData{jj}{:,1}]), abs([electrodePeakDataOff{jj}{1}])]);
            countGraph = countGraph + 1;
            subplot(length(electrodePeakData),2,countGraph); histogram(abs([electrodePeakData{jj}{:,1}]), edges, 'Normalization', 'probability'); hold on; histogram(abs([electrodePeakDataOff{jj}{1}]), edges, 'Normalization', 'probability');
            countGraph = countGraph + 1;
            [~, edges] = histcounts([abs([electrodePeakData{jj}{:,2}]), abs([electrodePeakDataOff{jj}{2}])]);
            subplot(length(electrodePeakData),2,countGraph); histogram(abs([electrodePeakData{jj}{:,2}]), edges, 'Normalization', 'probability'); hold on; histogram(abs([electrodePeakDataOff{jj}{2}]), edges, 'Normalization', 'probability');
        end
        sgtitle('ratios of on periodic power vs off periodic power');
        saveas(ratioCompare,[figDir '/' subjects{j} '/outFigs/ratioOnOff.png']); close();

        % electrodePeakDataPeriodicPower : electrodePeakDataPerioidicPowerOff
        periodicPower = figure;  countGraph = 0;
        for jj = 1:length(electrodePeakDataPeriodicPower)
            [~, edges] = histcounts([[electrodePeakDataPeriodicPower{jj}{:,1}], [electrodePeakDataPeriodicPowerOff{jj}{1}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodePeakData),2,countGraph); histogram(([electrodePeakDataPeriodicPower{jj}{:,1}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodePeakDataPeriodicPowerOff{jj}{1}]), edges, 'Normalization', 'probability');
            [~, edges] = histcounts([[electrodePeakDataPeriodicPower{jj}{:,2}], [electrodePeakDataPeriodicPowerOff{jj}{2}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodePeakData),2,countGraph); histogram(([electrodePeakDataPeriodicPower{jj}{:,2}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodePeakDataPeriodicPowerOff{jj}{2}]), edges, 'Normalization', 'probability');
        end
        sgtitle('on periodic power vs off periodic power');
        saveas(periodicPower,[figDir '/' subjects{j} '/outFigs/periodicPowerOnOff.png']); close();

        % electrodePeakDataPower : electrodePeakDataPowerOff
        powerOnOff = figure; title('on power to off power'); countGraph = 0;
        for jj = 1:length(electrodePeakDataPower)
            [~, edges] = histcounts([[electrodePeakDataPower{jj}{:,1}], [electrodePeakDataPowerOff{jj}{1}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodePeakData),2,countGraph); histogram(([electrodePeakDataPower{jj}{:,1}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodePeakDataPowerOff{jj}{1}]), edges, 'Normalization', 'probability');
            [~, edges] = histcounts([[electrodePeakDataPower{jj}{:,2}], [electrodePeakDataPowerOff{jj}{2}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodePeakData),2,countGraph); histogram(([electrodePeakDataPower{jj}{:,2}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodePeakDataPowerOff{jj}{2}]), edges, 'Normalization', 'probability');
        end
        sgtitle('on power vs off power');
        saveas(powerOnOff,[figDir '/' subjects{j} '/outFigs/onpowervsoffpower.png']); close();
    else
        error('stop');
    end
end

%%