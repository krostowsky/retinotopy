function quadrant_analysis_events()

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';
inFilePeriodicData = 'periodicData.mat';
inFilePercents = 'stimulusQuadrantsPercents.mat';
outputData = load('statsOut/onOff_LR_significance.mat');

subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

sig_LR_quadrant = struct('subjects', []);
sig_LR_quadrant.subjects = cell(length(subjects),1);

totSigQuadrant = 0;
totSigQuadrant_LT = 0;
totSigQuadrant_HT = 0;
totSigQuadrant_Gamma = 0;

for j = 1:length(subjects)
    mkdir([figDir '/' subjects{j} '/quadrant_stats']);
    percentData = load([indir '/' subjects{j} '/' inFilePercents], 'imagePercents').imagePercents;
    percentData = cat(1, percentData{:});
    periodicData = load([indir '/' subjects{j} '/' inFilePeriodicData]);

    currSigLT = outputData.sig_onOff_LR.subjects{j}.LR(:,1);
    currSigHT = outputData.sig_onOff_LR.subjects{j}.LR(:,2);
    currSigGamma = outputData.sig_onOff_LR.subjects{j}.LR(:,3);

    currSigLT_dir = outputData.sig_onOff_LR.subjects{j}.LRDirection(:,1);
    currSigHT_dir = outputData.sig_onOff_LR.subjects{j}.LRDirection(:,2);
    currSigGamma_dir = outputData.sig_onOff_LR.subjects{j}.LRDirection(:,3);

    sig_LR_quadrant.subjects{j}.('LR') = nan(length(periodicData.resultLabel), 3);
    sig_LR_quadrant.subjects{j}.('LRDirection') = nan(length(periodicData.resultLabel), 3);
    sig_LR_quadrant.subjects{j}.('quadrant') = cell(length(periodicData.resultLabel), 3);
    sig_LR_quadrant.subjects{j}.('quadrantDirection') = cell(length(periodicData.resultLabel), 3);
    sig_LR_quadrant.subjects{j}.('electrodeLabel') = periodicData.resultLabel;

    for k = 1:length(periodicData.resultLabel)
        if isempty(periodicData.electrodeOnPeriodicData{k,1})
            continue;
        end
        currDataLT_on = [periodicData.electrodeOnPeriodicData{k,1}];
        currDataLT_off = [periodicData.electrodeOffPeriodicData{k,1}];
        currDataLT = [currDataLT_on, currDataLT_off];

        currDataHT_on = [periodicData.electrodeOnPeriodicData{k,2}];
        currDataHT_off = [periodicData.electrodeOffPeriodicData{k,2}];
        currDataHT = [currDataHT_on, currDataHT_off];

        currDataGamma_on = [periodicData.electrodeOnGammaData{k,1}];
        currDataGamma_off = [periodicData.electrodeOffGammaData{k,1}];
        currDataGamma = [currDataGamma_on, currDataGamma_off];

        descriptorData_on = periodicData.modOnStuff{k};
        descriptorData_off = periodicData.modOffStuff{k};
        descriptorData = [descriptorData_on; descriptorData_off];
        [~, sortInd] = sort(descriptorData(:,3));
        descriptorData = descriptorData(sortInd, :);

        currDataLT = currDataLT(sortInd)';
        currDataHT = currDataHT(sortInd)';
        currDataGamma = currDataGamma(sortInd)';

        percentData_on = percentData(descriptorData_on(:,5), 2:end);
        percentData_off = zeros(length(descriptorData_off), size(percentData_on,2));
        percentData_total = [percentData_on; percentData_off];
        percentData_total = percentData_total(sortInd,:);

        % descriptorData(:,5) = percentData_total;

        if currSigLT(k) < 0.05
            sig_LR_quadrant.subjects{j}.('LR')(k,1) = currSigLT(k);
            sig_LR_quadrant.subjects{j}.('LRDirection')(k,1) = currSigLT_dir(k);

            % this fits each quadrant separately
            [~, ~, stats] = arrayfun(@(i) glmfit(percentData_total(:, i), currDataLT), ...
                1:size(percentData_total, 2), ...
                'UniformOutput', false);
            quad_pvals = cellfun(@(s) s.p(2), stats);
            quad_pvals_dir = cellfun(@(s) s.beta(2), stats);

            if any(quad_pvals < 0.05)
                temp_nan = nan(1,4);
                temp_nan(quad_pvals < 0.05) = quad_pvals(quad_pvals < 0.05);
                temp_nan_dir = nan(1,4);
                temp_nan_dir(quad_pvals < 0.05) = quad_pvals_dir(quad_pvals < 0.05);

                sig_LR_quadrant.subjects{j}.('quadrant'){k,1} = temp_nan;
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,1} = temp_nan_dir;
                %checkAssumptions([figDir '/' subjects{j} '/quadrant_stats'], stats.resid, 'quadrant', periodicData.resultLabel{k});
            end
        end

        if currSigHT(k) < 0.05
            sig_LR_quadrant.subjects{j}.('LR')(k,2) = currSigHT(k);
            sig_LR_quadrant.subjects{j}.('LRDirection')(k,2) = currSigHT_dir(k);

            [~, ~, stats] = arrayfun(@(i) glmfit(percentData_total(:, i), currDataHT), ...
                1:size(percentData_total, 2), ...
                'UniformOutput', false);
            quad_pvals = cellfun(@(s) s.p(2), stats);
            quad_pvals_dir = cellfun(@(s) s.beta(2), stats);

            if any(quad_pvals < 0.05)
                temp_nan = nan(1,4);
                temp_nan(quad_pvals < 0.05) = quad_pvals(quad_pvals < 0.05);
                temp_nan_dir = nan(1,4);
                temp_nan_dir(quad_pvals < 0.05) = quad_pvals_dir(quad_pvals < 0.05);

                sig_LR_quadrant.subjects{j}.('quadrant'){k,2} = temp_nan;
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,2} = temp_nan_dir;
                %checkAssumptions([figDir '/' subjects{j} '/quadrant_stats'], stats.resid, 'quadrant', periodicData.resultLabel{k});
            end
        end

        if currSigGamma(k) < 0.05
            sig_LR_quadrant.subjects{j}.('LR')(k,3) = currSigGamma(k);
            sig_LR_quadrant.subjects{j}.('LRDirection')(k,3) = currSigGamma_dir(k);
            [~, ~, stats] = arrayfun(@(i) glmfit(percentData_total(:, i), currDataGamma), ...
                1:size(percentData_total, 2), ...
                'UniformOutput', false);
            quad_pvals = cellfun(@(s) s.p(2), stats);
            quad_pvals_dir = cellfun(@(s) s.beta(2), stats);

            if any(quad_pvals < 0.05)
                temp_nan = nan(1,4);
                temp_nan(quad_pvals < 0.05) = quad_pvals(quad_pvals < 0.05);
                temp_nan_dir = nan(1,4);
                temp_nan_dir(quad_pvals < 0.05) = quad_pvals_dir(quad_pvals < 0.05);

                sig_LR_quadrant.subjects{j}.('quadrant'){k,3} = temp_nan;
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,3} = temp_nan_dir;
                %checkAssumptions([figDir '/' subjects{j} '/quadrant_stats'], stats.resid, 'quadrant', periodicData.resultLabel{k});
            end
        end
        if ~isempty([sig_LR_quadrant.subjects{j}.quadrant{k,:}])
            totSigQuadrant = totSigQuadrant + 1;
        end
    end
end
sig_LR_quadrant.('subjID') = subjects;
sig_LR_quadrant.('totSigOnOff') = totSigOnOff;
sig_LR_quadrant.('totSigOnOff_LT') = totSigOnOff_LT;
sig_LR_quadrant.('totSigOnOff_HT') = totSigOnOff_HT;
sig_LR_quadrant.('totSigOnOff_Gamma') = totSigOnOff_Gamma;
sig_LR_quadrant.('totSigLR') = totSigLR;
sig_LR_quadrant.('totSigLR_LT') = totSigLR_LT;
sig_LR_quadrant.('totSigLR_HT') = totSigLR_HT;
sig_LR_quadrant.('totSigLR_Gamma') = totSigLR_Gamma;
sig_LR_quadrant.('totElect') = totElect;

save(['statsOut/LR__quadrant_signifiance.mat'],'sig_LR_quadrant');

    function checkAssumptions(figDir, residuals, inLabel_freq, inLabel_contact)
        qqplot(residuals); saveas(gcf, [figDir '/qqplot-' inLabel_freq '-' inLabel_contact '.png']); close();
        indFig = figure('visible', 'Off'); scatter(1:length(residuals), residuals); title('residual independence plot'); xlabel('index'); ylabel('residual'); saveas(indFig, [figDir '/independence-' inLabel_freq '-' inLabel_contact '.png']); close();
    end

end