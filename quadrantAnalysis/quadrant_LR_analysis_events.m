function quadrant_LR_analysis_events()

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';
inFilePeriodicData = 'periodicData.mat';
inFilePercents = 'stimulusQuadrants_LR_Percents.mat';
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

            [~, ~, stats] = glmfit([percentData_total(:,1), percentData_total(:,2)], currDataLT);
            quad_pvals = stats.p(2:3);
            quad_pvals_dir = stats.beta(2:3);

            if any(quad_pvals < 0.05)
                totSigQuadrant_LT = totSigQuadrant_LT + 1;
                temp_nan = nan(1,2);
                temp_nan(quad_pvals < 0.05) = quad_pvals(quad_pvals < 0.05);
                temp_nan_dir = nan(1,2);
                temp_nan_dir(quad_pvals < 0.05) = quad_pvals_dir(quad_pvals < 0.05);

                sig_LR_quadrant.subjects{j}.('quadrant'){k,1} = temp_nan;
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,1} = temp_nan_dir;
                checkAssumptions([figDir '/' subjects{j} '/quadrant_stats'], stats.resid, 'lt', periodicData.resultLabel{k});
            else
                sig_LR_quadrant.subjects{j}.('quadrant'){k,1} = [];
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,1} = [];
            end
        end

        if currSigHT(k) < 0.05
            sig_LR_quadrant.subjects{j}.('LR')(k,2) = currSigHT(k);
            sig_LR_quadrant.subjects{j}.('LRDirection')(k,2) = currSigHT_dir(k);

            [~, ~, stats] = glmfit([percentData_total(:,1), percentData_total(:,2)], currDataHT);
            quad_pvals = stats.p(2:3);
            quad_pvals_dir = stats.beta(2:3);

            if any(quad_pvals < 0.05)
                totSigQuadrant_HT = totSigQuadrant_HT + 1;
                temp_nan = nan(1,2);
                temp_nan(quad_pvals < 0.05) = quad_pvals(quad_pvals < 0.05);
                temp_nan_dir = nan(1,2);
                temp_nan_dir(quad_pvals < 0.05) = quad_pvals_dir(quad_pvals < 0.05);

                sig_LR_quadrant.subjects{j}.('quadrant'){k,2} = temp_nan;
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,2} = temp_nan_dir;
                checkAssumptions([figDir '/' subjects{j} '/quadrant_stats'], stats.resid, 'HT', periodicData.resultLabel{k});
            else
                sig_LR_quadrant.subjects{j}.('quadrant'){k,2} = [];
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,2} = [];
            end
        end

        if currSigGamma(k) < 0.05
            sig_LR_quadrant.subjects{j}.('LR')(k,3) = currSigGamma(k);
            sig_LR_quadrant.subjects{j}.('LRDirection')(k,3) = currSigGamma_dir(k);

            [~, ~, stats] = glmfit([percentData_total(:,1), percentData_total(:,2)], currDataGamma);
            quad_pvals = stats.p(2:3);
            quad_pvals_dir = stats.beta(2:3);

            if any(quad_pvals < 0.05)
                totSigQuadrant_Gamma = totSigQuadrant_Gamma + 1;
                temp_nan = nan(1,2);
                temp_nan(quad_pvals < 0.05) = quad_pvals(quad_pvals < 0.05);
                temp_nan_dir = nan(1,2);
                temp_nan_dir(quad_pvals < 0.05) = quad_pvals_dir(quad_pvals < 0.05);

                sig_LR_quadrant.subjects{j}.('quadrant'){k,3} = temp_nan;
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,3} = temp_nan_dir;
                checkAssumptions([figDir '/' subjects{j} '/quadrant_stats'], stats.resid, 'gamma', periodicData.resultLabel{k});
            else
                sig_LR_quadrant.subjects{j}.('quadrant'){k,3} = [];
                sig_LR_quadrant.subjects{j}.('quadrantDirection'){k,3} = [];
            end
        end

        if ~isempty([sig_LR_quadrant.subjects{j}.quadrant{k,:}])
            totSigQuadrant = totSigQuadrant + 1;
        end

    end
end

%%
sig_LR_quadrant.('subjID') = subjects;
sig_LR_quadrant.('totSigQuadrant') = totSigQuadrant;
sig_LR_quadrant.('totSigQuadrant_LT') = totSigQuadrant_LT;
sig_LR_quadrant.('totSigQuadrant_HT') = totSigQuadrant_HT;
sig_LR_quadrant.('totSigQuadrant_Gamma') = totSigQuadrant_Gamma;

save(['statsOut/LR_quadrant_signifiance.mat'], 'sig_LR_quadrant');

end

function checkAssumptions(figDir, residuals, inLabel_freq, inLabel_contact)
% make plots to check independence and normality
figure('visible', 'Off');
subplot(2,1,1); qqplot(residuals);
subplot(2,1,2); scatter(1:length(residuals), residuals); title('residual independence plot'); xlabel('index'); ylabel('residual');
saveas(gcf, [figDir '/normality-independence-assumptions-' inLabel_freq '-' inLabel_contact '-quadrant-LR' '.png']); close();
close();
end