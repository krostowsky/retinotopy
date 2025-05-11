function tfAnalysis_jk_kr(indir, subject)
indir2 = indir;
addpath(genpath('/project/joelvoss/tmp-rostowsky/bosc/'));

load([indir2 '/' subject '/hpcElectrodeData.mat'], '-regexp', '^(?!indir$)...');

fs = datCell.sr;
numFreqs = 80;
% lowerFreqBound = 0.2;
lowerFreqBound = 1;
upperFreqBound = 200;

%% 
[res, ~] = tfAnalysis(hippocampalDataReref, fs, lowerFreqBound, upperFreqBound, indir2, numFreqs, resultLabel, figDir, subject);
close();

%%
fprintf([indir2 '/n']);
save([indir2 '/tfAnalysis.mat' ], '-v7.3');
%delete([indir2 '/centerData.mat']);
fprintf('\nfinished\n');
end

%%
function [res, T] = tfAnalysis(eeg, fs, lowerFreqBound, upperFreqBound, indir, numFreqs, resultLabel, figDir, subject)

    % define sampling rate and frequencies to analyze
    freqs = logspace(log10(lowerFreqBound), log10(upperFreqBound), numFreqs);
    
    % for every contact within the hippocampus, run FOOOF analysis
    res.B = [];

    for currentContact = 1:size(eeg,1)
    
        %% time-frequency analysis
        [B, T, F] = BOSC_tf(eeg(currentContact,:), ... % signal
            freqs, ... %freqs
            round(fs), ... %samplerate
            6); % number of wavelet cycles
    
        %%
        % padding for edge effects
        edge=ceil(6*fs/min(F));
    
        exclude_mask = false(size(eeg(currentContact,:)));
        exclude_mask(1:edge) = true;
        exclude_mask(end-edge:end) = true;
    
        %%
        ps = nanmean(B(:, ~exclude_mask),2);
        ps = log10(ps');
        
        diffUpperBound = F - 40;
        upperBoundInd = find(diffUpperBound == min(diffUpperBound(diffUpperBound > 0 )));
        
        ap_guess = [nan, ps(1)];

        [ap_params, ap_ps] = robust_ap_fit(freqs(1:upperBoundInd), ps(1:upperBoundInd), ap_guess);

        fitFig = figure; plot(F(1:length(ap_ps)), ap_ps(:), 'LineStyle', '-'); hold on; plot(F, ps, 'LineStyle', '--'); 
        axes('position', [.6 .6 .3 .3]);
        box on; plot(freqs(1:57), ap_ps(1:57)); hold on; plot(freqs(1:57), ps(1:57)); legend('ap fit', 'power spectrum');
        saveas(fitFig, [figDir '/' subject '/' resultLabel{currentContact} '-plottedFit.png']);
        close();

        %%      
        res.freqs = F;
        res.ap_params(:,currentContact) = ap_params;
        res.ap_ps(:,currentContact) = (ap_ps);
        res.B = cat(3, res.B, B);
        res.T = T;
        res.ps(:,currentContact) = ps;
        legendEntries{currentContact} = sprintf('%s', resultLabel{currentContact});
    end

    %% peak detection on periodic activity
    res.diff = res.ps(1:upperBoundInd,:) - res.ap_ps;
    res.peakFreq = cell(size(res.diff,2),2);
    
    peakFig = figure; hold on; plot(F(1:upperBoundInd), res.diff); 
    for j = 1:size(res.diff,2)
        [peakVal, peakFreq] = findpeaks(res.diff(:,j), freqs(1:upperBoundInd));
        peakFreq = peakFreq(peakVal > 0);
        peakVal = peakVal(peakVal > 0);
        res.peakFreq{j,1} = peakFreq;
        res.peakFreq{j,2} = peakVal;
        plot(peakFreq, peakVal, '*r')
    end
    
    saveas(peakFig, [figDir '/' subject '/periodicActivityPeaks.png']);
    close();
    hold off
    
    % %% plotting
    % fooofFig = figure; hold on; plot(F(1:upperBoundInd), res.ap_ps); plot(F, res.ps,'--'); 
    % xlabel('Frequency (Hz)'); ylabel('Power (log10))');
    % saveas(fooofFig, [figDir '/' subject '/fooofFig.png']);
    % close();
end