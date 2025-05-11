function [fixedTimestamps, adjBehavInd] = correctBehav(indir, subjid, behavFiles, Timestamps, exampleSamples)

Timestamps = Timestamps(1:11);

% this is the collection times of samples based on the stimuli, you only
% need the last recorded timeframe because each behavior file starts
% recording time at 0, thus the last timeframe is the length of collection
collectionTimesStimuli = zeros(1, length(behavFiles));
for q = 1:length(collectionTimesStimuli)
    collectionTimesStimuli(q) = (load([indir '/' subjid '/behav/retinotopy/' behavFiles{q}]).timeframes(end)) * (10^6);
end
collectionTimesBlockOne = load([indir '/' subjid '/behav/retinotopy/' behavFiles{1}]).timeframes;

% this is the length of collection based on the TTL samples
% these should be 1 less than the number of behavior files
% this all assumes that the first recorded TTL is an "off" TTL so it is
% skipped here
collectionTimesSamples = zeros(1, length(behavFiles)-1);
for q = 1:length(collectionTimesSamples)
    collectionTimesSamples(q) = abs(Timestamps((q * 2)) - Timestamps((q * 2) + 1));
end

% there is time between the TTL being sent and actually taking
% effect, during this time samples are still recording
% to minimize the data lost (estimated ~1300 samples) I take a
% mean of the time difference between stimulus on periods and
% their calculated collection times
correctionFactor = mean(collectionTimesSamples - collectionTimesStimuli(2:end));
if correctionFactor < 0
    error('Something went terribly wrong with the correction factor')
end

% exampleSamples(1) is the time of the first ever sample, this assumes that
% the stimulus is being shown and a TTL was simply not sent
fixedTimestamps = [exampleSamples(1), Timestamps];

% recordedTime is the length of the remade block
% here I apply the correction factor to correct for the possible time
% between TTL signal and recording
recordedTime = abs(exampleSamples(1) - Timestamps(1));
estimatedStimTime = (recordedTime - correctionFactor) / 10^6;

% the recording is missing stimTimeStart amounts of information, so we only
% want the behavioral data from stimTimeStart forward
stimTimeStart = collectionTimesSamples(1)/10^6 - estimatedStimTime; 
adjBehavInd = find((collectionTimesBlockOne >= stimTimeStart) & (collectionTimesBlockOne <= (collectionTimesStimuli(1) / 10^6)));
adjBehavInd = [adjBehavInd(1), adjBehavInd(end)];

end