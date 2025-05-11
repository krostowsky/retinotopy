function [fixedTTLs, adjBehavInd] = correctBehav(indir, subjid, behavFiles, Timestamps, exampleSamples)

collectionTimesStimuli = zeros(1, length(behavFiles));
for q = 1:length(collectionTimesStimuli)
    collectionTimesStimuli(q) = (load([indir '/' subjid '/behav/retinotopy/' behavFiles{q}]).timeframes(end)) * (10^6);
end
collectionTimesBlockOne = load([indir '/' subjid '/behav/retinotopy/' behavFiles{q}]).timeframes;

% this is the length of collection based on the TTL samples
% these should be 1 less than the n
%TimestampsOrig = Timestamps;
%Timestamps = Timestamps(2:end);
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

recordedTime = abs(exampleSamples(1) - Timestamps(1));
estimatedStimTime = (recordedTime - correctionFactor) / 10^6;

stimTimeStart = collectionTimesSamples(1)/10^6 - estimatedStimTime;
adjBehavInd = find((collectionTimesBlockOne >= stimTimeStart) & (collectionTimesBlockOne <= (collectionTimesStimuli(1) / 10^6)));
adjBehavInd = [adjBehavInd(1), adjBehavInd(end)];

end