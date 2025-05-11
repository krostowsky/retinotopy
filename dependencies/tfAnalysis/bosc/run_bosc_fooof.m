function res = run_bosc_fooof(dat)

for t = 1:length(dat.trial)

    eeg = squeeze(dat.trial{t});

    for e = 1:size(eeg,1) % contacts
        fs = dat.fsample;
        freqs = logspace(log10(2), log10(40), 50); % was 0.5
        [B, T, F] = BOSC_tf(eeg(e,:), ... % signal
            freqs, ... %freqs
            round(fs), ... %samplerate
            5); % number of wavelet cycles

        % padding for edge effects
        edge=ceil(5*fs/min(F));

        exclude_mask = false(size(eeg(1,:)));
        exclude_mask(1:edge) = true;
        exclude_mask(end-edge:end) = true;

        ps = nanmean(B(:, ~exclude_mask),2);
        ps = log10(ps');

        ap_guess = [nan, ps(1), 0];
        [ap_params, ap_ps] = robust_ap_fit(freqs, ps, [nan, ps(1), 0]);

        ap_ps = 10.^ap_ps;

        [powthresh,durthresh] = BOSC_thresholds(fs, ...
            0.95, ... %
            3, ...
            F, ...
            ap_ps);
        is_osc = nan(size(B));
        for fr = 1:size(B,1)
            is_osc(fr,:) = BOSC_detect(B(fr,:), ...
                powthresh(fr), ...
                durthresh(fr), ...
                fs);
        end

        theta_freq = freqs > 4 & freqs < 9;

        is_theta_burst = any(is_osc(theta_freq,:));% & ~any(is_osc(freqs < 6,:));

        if any(is_theta_burst)
        res.m(t,:,e) = nanmean(log10(B(:, is_theta_burst)),2);
        res.e(t,:,e) = nanstd(log10(B(:, is_theta_burst)),[],2);
        else
        res.m(t,:,e) = nan(size(B,1),1);
        res.e(t,:,e) = nan(size(B,1),1);
        end
        res.freqs = F;

        res.ap_params(t,:,e) = ap_params;
        res.ap_ps(t,:,e) = log10(ap_ps);
        res.ps(t,:,e) = ps;
    end
end