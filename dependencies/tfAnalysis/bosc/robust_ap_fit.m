function [ap_params, ap_ps] = robust_ap_fit(freqs, ps, ap_guess)
% ROBUST_AP_FIT robust aperiodic fit to replace BOSC_bgfit
%
% INPUTS:  freqs - double, frequencies contained in the power spectrum
%
%          freqs - double, frequencies contained in the power spectrum
%
% OUTPUTS: none, this function generates a figure for display

% Copyright (C) 2021  James Kragel

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% ps should be in log10 scale, as is the ap_ps output

% some parameters - could be inputs
% ap_amp_thresh = 2.5;
ap_amp_thresh = 0.2;

if ~exist('ap_guess','var')
    ap_guess = [nan, ps(1)];
end

% initial aperiodic fit
init_params = simple_ap_fit(freqs, ps, ap_guess);
init_params = real(init_params);                        
init_ps = gen_ap(init_params, freqs);

% flatten power spectrum from initial
flat_ps = ps - init_ps;
flat_ps(flat_ps<0)=0;

% amplitude threshold
p_thresh = prctile(flat_ps, ap_amp_thresh);
mask = flat_ps <= p_thresh;
freqs_use= freqs(mask); % frequencies in aperiodic range
ps_use = ps(mask);

% find the function to use
if length(init_params) == 2
    ap_func = @exp_nk_func;
elseif length(init_params) == 3
    ap_func = @exp_func;
end

% number of iterations, bounds
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','off');
options.OptimalityTolerance = 1e-10;
options.StepTolerance = 1e-10;
options.FunctionTolerance = 1e-10;

ap_params = lsqcurvefit(ap_func, ...
    init_params, ...
    freqs_use, ...
    ps_use, ...
    [], ...
    [], ...
    options);

ap_ps = gen_ap(ap_params, freqs);


end


function ap_params = simple_ap_fit(freqs, ps, ap_guess)

if isnan(ap_guess(1))
    ap_guess(1) = ps(1);
end

if length(ap_guess) == 2
    ap_func = @exp_nk_func;
elseif length(ap_guess) == 3
    ap_func = @exp_func;
end

% for now, assume no bounds
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','off');
options.OptimalityTolerance = 1e-10;
options.StepTolerance = 1e-10;
options.FunctionTolerance = 1e-10;

[ap_params, resnorm, resid, exflag] = lsqcurvefit(ap_func, ...
                        ap_guess, ...
                        freqs, ...
                        ps, ...
                        [], ... % lower bound
                        [], ... % upper bound
                        options); 

end

function ys = exp_func(params, xs)
% EXP_FUNC exponential function to return ap shape
% xs - frequncies
% params - 3 set of parameters (offset, knee, exp)

ys = zeros(size(xs));
ys = ys + params(1) - log10(params(2) + xs.^(params(3)));

end

function ys = exp_nk_func(params, xs)
% EXP_NK_FUNC exponential function to return ap shape
% xs - frequncies
% params - 2 set of parameters (offset, exp)

ys = zeros(size(xs));
ys = ys + params(1) - log10(xs.^(params(2)));

end

% function pow = gen_ps(freqs, ap_params, gauss_params, nlv)
% 
% ap = gen_ap(ap_params, freqs);
% peaks = gen_peaks(gauss_params, freqs);
% noise = normrnd(0, nlv, size(freqs));
% 
% pow = 10.^(ap+peaks+noise);
% 
% end

function ap = gen_ap(ap_params, freqs)

if length(ap_params)==2 % no knee
    ap = exp_nk_func(ap_params, freqs);
elseif length(ap_params)==3 % knee
    ap = exp_func(ap_params, freqs);
end

end