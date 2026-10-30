function export_plot_csv(resultsDir)
% EXPORT_PLOT_CSV  Write plot-ready CSV files for every figure of the paper,
% so the results can be plotted in Excel, Origin, Python, etc.
% The sub-bands are taken from the saved results (paper version: 6-12 GHz,
% six 1-GHz sub-bands 6.5 ... 11.5 GHz), so no band settings are needed here.
%
% Usage (after fr3_inf_analysis has run once):
%   export_plot_csv('/Users/faizanrauf/Desktop/FR3 ICC PAPER DATA /VNA_data/results')
% It is also called automatically at the end of fr3_inf_analysis.
%
% Output files (in <resultsDir>/csv/), one row per point / band:
%   csv_pl_vs_distance.csv     measured PL per point and band     -> Fig. path loss scatter
%   csv_pl_model_curves.csv    FSPL, CI fit and 3GPP vs distance   -> lines in the same figure
%   csv_ple_sigma_vs_freq.csv  n, sigma, alpha, beta per band      -> Fig. PLE / sigma vs frequency
%   csv_pl_minus_3gpp.csv      PL - 3GPP InF-LOS per point         -> raw data for the CDF
%   csv_pl_minus_3gpp_cdf.csv  sorted values + CDF probability     -> Fig. CDF (plot x vs y directly)
%   csv_ds_per_point.csv       RMS delay spread per point and band -> raw data for the box plot
%   csv_ds_boxplot_stats.csv   10/25/50/75/90 % per band           -> Fig. delay-spread box plot

R = load(fullfile(resultsDir, 'fr3_inf_results.mat'));
out = fullfile(resultsDir, 'csv');
if ~exist(out, 'dir'), mkdir(out); end

fc  = R.fc(:)';  fG = fc/1e9;  nB = numel(fc);
P   = R.P;  ok = R.ok(:);  PL = R.PL;  DS = R.DS*1e9;  res = R.res;
D0  = R.cfg.D0;
FS1m = 32.4 + 20*log10(fG);
bandName = arrayfun(@(f) sprintf('%.1fGHz', f), fG, 'UniformOutput', false);
fitm = ok & P.dist(:) >= R.cfg.D_MIN_FIT;

% 1) measured path loss per point
hdr = [{'pt','distance_m','los','used_in_fit'}, strcat('PL_dB_', bandName)];
write_csv(fullfile(out,'csv_pl_vs_distance.csv'), hdr, ...
          [P.pt(:) P.dist(:) P.los(:) fitm PL]);

% 2) model curves vs distance (FSPL, CI fit, 3GPP InF-LOS)
d = (1:0.1:20)';
M = d;  hdr = {'distance_m'};
for k = 1:nB
    fs  = FS1m(k) + 20*log10(d/D0);
    ci  = FS1m(k) + 10*res.LOS.n(k)*log10(d/D0);
    g3  = 31.84 + 21.5*log10(d) + 19*log10(fG(k));
    M   = [M fs ci g3];                                         %#ok<AGROW>
    hdr = [hdr, {['FSPL_dB_' bandName{k}], ['CI_fit_dB_' bandName{k}], ['3GPP_InF_LOS_dB_' bandName{k}]}]; %#ok<AGROW>
end
write_csv(fullfile(out,'csv_pl_model_curves.csv'), hdr, M);

% 3) parameters vs frequency
hdr = {'fc_GHz','N','n_CI','sigma_CI_dB','alpha_FI_dB','beta_FI','sigma_FI_dB', ...
       'n_freespace','n_3GPP','sigma_3GPP_dB'};
write_csv(fullfile(out,'csv_ple_sigma_vs_freq.csv'), hdr, ...
    [fG' res.LOS.N(:) res.LOS.n(:) res.LOS.sig(:) res.LOS.alpha(:) res.LOS.beta(:) ...
     res.LOS.sigFI(:) 2*ones(nB,1) 2.15*ones(nB,1) 4.3*ones(nB,1)]);

% 4) PL - 3GPP InF-LOS per point, and its CDF
dPL = NaN(size(PL));
for k = 1:nB
    dPL(:,k) = PL(:,k) - (31.84 + 21.5*log10(P.dist(:)) + 19*log10(fG(k)));
end
dPL(~fitm,:) = NaN;
write_csv(fullfile(out,'csv_pl_minus_3gpp.csv'), ...
          [{'pt','distance_m'}, strcat('dPL_dB_', bandName)], [P.pt(:) P.dist(:) dPL]);

nMax = size(PL,1);  C = NaN(nMax, 2*nB);  hdr = {};
for k = 1:nB
    x = sort(dPL(~isnan(dPL(:,k)),k));  n = numel(x);
    C(1:n, 2*k-1) = x;  C(1:n, 2*k) = (1:n)'/n;
    hdr = [hdr, {['dPL_dB_' bandName{k}], ['CDF_' bandName{k}]}]; %#ok<AGROW>
end
% 3GPP model reference: Gaussian with sigma = 4.3 dB
xg = linspace(-15, 15, nMax)';  cg = 0.5*erfc(-xg/(4.3*sqrt(2)));
C = [C xg cg];  hdr = [hdr, {'dPL_dB_3GPP_model','CDF_3GPP_model'}];
write_csv(fullfile(out,'csv_pl_minus_3gpp_cdf.csv'), hdr, C);

% 5) delay spread per point and box-plot statistics
DSo = DS;  DSo(~fitm,:) = NaN;
write_csv(fullfile(out,'csv_ds_per_point.csv'), ...
          [{'pt','distance_m'}, strcat('DS_ns_', bandName)], [P.pt(:) P.dist(:) DSo]);
S = NaN(nB, 8);
for k = 1:nB
    x = sort(DSo(~isnan(DSo(:,k)),k));
    if isempty(x), continue; end
    S(k,:) = [fG(k) numel(x) pct(x,10) pct(x,25) median(x) pct(x,75) pct(x,90) mean(x)];
end
S = [S repmat(R.ds3gpp.LOS, nB, 1)];
write_csv(fullfile(out,'csv_ds_boxplot_stats.csv'), ...
    {'fc_GHz','N','p10_ns','p25_ns','median_ns','p75_ns','p90_ns','mean_ns','3GPP_median_ns'}, S);

fprintf('Plot-ready CSV files for %d sub-bands (%.1f-%.1f GHz centres) written to "%s".\n', ...
        nB, fG(1), fG(end), out);
end

% ---------------------------------------------------------------------------
function q = pct(x, p)
% percentile with linear interpolation (no toolbox needed); x sorted
n = numel(x);
if n == 1, q = x; return; end
r = 1 + (n-1)*p/100;  lo = floor(r);  hi = min(lo+1, n);
q = x(lo) + (r-lo)*(x(hi)-x(lo));
end

function write_csv(fn, hdr, M)
fid = fopen(fn, 'w');
fprintf(fid, '%s', hdr{1});  fprintf(fid, ',%s', hdr{2:end});  fprintf(fid, '\n');
for i = 1:size(M,1)
    s = sprintf('%.4g,', M(i,:));
    s = strrep(s(1:end-1), 'NaN', '');       % empty cell instead of NaN for Excel
    fprintf(fid, '%s\n', s);
end
fclose(fid);
end
