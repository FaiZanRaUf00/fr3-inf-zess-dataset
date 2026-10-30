function fr3_inf_analysis()
% FR3_INF_ANALYSIS  Sub-band-resolved FR3 (measured 6-17 GHz, analysed 6-12 GHz) indoor-factory channel analysis
%
% Pipeline (per measurement point ptN.s2p):
%   1) load S21 (complex)                     5) per 1 GHz sub-band: path loss
%   2) antenna-gain correction per frequency  6) per sub-band: PDP + RMS delay spread
%   3) optional thru normalisation            7) CI / FI path-loss fits per band (LOS/NLOS)
%   4) SNR mask from noise-floor measurement  8) 3GPP TR 38.901 InF references + figures
%
% Inputs (in DATA_DIR):
%   ptN.s2p            measurement files (N = point number)
%   positions.csv      columns: pt, distance_m, los        (los: 1 = LOS, 0 = NLOS)
%                      or:      pt, x_m, y_m, los          (distance computed from TX_XY)
%   gain_tx.csv        columns: freq_GHz, gain_dBi   (BBHA 9120 D calibration sheet, Tx)
%   gain_rx.csv        columns: freq_GHz, gain_dBi   (Rx)
%   noise.s2p          (optional) noise-floor measurement, same settings
%
% NLOS (optional, cfg.NLOS_DIR): ptN.s2p (+ ptNbest.s2p) and positions_nlos.csv in their own folder;
%   NLOS points are added to the fits, subplots and Table II (with their own Tx-Rx distances).
%
% Outputs (in OUT_DIR): per-point CSV, fitted-parameter CSV + LaTeX table, PDF figures.
%
% Works in MATLAB (no toolbox needed) and GNU Octave.

%% ======================= CONFIGURATION ==================================
REPO_DIR       = fullfile(fileparts(mfilename('fullpath')), '..');   % repository root
cfg.DATA_DIR   = fullfile(REPO_DIR, 'data', 'LOS');   % LOS ptN.s2p + positions.csv
cfg.OUT_DIR    = fullfile(REPO_DIR, 'results');
cfg.PREFIX     = 'pt';

cfg.F_START    = 6e9;                 % analysis range [Hz]
cfg.F_STOP     = 12e9;                % analysis up to 12 GHz -> 6 sub-bands (paper version)
                                       % (12-17 GHz failed the 1-m check / large scatter;
                                       %  set to 17e9 once the remeasured data passes check_1m)
cfg.BAND_W     = 1e9;                 % sub-band width [Hz]  -> 6 bands for 6-12 GHz

cfg.GAIN_TX    = fullfile('..', 'antenna', 'gain.csv');          % BBHA 9120 L datasheet, 3 m gain ('' = no correction)
cfg.GAIN_RX    = fullfile('..', 'antenna', 'gain.csv');          % same antenna model for Tx and Rx
cfg.THRU_FILE  = '';                  % e.g. 'thru_day1.s2p' to normalise, '' = off
cfg.NOISE_FILE = '';                  % e.g. 'noise.s2p'; '' = no SNR mask
cfg.PT_SELECT  = [];                  % process ONLY these points, e.g. 46:111 ([] = all points)
                                       % 46-111 = pure LOS positions
cfg.ALL_LOS    = true;                % true = treat ALL selected points as LOS
                                       % los column: 1 = LOS, 0 = NLOS (steel), 2 = OLOS (plastic fence)
                                       % OLOS points are processed but kept out of the LOS/NLOS fits
cfg.SNR_MIN_DB = 10;                  % per-frequency-point SNR threshold [dB]
cfg.MIN_VALID  = 0.8;                 % min. fraction of valid points per band, else NaN

cfg.TX_XY      = [26.737 11.147];      % Tx1 position [m] (origin: bottom-left hall corner)
cfg.D_MIN_FIT  = 1.0;                 % exclude d < D_MIN_FIT from fits (near field) [m]
cfg.D0         = 1.0;                 % CI reference distance [m]

cfg.DS_THR_DB  = 25;                  % PDP dynamic range for RMS delay spread [dB]
cfg.DS_NOISE_DB = 5;                  % ... and at least this far above the noise floor [dB] (as NYU)
cfg.ZP         = 4;                   % zero-padding factor for IFFT (interpolation only)

cfg.HALL_LWH   = [32.73 23.29 10.0];  % hall length, width, height [m] -> 3GPP DS (V/S)
cfg.EXCLUDE_BANDS = [];              % sub-band centres [GHz] to drop ([] = none), e.g. bands where the
                                       % 1-m reference deviates > 3 dB from free space
cfg.REP_BANDS  = [6.5 7.5 8.5 9.5 10.5 11.5]; % band centres [GHz] shown in the PL-vs-d and CDF figures
cfg.D_PLOT_MAX = 0;                   % subplot figure: 0 = only the measured range (no extrapolation)
                                       % e.g. 40 = extend fits dashed to 40 m
cfg.N_REF      = [];                  % reference PLE lines, e.g. [1.5 2 2.5 3]; [] = none

% ---- NLOS campaign (robot blockage, separate Tx position) ----
cfg.NLOS_DIR      = fullfile(REPO_DIR, 'data', 'NLOS');  % NLOS .s2p files + positions_nlos.csv ('' = no NLOS data)
cfg.NLOS_POS      = 'positions_nlos.csv';  % in NLOS_DIR; columns pt,x_m,y_m,distance_m,los,tx
cfg.NLOS_FILE_FMT = 'pt%d.s2p';       % boresight file name of NLOS point N
cfg.NLOS_BEST_FMT = '';               % "NLOS best" file name, e.g. 'pt%dbest.s2p' ('' = not measured)
cfg.NLOS_PT_SELECT = [];              % use only these NLOS points ([] = all)
%% ========================================================================

if isempty(cfg.DATA_DIR), cfg.DATA_DIR = pwd; cfg.OUT_DIR = fullfile(pwd, 'results'); end
fprintf('Data folder:   %s\nResults folder: %s\n', cfg.DATA_DIR, cfg.OUT_DIR);
if ~exist(fullfile(cfg.DATA_DIR, 'positions.csv'), 'file')
    error(['positions.csv not found in "%s".\nPut fr3_inf_analysis.m in the same folder as ' ...
           'the ptN.s2p files, or set cfg.DATA_DIR to that folder.'], cfg.DATA_DIR);
end
if ~exist(cfg.OUT_DIR, 'dir')
    [okm, msg] = mkdir(cfg.OUT_DIR);
    if ~okm
        error('Cannot create "%s" (%s). Set cfg.OUT_DIR to a writable folder, e.g. your Desktop.', ...
              cfg.OUT_DIR, msg);
    end
end

% ---------- positions ----------
P = read_positions(fullfile(cfg.DATA_DIR, 'positions.csv'), cfg);
if ~isempty(cfg.PT_SELECT)                    % keep only the selected points
    keep = ismember(P.pt, cfg.PT_SELECT);
    for fn = fieldnames(P)'
        P.(fn{1}) = P.(fn{1})(keep);
    end
    fprintf('PT_SELECT: using %d of %d points (pt %d-%d).\n', sum(keep), numel(keep), ...
            min(P.pt), max(P.pt));
end
if cfg.ALL_LOS
    P.los(:) = 1;
    fprintf('ALL_LOS = true: all points treated as LOS.\n');
end
fprintf('Positions: %d points (%d LOS, %d NLOS, %d OLOS)\n', numel(P.pt), ...
        sum(P.los==1), sum(P.los==0), sum(P.los==2));

% ---------- band edges ----------
edges = cfg.F_START:cfg.BAND_W:cfg.F_STOP;
nB    = numel(edges) - 1;
fc    = (edges(1:end-1) + edges(2:end)) / 2;            % band centres [Hz]
fprintf('Analysis range %.0f-%.0f GHz -> %d sub-bands, centres: %s GHz\n', ...
        cfg.F_START/1e9, cfg.F_STOP/1e9, nB, mat2str(fc/1e9));

% ---------- reference data (gain, thru, noise) -- interpolated per file ----------
gTx = []; gRx = [];
if ~isempty(cfg.GAIN_TX), gTx = read_gain(fullfile(cfg.DATA_DIR, cfg.GAIN_TX)); end
if ~isempty(cfg.GAIN_RX), gRx = read_gain(fullfile(cfg.DATA_DIR, cfg.GAIN_RX)); end
if isempty(gTx) || isempty(gRx)
    warning('No antenna gain correction applied -> results include antenna gains!');
end
thru = []; noise = [];
if ~isempty(cfg.THRU_FILE),  [thru.f,  thru.H ] = read_s21(fullfile(cfg.DATA_DIR, cfg.THRU_FILE));  end
if ~isempty(cfg.NOISE_FILE) && exist(fullfile(cfg.DATA_DIR, cfg.NOISE_FILE), 'file')
    [noise.f, noise.H] = read_s21(fullfile(cfg.DATA_DIR, cfg.NOISE_FILE));
else
    if ~isempty(cfg.NOISE_FILE), warning('Noise file not found -> no SNR mask.'); end
end

% ---------- per-point processing (LOS campaign) ----------
P.set = ones(numel(P.pt), 1);            % 1 = LOS campaign, 2 = NLOS campaign
[PL, DS, DSfb, VAL, ok, KF] = process_points(cfg, cfg.DATA_DIR, ...
    @(n) sprintf('%s%d.s2p', cfg.PREFIX, n), P.pt, edges, gTx, gRx, thru, noise);
fprintf('LOS campaign: processed %d / %d points.\n', sum(ok), numel(P.pt));

% ---------- NLOS campaign (own folder, own positions, all los = 0) ----------
nL   = numel(P.pt);
PLb  = nan(nL, nB);  DSb = nan(nL, nB);  DSfbb = nan(nL, 1);  okb = false(nL, 1);
if ~isempty(cfg.NLOS_DIR)
    Pn = read_positions(fullfile(cfg.NLOS_DIR, cfg.NLOS_POS), cfg);
    if ~isempty(cfg.NLOS_PT_SELECT)
        keep = ismember(Pn.pt, cfg.NLOS_PT_SELECT);
        for fnm = fieldnames(Pn)', Pn.(fnm{1}) = Pn.(fnm{1})(keep); end
    end
    Pn.los(:) = 0;
    [PLn, DSn, DSfbn, VALn, okn, KFn] = process_points(cfg, cfg.NLOS_DIR, ...
        @(n) sprintf(cfg.NLOS_FILE_FMT, n), Pn.pt, edges, gTx, gRx, thru, noise);
    fprintf('NLOS campaign (boresight): processed %d / %d points.\n', sum(okn), numel(Pn.pt));
    nN = numel(Pn.pt);
    PLnb = nan(nN, nB);  DSnb = nan(nN, nB);  DSfbnb = nan(nN, 1);  oknb = false(nN, 1);
    if ~isempty(cfg.NLOS_BEST_FMT)
        [PLnb, DSnb, DSfbnb, ~, oknb] = process_points(cfg, cfg.NLOS_DIR, ...
            @(n) sprintf(cfg.NLOS_BEST_FMT, n), Pn.pt, edges, gTx, gRx, thru, noise);
        fprintf('NLOS campaign (best direction): processed %d / %d points.\n', sum(oknb), nN);
    end
    % append the NLOS points to the LOS points (x/y only kept for the LOS floor plan)
    P.pt   = [P.pt;   Pn.pt];    P.dist = [P.dist; Pn.dist];
    P.los  = [P.los;  Pn.los];   P.set  = [P.set;  2*ones(nN, 1)];
    if isfield(P, 'x'), P.x = [P.x; nan(nN,1)]; P.y = [P.y; nan(nN,1)]; end
    PL  = [PL;  PLn];  DS  = [DS;  DSn];  DSfb  = [DSfb;  DSfbn];  VAL = [VAL; VALn];  ok = [ok; okn];
    KF  = [KF;  KFn];
    PLb = [PLb; PLnb]; DSb = [DSb; DSnb]; DSfbb = [DSfbb; DSfbnb]; okb = [okb; oknb];
end
nP = numel(P.pt);

% ---------- sanity check: nearest point vs free space ----------
[dmin, i0] = min(P.dist + 1e9*~(ok & P.set == 1));   % LOS campaign only
FSd = 32.4 + 20*log10(fc/1e9) + 20*log10(dmin);
fprintf('\nSanity check, nearest point pt%d (d = %.2f m): measured PL - free-space PL [dB]\n', P.pt(i0), dmin);
fprintf('  %5.1f GHz: %+5.1f\n', [fc/1e9; PL(i0,:) - FSd]);
fprintf('  -> within about +-3 dB is fine; larger offsets point to gain, cable or calibration issues.\n\n');

% ---------- exclude sub-bands (e.g. failed 1-m reference check) ----------
if ~isempty(cfg.EXCLUDE_BANDS)
    keepB = ~ismember(round(fc/1e8), round(cfg.EXCLUDE_BANDS*10));
    fprintf('Excluding sub-bands centred at %s GHz -> %d sub-bands remain.\n\n', ...
            mat2str(fc(~keepB)/1e9), sum(keepB));
    fc = fc(keepB);  PL = PL(:,keepB);  DS = DS(:,keepB);  VAL = VAL(:,keepB);
    PLb = PLb(:,keepB);  DSb = DSb(:,keepB);  KF = KF(:,keepB);
    nB = numel(fc);
end

% ---------- per-band fits ----------
d = P.dist;
fitm = d >= cfg.D_MIN_FIT & ok;
FS1m = 32.4 + 20*log10(fc/1e9);           % FSPL at 1 m [dB]
res  = struct();
res.LOS      = fit_class(PL,  DS,  DSfb,  d, fitm & P.los == 1,  FS1m, cfg);
res.NLOS     = fit_class(PL,  DS,  DSfb,  d, fitm & P.los == 0,  FS1m, cfg);
res.NLOSbest = fit_class(PLb, DSb, DSfbb, d, d >= cfg.D_MIN_FIT & okb & P.los == 0, FS1m, cfg);

% excess loss of the NLOS points over the LOS CI fit at the same distance (= blockage loss)
XL = nan(nP, nB);  XLb = nan(nP, nB);
for k = 1:nB
    plLOS    = FS1m(k) + 10*res.LOS.n(k)*log10(d/cfg.D0);
    XL(:,k)  = PL(:,k)  - plLOS;
    XLb(:,k) = PLb(:,k) - plLOS;
end
mN = ok & P.los == 0;  mNb = okb & P.los == 0;
for k = 1:nB
    res.NLOS.xl_mean(k)     = mean_or_nan(XL(mN & ~isnan(XL(:,k)), k));
    res.NLOS.xl_std(k)      = std_or_nan(XL(mN & ~isnan(XL(:,k)), k));
    res.NLOSbest.xl_mean(k) = mean_or_nan(XLb(mNb & ~isnan(XLb(:,k)), k));
    res.NLOSbest.xl_std(k)  = std_or_nan(XLb(mNb & ~isnan(XLb(:,k)), k));
end
% Rician K-factor statistics (moment method, Greenstein et al. 1999), in dB
for k = 1:nB
    for c = {'LOS', 'NLOS'}
        if strcmp(c{1}, 'LOS'), mk = ok & P.los == 1; else, mk = ok & P.los == 0; end
        x = KF(mk, k);  x = x(~isnan(x));
        res.(c{1}).kf_mean(k) = mean_or_nan(x);
        res.(c{1}).kf_std(k)  = std_or_nan(x);
        res.(c{1}).kf_med(k)  = median_or_nan(x);
    end
end
fprintf('\nRician K-factor (moment method) [dB]:  LOS mean +- std (median) | NLOS mean +- std (median)\n');
for k = 1:nB
    fprintf('  %5.1f GHz:   %5.1f +- %4.1f (%5.1f)   |   %5.1f +- %4.1f (%5.1f)\n', fc(k)/1e9, ...
        res.LOS.kf_mean(k), res.LOS.kf_std(k), res.LOS.kf_med(k), ...
        res.NLOS.kf_mean(k), res.NLOS.kf_std(k), res.NLOS.kf_med(k));
end
if any(mN)
    fprintf('\nNLOS results per sub-band:\n');
    fprintf('  band [GHz]   N   n_NLOS  sig   | n_best  sig_best | excess loss bore / best [dB] | DS mean / 90%% [ns]\n');
    for k = 1:nB
        fprintf('  %5.1f      %3d   %5.2f  %4.1f  |  %5.2f   %4.1f   |   %5.1f +- %4.1f / %5.1f     |  %5.1f / %5.1f\n', ...
            fc(k)/1e9, res.NLOS.N(k), res.NLOS.n(k), res.NLOS.sig(k), res.NLOSbest.n(k), res.NLOSbest.sig(k), ...
            res.NLOS.xl_mean(k), res.NLOS.xl_std(k), res.NLOSbest.xl_mean(k), res.NLOS.ds_mean(k), res.NLOS.ds_p90(k));
    end
end

fprintf('\nRMS delay spread, LOS [ns]:      mean / median / 90%% / 10^E[log10]\n');
for k = 1:nB
    fprintf('  %5.1f GHz (res. %.1f ns): %6.1f / %5.1f / %5.1f / %5.1f\n', fc(k)/1e9, 1e9/cfg.BAND_W, ...
            res.LOS.ds_mean(k), res.LOS.ds_med(k), res.LOS.ds_p90(k), res.LOS.ds_lgmean(k));
end
fprintf('  full band %g-%g GHz (res. %.2f ns): %6.1f / %5.1f / %5.1f / %5.1f\n\n', ...
        cfg.F_START/1e9, cfg.F_STOP/1e9, 1e9/(cfg.F_STOP-cfg.F_START), ...
        res.LOS.dsfb_mean, res.LOS.dsfb_med, res.LOS.dsfb_p90, res.LOS.dsfb_lgmean);

% ---------- 3GPP TR 38.901 InF references ----------
L = cfg.HALL_LWH;  V = prod(L);  S = 2*(L(1)*L(2) + L(1)*L(3) + L(2)*L(3));
VS = V / S;
ds3gpp.LOS  = 10^(log10(26*VS + 14) - 9.35) * 1e9;      % median DS [ns]
ds3gpp.NLOS = 10^(log10(30*VS + 32) - 9.44) * 1e9;
fprintf('3GPP InF median DS (V/S = %.2f m): LOS %.1f ns, NLOS %.1f ns\n', VS, ds3gpp.LOS, ds3gpp.NLOS);

% ---------- outputs ----------
save_point_table(cfg, P, ok, fc, PL, DS, DSfb, VAL, KF);
save_param_table(cfg, fc, res, ds3gpp, P, ok);
save_nlos_table(cfg, fc, P, ok, okb, PL, PLb, XL, XLb, DS, DSb);
make_figures(cfg, P, ok, fc, PL, DS, res, ds3gpp, FS1m, PLb, okb);
save(fullfile(cfg.OUT_DIR, 'fr3_inf_results.mat'), 'cfg', 'P', 'ok', 'fc', 'edges', ...
     'PL', 'DS', 'DSfb', 'VAL', 'res', 'ds3gpp', 'PLb', 'DSb', 'okb', 'XL', 'XLb', 'KF');
export_plot_csv(cfg.OUT_DIR);             % plot-ready CSVs in results/csv/
fprintf('Done. Results in "%s".\n', cfg.OUT_DIR);
end

%% ========================================================================
function [PL, DS, DSfb, VAL, ok, KF] = process_points(cfg, dataDir, nameFcn, pts, edges, gTx, gRx, thru, noise)
% Path loss and RMS delay spread per sub-band for a list of points
nP = numel(pts);  nB = numel(edges) - 1;
PL = nan(nP, nB);  DS = nan(nP, nB);  DSfb = nan(nP, 1);  VAL = nan(nP, nB);  ok = false(nP, 1);
KF = nan(nP, nB);                                % Rician K-factor per band [dB]
for i = 1:nP
    fn = fullfile(dataDir, nameFcn(pts(i)));
    if ~exist(fn, 'file'), warning('Missing %s', fn); continue; end
    [f, H] = read_s21(fn);
    sel = f >= cfg.F_START & f <= cfg.F_STOP;
    f = f(sel); H = H(sel);
    valid = true(size(H));                       % SNR mask (raw S21, before gain correction)
    if ~isempty(noise)
        Pn = interp1(noise.f, abs(noise.H).^2, f, 'linear', 'extrap');
        valid = abs(H).^2 > Pn * 10^(cfg.SNR_MIN_DB/10);
    end
    if ~isempty(thru)                            % thru normalisation (complex)
        H = H ./ interp1(thru.f, thru.H, f, 'linear', 'extrap');
    end
    if ~isempty(gTx) && ~isempty(gRx)            % antenna gain correction: H / sqrt(Gt*Gr)
        G_dB = interp1(gTx.f, gTx.g, f, 'linear', 'extrap') + ...
               interp1(gRx.f, gRx.g, f, 'linear', 'extrap');
        H = H ./ 10.^(G_dB/20);
    end
    for k = 1:nB
        if k < nB, idx = f >= edges(k) & f <  edges(k+1);
        else,      idx = f >= edges(k) & f <= edges(k+1); end
        VAL(i,k) = mean(valid(idx));
        if VAL(i,k) < cfg.MIN_VALID, continue; end
        Hk = H(idx);
        PL(i,k) = -10*log10(mean(abs(Hk(valid(idx))).^2));
        DS(i,k) = rms_ds(Hk, f(idx), cfg);
        KF(i,k) = kfactor_mom(Hk(valid(idx)));
    end
    DSfb(i) = rms_ds(H, f, cfg);
    ok(i) = true;
end
end

function r = fit_class(PL, DS, DSfb, d, m, FS1m, cfg)
% CI and FI fits plus delay-spread statistics for the points selected by mask m
nB = size(PL, 2);
for k = 1:nB
    y = PL(m, k); x = d(m);
    g = ~isnan(y); y = y(g); x = x(g);
    r.N(k) = numel(y);
    if numel(y) < 3
        [r.n(k), r.sig(k), r.alpha(k), r.beta(k), r.sigFI(k)] = deal(NaN);
    else
        D = 10*log10(x/cfg.D0);
        A = y - FS1m(k);
        n = sum(D.*A) / sum(D.^2);               % CI model
        r.n(k)   = n;
        r.sig(k) = std(A - n*D, 1);
        p = [ones(size(D)) D] \ y;               % FI (alpha-beta) model
        r.alpha(k) = p(1);
        r.beta(k)  = p(2);                       % D is already 10log10(d)
        r.sigFI(k) = std(y - [ones(size(D)) D]*p, 1);
    end
    ds = DS(m, k) * 1e9;  ds = ds(~isnan(ds));
    r.ds_med(k)    = median_or_nan(ds);
    r.ds_p10(k)    = prctile_or_nan(ds, 10);
    r.ds_p90(k)    = prctile_or_nan(ds, 90);
    r.ds_mean(k)   = mean_or_nan(ds);            % mean = statistic used by NYU
    r.ds_lgmean(k) = 10^mean_or_nan(log10(ds));  % 10^E[log10 DS] = 3GPP lgDS convention
end
dsf = DSfb(m) * 1e9;  dsf = dsf(~isnan(dsf));   % full-band delay spread
r.dsfb_mean   = mean_or_nan(dsf);
r.dsfb_lgmean = 10^mean_or_nan(log10(dsf));
r.dsfb_med    = median_or_nan(dsf);
r.dsfb_p90    = prctile_or_nan(dsf, 90);
end

function save_nlos_table(cfg, fc, P, ok, okb, PL, PLb, XL, XLb, DS, DSb)
% per-point NLOS results (boresight and best direction) incl. excess loss over the LOS fit
m = P.los == 0 & (ok | okb);
if ~any(m), return; end
fid = fopen(fullfile(cfg.OUT_DIR, 'per_point_nlos.csv'), 'w');
fprintf(fid, 'pt,distance_m');
for k = 1:numel(fc), fprintf(fid, ',PL_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',PLbest_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',XL_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',XLbest_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',DS_ns_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',DSbest_ns_%.1fGHz', fc(k)/1e9); end
fprintf(fid, '\n');
for i = find(m)'
    fprintf(fid, '%d,%.3f', P.pt(i), P.dist(i));
    fprintf(fid, ',%.2f', PL(i,:), PLb(i,:), XL(i,:), XLb(i,:), DS(i,:)*1e9, DSb(i,:)*1e9);
    fprintf(fid, '\n');
end
fclose(fid);
end

function y = std_or_nan(x), if numel(x) < 2, y = NaN; else, y = std(x); end, end

function KdB = kfactor_mom(H)
% Rician K-factor by the moment method (Greenstein, Michelson, Erceg, IEEE Comm. Lett. 1999),
% using the frequency samples of one sub-band (as in Runtti et al., EuCAP 2025):
%   Ga = E[|H|^2],  Gv = std(|H|^2),  K = sqrt(Ga^2-Gv^2) / (Ga - sqrt(Ga^2-Gv^2))
p  = abs(H(:)).^2;
if numel(p) < 8, KdB = NaN; return; end
Ga = mean(p);  Gv = std(p);
if Gv >= Ga, KdB = -20; return; end               % Rayleigh-like (K -> 0): floor at -20 dB
r  = sqrt(Ga^2 - Gv^2);
KdB = max(10*log10(r / (Ga - r)), -20);            % same lower floor as Runtti et al.
end

function tau_rms = rms_ds(Hk, fk, cfg)
% RMS delay spread from a band-limited frequency response (Hann window, IFFT)
N   = numel(Hk);
if N < 8, tau_rms = NaN; return; end
w   = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));               % Hann (no toolbox)
L   = cfg.ZP * N;
h   = ifft(Hk(:) .* w, L);
df  = mean(diff(fk));
tau = (0:L-1)' / (L*df);
pdp = abs(h).^2;
pdp = pdp(1:floor(L/2));  tau = tau(1:floor(L/2));      % drop wrap-around half
% threshold as in NYU (Ying et al., ICC 2025): 25 dB below the peak OR 5 dB above
% the noise floor, whichever is higher. Noise floor estimated from the last 25 %
% of the delay window (> ~375 ns excess delay, where no paths arrive in the hall).
nTail = pdp(round(0.75*numel(pdp)):end);
noise = median(nTail);
thr   = max(max(pdp) * 10^(-cfg.DS_THR_DB/10), noise * 10^(cfg.DS_NOISE_DB/10));
m   = pdp >= thr;
p   = pdp(m); t = tau(m);
t0  = sum(t.*p) / sum(p);
tau_rms = sqrt(sum((t - t0).^2 .* p) / sum(p));
end

%% ------------------------------------------------------------------------
function [f, S21] = read_s21(fn)
% Minimal Touchstone .s2p reader (RI / MA / DB, any frequency unit)
fid = fopen(fn, 'r');  if fid < 0, error('Cannot open %s', fn); end
unit = 1e9; fmt = 'MA'; data = [];
while true
    ln = fgetl(fid);  if ~ischar(ln), break; end
    ln = strtrim(ln);
    if isempty(ln) || ln(1) == '!', continue; end
    if ln(1) == '#'
        tok = upper(strsplit(strtrim(ln(2:end))));
        for t = tok
            switch t{1}
                case 'HZ',  unit = 1;   case 'KHZ', unit = 1e3;
                case 'MHZ', unit = 1e6; case 'GHZ', unit = 1e9;
                case {'RI','MA','DB'}, fmt = t{1};
            end
        end
        continue
    end
    c = strfind(ln, '!'); if ~isempty(c), ln = ln(1:c(1)-1); end
    v = sscanf(ln, '%f')';
    data = [data; v(1:min(9,end)), nan(1, 9-numel(v))]; %#ok<AGROW>
end
fclose(fid);
f = data(:,1) * unit;
a = data(:,4); b = data(:,5);              % S21 = columns 4,5 (S11 S21 S12 S22 order)
switch fmt
    case 'RI', S21 = a + 1j*b;
    case 'MA', S21 = a .* exp(1j*b*pi/180);
    case 'DB', S21 = 10.^(a/20) .* exp(1j*b*pi/180);
end
end

%% ------------------------------------------------------------------------
function G = read_gain(fn)
M = read_numeric_csv(fn);
G.f = M(:,1) * 1e9;   % freq_GHz -> Hz
G.g = M(:,2);         % gain_dBi
end

function P = read_positions(fn, cfg)
fid = fopen(fn, 'r'); if fid < 0, error('Cannot open %s', fn); end
hdr = lower(strtrim(strsplit(fgetl(fid), ',', 'CollapseDelimiters', false)));
fclose(fid);
M = read_numeric_csv(fn);
col = @(name) find(strcmp(hdr, name), 1);
P.pt  = M(:, col('pt'));
P.los = M(:, col('los'));
if ~isempty(col('distance_m'))
    P.dist = M(:, col('distance_m'));
else
    if any(isnan(cfg.TX_XY)), error('positions.csv has no distance_m -> set cfg.TX_XY'); end
    P.x = M(:, col('x_m')); P.y = M(:, col('y_m'));
    P.dist = hypot(P.x - cfg.TX_XY(1), P.y - cfg.TX_XY(2));
end
if ~isempty(col('x_m')), P.x = M(:, col('x_m')); P.y = M(:, col('y_m')); end
end

function M = read_numeric_csv(fn)
fid = fopen(fn, 'r'); fgetl(fid);            % skip header
M = [];
while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    if isempty(strtrim(ln)), continue; end
    v = str2double(strsplit(ln, ',', 'CollapseDelimiters', false));
    M = [M; v]; %#ok<AGROW>
end
fclose(fid);
end

function y = median_or_nan(x), if isempty(x), y = NaN; else, y = median(x); end, end
function y = mean_or_nan(x),   if isempty(x), y = NaN; else, y = mean(x);   end, end
function y = prctile_or_nan(x, p)
if isempty(x), y = NaN; return; end
x = sort(x(:)); n = numel(x);
r = 1 + (n-1)*p/100; lo = floor(r); hi = ceil(r);
y = x(lo) + (r-lo)*(x(hi)-x(lo));
end

%% ------------------------------------------------------------------------
function save_point_table(cfg, P, ok, fc, PL, DS, DSfb, VAL, KF)
fid = fopen(fullfile(cfg.OUT_DIR, 'per_point_results.csv'), 'w');
fprintf(fid, 'pt,distance_m,los');
for k = 1:numel(fc), fprintf(fid, ',PL_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',DS_ns_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',valid_%.1fGHz', fc(k)/1e9); end
for k = 1:numel(fc), fprintf(fid, ',K_dB_%.1fGHz', fc(k)/1e9); end
fprintf(fid, ',DS_ns_fullband\n');
for i = 1:numel(P.pt)
    if ~ok(i), continue; end
    fprintf(fid, '%d,%.3f,%d', P.pt(i), P.dist(i), P.los(i));
    fprintf(fid, ',%.2f', PL(i,:));
    fprintf(fid, ',%.2f', DS(i,:)*1e9);
    fprintf(fid, ',%.2f', VAL(i,:));
    fprintf(fid, ',%.2f', KF(i,:));
    fprintf(fid, ',%.2f\n', DSfb(i)*1e9);
end
fclose(fid);
end

function save_param_table(cfg, fc, res, ds3gpp, P, ok)
% ---- fitted_parameters.csv (all classes) ----
fid = fopen(fullfile(cfg.OUT_DIR, 'fitted_parameters.csv'), 'w');
cls = {'LOS', 'NLOS', 'NLOSbest'};
fprintf(fid, 'band_GHz,fc_GHz');
for c = cls
    fprintf(fid, '%s', strrep(',N_X,n_X,sigma_X_dB,alpha_X,beta_X,DSmed_X_ns,DSmean_X_ns,DSp90_X_ns', 'X', c{1}));
end
fprintf(fid, ',XLmean_NLOS_dB,XLstd_NLOS_dB,XLmean_NLOSbest_dB\n');
for k = 1:numel(fc)
    fprintf(fid, '%g-%g,%.1f', (fc(k)-cfg.BAND_W/2)/1e9, (fc(k)+cfg.BAND_W/2)/1e9, fc(k)/1e9);
    for c = cls
        r = res.(c{1});
        fprintf(fid, ',%d,%.2f,%.2f,%.2f,%.2f,%.1f,%.1f,%.1f', r.N(k), r.n(k), r.sig(k), ...
                r.alpha(k), r.beta(k), r.ds_med(k), r.ds_mean(k), r.ds_p90(k));
    end
    fprintf(fid, ',%.1f,%.1f,%.1f\n', getk(res.NLOS, 'xl_mean', k), getk(res.NLOS, 'xl_std', k), ...
            getk(res.NLOSbest, 'xl_mean', k));
end
fclose(fid);

% ---- Table II (table*): this work + literature + 3GPP, layout of the paper ----
% needs \usepackage{multirow}
nB   = numel(fc);
nLOS = sum(ok & P.los == 1);  nNL = sum(ok & P.los == 0);
dmin = min(P.dist(ok));  dmax = max(P.dist(ok));
num  = @(x, fmt) ternstr(isnan(x), '--', sprintf(fmt, x));
fid = fopen(fullfile(cfg.OUT_DIR, 'table_params.tex'), 'w');
fprintf(fid, '%% Auto-generated by fr3_inf_analysis.m -- needs \\usepackage{multirow}\n');
fprintf(fid, '\\begin{table*}[t]\n\\centering\n');
fprintf(fid, ['\\caption{Directional path-loss and delay-spread parameters per sub-band ($d_0=%g$~m), ' ...
              'compared with other indoor-factory measurements and 3GPP}\n'], cfg.D0);
fprintf(fid, '\\label{tab:params}\n\\footnotesize\\setlength{\\tabcolsep}{3pt}\n');
fprintf(fid, '\\begin{tabular}{|l|c|c|c|c|ccc|ccc|cc|cc|}\n\\hline\n');
fprintf(fid, ['\\multirow{3}{*}{Campaign} & Band & Distance & Links & HPBW\n' ...
              ' & \\multicolumn{3}{c|}{LOS} & \\multicolumn{3}{c|}{NLOS} ' ...
              '& \\multicolumn{2}{c|}{LOS DS [ns]} & \\multicolumn{2}{c|}{NLOS DS [ns]} \\\\\n\\cline{6-15}\n']);
fprintf(fid, ' & & & & & \\multicolumn{2}{c}{CI} & FI & \\multicolumn{2}{c}{CI} & FI & & & & \\\\\n');
fprintf(fid, [' & [GHz] & [m] & (LOS/NLOS) & (Tx/Rx) & $n$ & $\\sigma$ [dB] & $\\beta$ & $n$ & $\\sigma$ [dB] & $\\beta$ ' ...
              '& mean & 90\\,\\%% & mean & 90\\,\\%% \\\\\n\\hline\n']);
if nNL > 0, links = sprintf('%d/%d', nLOS, nNL); else, links = sprintf('%d/--', nLOS); end
for k = 1:nB
    if k == 1
        c1 = sprintf('\\multirow{%d}{*}{\\textbf{This work}}', nB);
        c3 = sprintf('\\multirow{%d}{*}{%s--%s}', nB, trimnum(dmin), trimnum(dmax));
        c4 = sprintf('\\multirow{%d}{*}{%s}', nB, links);
        c5 = sprintf('\\multirow{%d}{*}{(34$^\\circ$/34$^\\circ$)$^\\ast$}', nB);
    else
        c1 = ''; c3 = ''; c4 = ''; c5 = '';
    end
    fprintf(fid, ' %s & %g--%g & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
        c1, (fc(k)-cfg.BAND_W/2)/1e9, (fc(k)+cfg.BAND_W/2)/1e9, c3, c4, c5, ...
        num(res.LOS.n(k), '%.2f'),  num(res.LOS.sig(k), '%.2f'),  num(res.LOS.beta(k), '%.2f'), ...
        num(res.NLOS.n(k), '%.2f'), num(res.NLOS.sig(k), '%.2f'), num(res.NLOS.beta(k), '%.2f'), ...
        num(res.LOS.ds_mean(k), '%.1f'),  num(res.LOS.ds_p90(k), '%.1f'), ...
        num(res.NLOS.ds_mean(k), '%.1f'), num(res.NLOS.ds_p90(k), '%.1f'));
end
fprintf(fid, '\\hline\n');
fprintf(fid, ['\\multirow{2}{*}{\\cite{ying2025inf}}\n' ...
  ' & 6.75  & \\multirow{2}{*}{9--38} & \\multirow{2}{*}{5/7} & (30$^\\circ$/30$^\\circ$) & 1.64 & 1.44 & -- & 2.92 & 7.69  & -- & 11.2 & -- & 11.3 & -- \\\\\n' ...
  ' & 16.95 & & & (15$^\\circ$/15$^\\circ$) & 2.09 & 6.20 & -- & 3.81 & 10.42 & -- & 13.7 & -- & 9.2 & -- \\\\\n\\hline\n' ...
  '\\cite{ju2023factory} & 142 & 5--87 & 30/48 & (8$^\\circ$/8$^\\circ$) & 2.10 & 1.80 & -- & 4.92 & 11.30 & -- & 1.8 & -- & 3.2 & -- \\\\\n\\hline\n']);
fprintf(fid, ['3GPP InF~\\cite{3gpp38901}$^\\dagger$ & 0.5--100 & 1--600 & -- & omni & -- & 4.30 & 2.15 & -- & 5.7--7.2 & 2.55--3.57 ' ...
              '& %.1f$^\\ddagger$ & -- & %.1f$^\\ddagger$ & -- \\\\\n\\hline\n'], ds3gpp.LOS, ds3gpp.NLOS);
fprintf(fid, '\\end{tabular}\n\\\\[3pt]\n\\begin{minipage}{0.95\\textwidth}\\footnotesize\n');
fprintf(fid, '$^\\ast$Typical 3-dB beamwidth $32^\\circ$ (E-plane) / $36^\\circ$ (H-plane) at both ends. ');
fprintf(fid, ['$^\\dagger$Omnidirectional; floating-intercept models (slope given as $\\beta$); ' ...
              'NLOS for the low-BS sub-scenarios InF-SL / InF-DL. ']);
fprintf(fid, '$^\\ddagger$Median for the measured hall (depends on the volume-to-surface ratio).\n');
fprintf(fid, '\\end{minipage}\n\\end{table*}\n');
fclose(fid);
end

function v = getk(r, f, k), if isfield(r, f), v = r.(f)(k); else, v = NaN; end, end
function s = ternstr(c, a, b), if c, s = a; else, s = b; end, end
function s = trimnum(x)
% 1 -> '1', 2.5 -> '2.5', 17.26 -> '17' (distance range for the table)
if abs(x - round(x)) < 0.05, s = sprintf('%d', round(x)); else, s = sprintf('%.1f', x); end
end

%% ------------------------------------------------------------------------
function make_figures(cfg, P, ok, fc, PL, DS, res, ds3gpp, FS1m, PLb, okb)
W = 3.5; H = 2.6;                        % IEEE single column [in]
fGHz = fc/1e9;
cL = [0 0.447 0.741]; cN = [0.850 0.325 0.098];

% ---- Fig: PL vs distance for representative bands ----
fig = newfig(W, H*1.25);
kk = arrayfun(@(b) find(abs(fGHz - b) == min(abs(fGHz - b)), 1), cfg.REP_BANDS);
mk = {'o','s','^','d','v','p','h','>'};   % one marker per representative band
dd = logspace(log10(max(cfg.D_MIN_FIT, 0.8)), log10(max(P.dist)*1.05), 50);
hold on; leg = {}; hh = [];
for j = 1:numel(kk)
    k = kk(j); sh = (j-1)/(numel(kk)-1+eps);
    col = (1-sh)*[0.2 0.4 0.9] + sh*[0.9 0.2 0.1];
    mL = ok & P.los == 1; mN = ok & P.los == 0;
    h1 = semilogx(P.dist(mL), PL(mL,k), mk{j}, 'Color', col, 'MarkerSize', 3);
    semilogx(P.dist(mN), PL(mN,k), mk{j}, 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 3);
    if ~isnan(res.LOS.n(k))
        semilogx(dd, FS1m(k) + 10*res.LOS.n(k)*log10(dd), '-', 'Color', col, 'LineWidth', 1);
    end
    if ~isnan(res.NLOS.n(k))
        semilogx(dd, FS1m(k) + 10*res.NLOS.n(k)*log10(dd), '--', 'Color', col, 'LineWidth', 1);
    end
    % 3GPP TR 38.901 InF-LOS (Table 7.4.1-1)
    h3 = semilogx(dd, pl_3gpp_inf_los(dd, fGHz(k)), ':', 'Color', col, 'LineWidth', 1.2);
    hh(end+1) = h1; leg{end+1} = sprintf('%.1f GHz', fGHz(k)); %#ok<AGROW>
end
hh(end+1) = h3; leg{end+1} = '3GPP InF-LOS';
set(gca, 'XScale', 'log'); xlim([min(dd)*0.95, max(dd)]);
set(gca, 'XTick', [1 2 5 10 20 30]);
xlabel('Tx-Rx distance d [m]'); ylabel('Path loss [dB]');
legend(hh, leg, 'Location', 'northwest', 'FontSize', 7);
title('markers: measured;  solid: CI fit;  dotted: 3GPP InF-LOS', 'FontSize', 7, 'FontWeight', 'normal');
grid on; box on;
savefig_pdf(fig, cfg, 'fig_pl_vs_distance');

% ---- Fig: PL vs distance, one subplot per sub-band (full-width figure*) ----
nB   = numel(fGHz);
if nB > 6, nCol = 4; else, nCol = min(3, nB); end   % 6 bands -> 2 x 3 grid, 7-8 bands -> 2 x 4
nRow = ceil(nB / nCol);
Wf   = 7.16;                                  % IEEE double-column width [in]
fig  = newfig(Wf, 1.9*nRow + 0.5);             % extra top margin for the legend row
axAll = cell(1, nB);
cO   = [0.47 0.67 0.19];                      % OLOS colour
mL = ok & P.los == 1;  mN = ok & P.los == 0;  mO = ok & P.los == 2;  mNb = okb & P.los == 0;
cNb = [0.494 0.184 0.556];                    % NLOS-best colour
dMeas = max(P.dist(ok));                           % largest measured distance
dIn   = logspace(log10(max(cfg.D_MIN_FIT, 0.9)), log10(dMeas), 60);   % measured range
dOut  = logspace(log10(dMeas), log10(max(cfg.D_PLOT_MAX, dMeas)), 40); % extrapolation
dAll  = [dIn dOut(2:end)];
yAll = [PL(ok, :); PLb(okb, :)];  yAll = yAll(:);  yAll = yAll(~isnan(yAll));
yTop = max(FS1m) + 10*max([cfg.N_REF 2.2])*log10(max(dAll));          % top of the reference lines
if isempty(yAll), yl = [40 100]; else, yl = [floor(min(yAll)/5)*5-2, ceil(max([yAll; yTop])/5)*5+2]; end
cGr = [0.6 0.6 0.6];
for k = 1:nB
    subplot(nRow, nCol, k); hold on;
    hp = []; lp = {};
    % reference PLE lines (light grey, labelled at the right edge)
    for nr = cfg.N_REF
        semilogx(dAll, FS1m(k) + 10*nr*log10(dAll), ':', 'Color', cGr, 'LineWidth', 0.6);
        text(max(dAll)*0.97, FS1m(k) + 10*nr*log10(max(dAll)*0.97), sprintf('n=%.1f', nr), ...
             'FontSize', 5.5, 'Color', [0.45 0.45 0.45], 'HorizontalAlignment', 'right', ...
             'VerticalAlignment', 'bottom');
    end
    if any(mL)
        hp(end+1) = semilogx(P.dist(mL), PL(mL,k), 'o', 'Color', cL, 'MarkerSize', 2.5); lp{end+1} = 'LOS';
    end
    if any(mO)
        hp(end+1) = semilogx(P.dist(mO), PL(mO,k), '^', 'Color', cO, 'MarkerSize', 2.5); lp{end+1} = 'OLOS';
    end
    if any(mN)
        hp(end+1) = semilogx(P.dist(mN), PL(mN,k), 's', 'Color', cN, 'MarkerSize', 2.5); lp{end+1} = 'NLOS';
    end
    if any(mNb)
        hp(end+1) = semilogx(P.dist(mNb), PLb(mNb,k), 'd', 'Color', cNb, 'MarkerSize', 2.5); lp{end+1} = 'NLOS best';
    end
    % CI fits: solid in the measured range, dashed beyond (extrapolation)
    if ~isnan(res.LOS.n(k))
        hp(end+1) = semilogx(dIn,  FS1m(k) + 10*res.LOS.n(k)*log10(dIn),  '-',  'Color', cL, 'LineWidth', 1.1);
        lp{end+1} = 'CI LOS';
        if dOut(end) > dMeas
            semilogx(dOut, FS1m(k) + 10*res.LOS.n(k)*log10(dOut), '--', 'Color', cL, 'LineWidth', 1.1);
        end
    end
    if ~isnan(res.NLOS.n(k))
        hp(end+1) = semilogx(dIn,  FS1m(k) + 10*res.NLOS.n(k)*log10(dIn),  '-',  'Color', cN, 'LineWidth', 1.1);
        lp{end+1} = 'CI NLOS';
        if dOut(end) > dMeas
            semilogx(dOut, FS1m(k) + 10*res.NLOS.n(k)*log10(dOut), '--', 'Color', cN, 'LineWidth', 1.1);
        end
    end
    if ~isnan(res.NLOSbest.n(k))
        hp(end+1) = semilogx(dIn,  FS1m(k) + 10*res.NLOSbest.n(k)*log10(dIn), '--', 'Color', cNb, 'LineWidth', 1.1);
        lp{end+1} = 'CI NLOS best';
    end
    hp(end+1) = semilogx(dAll, pl_3gpp_inf_los(dAll, fGHz(k)), 'k:', 'LineWidth', 1.1); lp{end+1} = '3GPP InF-LOS';
    if any(mN)
        plSL = max(33 + 25.5*log10(dAll) + 20*log10(fGHz(k)), pl_3gpp_inf_los(dAll, fGHz(k)));
        hp(end+1) = semilogx(dAll, plSL, '-.', 'Color', [0.45 0.45 0.45], 'LineWidth', 1); lp{end+1} = '3GPP InF-SL';
    end
    % mark the end of the measured range
    if dOut(end) > dMeas
        plot([dMeas dMeas], yl, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    end
    set(gca, 'XScale', 'log', 'XTick', [1 2 5 10 20 40], 'FontSize', 7);
    xlim([min(dAll)*0.95, max(dAll)]); ylim(yl); grid on; box on;
    title(sprintf('%g-%g GHz', (fc(k)-cfg.BAND_W/2)/1e9, (fc(k)+cfg.BAND_W/2)/1e9), ...
          'FontSize', 8, 'FontWeight', 'normal');
    % parameters in the lower-right corner
    s = sprintf('LOS: n = %.2f, \\sigma = %.1f dB', res.LOS.n(k), res.LOS.sig(k));
    if ~isnan(res.NLOS.n(k))
        s = sprintf('%s\nNLOS: n = %.2f, \\sigma = %.1f dB', s, res.NLOS.n(k), res.NLOS.sig(k));
    end
    if ~isnan(res.NLOSbest.n(k))
        s = sprintf('%s\nNLOS best: n = %.2f, \\sigma = %.1f dB', s, res.NLOSbest.n(k), res.NLOSbest.sig(k));
    end
    text(0.97, 0.05, s, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', 'FontSize', 6, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);
    if mod(k-1, nCol) == 0, ylabel('Path loss [dB]', 'FontSize', 8); end
    if k > nB - nCol,       xlabel('Distance d [m]', 'FontSize', 8); end
    axAll{k} = gca;
    if k == 1, hpLeg = hp; lpLeg = lp; end
end
% one horizontal legend above all panels (created last, so MATLAB does not move it):
% shrink every panel slightly towards the bottom to free a strip at the top
for k = 1:nB
    pa = get(axAll{k}, 'Position');
    set(axAll{k}, 'Position', [pa(1), pa(2)*0.90, pa(3), pa(4)*0.90]);
end
% legend in its own invisible axes in the top strip -> never inside a panel
legAx = axes('Position', [0.05 0.935 0.90 0.06], 'Visible', 'off');  hold(legAx, 'on');
hd = zeros(1, numel(hpLeg));
for j = 1:numel(hpLeg)                           % dummy copies of the plotted styles
    hd(j) = plot(legAx, NaN, NaN, 'LineStyle', get(hpLeg(j), 'LineStyle'), ...
        'Marker', get(hpLeg(j), 'Marker'), 'Color', get(hpLeg(j), 'Color'), ...
        'LineWidth', get(hpLeg(j), 'LineWidth'), 'MarkerSize', get(hpLeg(j), 'MarkerSize'));
end
lgd = legend(legAx, hd, lpLeg, 'Orientation', 'horizontal', 'Location', 'north', 'FontSize', 7);
set(lgd, 'Box', 'off');
savefig_pdf(fig, cfg, 'fig_pl_vs_distance_subplots');

% ---- Fig: PLE and sigma vs frequency (key novelty figure) ----
% one shared legend in a strip above both panels; 3GPP references dashed in both panels
fig = newfig(W, H*1.3 + 0.35);
cNb = [0.494 0.184 0.556];
hasN  = any(~isnan(res.NLOS.n));               % only show NLOS if there is NLOS data
hasNb = any(~isnan(res.NLOSbest.n));
fx  = [fGHz(1) fGHz(end)];
ax1 = subplot(2,1,1); hold on;
plot(gx(fGHz), gy(fGHz, res.LOS.n), '-o', 'Color', cL, 'MarkerFaceColor', cL, 'MarkerSize', 3);
if hasN,  plot(gx(fGHz), gy(fGHz, res.NLOS.n), '-s', 'Color', cN, 'MarkerFaceColor', cN, 'MarkerSize', 3); end
if hasNb, plot(gx(fGHz), gy(fGHz, res.NLOSbest.n), '-d', 'Color', cNb, 'MarkerFaceColor', cNb, 'MarkerSize', 3); end
plot(fx, [2 2], 'k:');
plot(fx, [2.15 2.15], '--', 'Color', cL);
if hasN, plot(fx, [2.55 2.55], '--', 'Color', cN); end
ylabel('PLE n'); grid on; box on; xlim([fGHz(1)-0.5 fGHz(end)+0.5]);
ax2 = subplot(2,1,2); hold on;
plot(gx(fGHz), gy(fGHz, res.LOS.sig), '-o', 'Color', cL, 'MarkerFaceColor', cL, 'MarkerSize', 3);
if hasN,  plot(gx(fGHz), gy(fGHz, res.NLOS.sig), '-s', 'Color', cN, 'MarkerFaceColor', cN, 'MarkerSize', 3); end
if hasNb, plot(gx(fGHz), gy(fGHz, res.NLOSbest.sig), '-d', 'Color', cNb, 'MarkerFaceColor', cNb, 'MarkerSize', 3); end
plot(fx, [4.3 4.3], '--', 'Color', cL);
if hasN, plot(fx, [5.7 5.7], '--', 'Color', cN); end
xlabel('Sub-band centre frequency [GHz]'); ylabel('\sigma_{SF} [dB]');
grid on; box on; xlim([fGHz(1)-0.5 fGHz(end)+0.5]);
% make room at the top and put the shared legend there (own invisible axes)
for axc = {ax1, ax2}
    pa = get(axc{1}, 'Position');
    set(axc{1}, 'Position', [pa(1), pa(2)*0.92, pa(3), pa(4)*0.92]);
end
legAx = axes('Position', [0.05 0.90 0.90 0.09], 'Visible', 'off');  hold(legAx, 'on');
hd = []; ld = {};
hd(end+1) = plot(legAx, NaN, NaN, '-o', 'Color', cL, 'MarkerFaceColor', cL, 'MarkerSize', 3); ld{end+1} = 'LOS';
if hasN,  hd(end+1) = plot(legAx, NaN, NaN, '-s', 'Color', cN, 'MarkerFaceColor', cN, 'MarkerSize', 3); ld{end+1} = 'NLOS'; end
if hasNb, hd(end+1) = plot(legAx, NaN, NaN, '-d', 'Color', cNb, 'MarkerFaceColor', cNb, 'MarkerSize', 3); ld{end+1} = 'NLOS best'; end
hd(end+1) = plot(legAx, NaN, NaN, 'k:');             ld{end+1} = 'free space';
hd(end+1) = plot(legAx, NaN, NaN, '--', 'Color', cL); ld{end+1} = '3GPP InF-LOS';
if hasN, hd(end+1) = plot(legAx, NaN, NaN, '--', 'Color', cN); ld{end+1} = '3GPP InF-SL'; end
lgd = legend(legAx, hd, ld, 'Orientation', 'horizontal', 'Location', 'north', 'FontSize', 6);
set(lgd, 'Box', 'off');
try, set(lgd, 'NumColumns', 3); catch, end            % 2 rows x 3 entries (MATLAB >= R2018a)
savefig_pdf(fig, cfg, 'fig_ple_sigma_vs_freq');

% ---- Fig: RMS delay spread vs frequency ----
fig = newfig(W, H); hold on;
errbar(fGHz-0.08, res.LOS.ds_med,  res.LOS.ds_p10,  res.LOS.ds_p90,  cL, 'o');
errbar(fGHz+0.08, res.NLOS.ds_med, res.NLOS.ds_p10, res.NLOS.ds_p90, cN, 's');
plot([fGHz(1) fGHz(end)], ds3gpp.LOS*[1 1],  ':', 'Color', cL);
plot([fGHz(1) fGHz(end)], ds3gpp.NLOS*[1 1], ':', 'Color', cN);
xlabel('Sub-band centre frequency [GHz]'); ylabel('RMS delay spread [ns]');
title('median, 10-90 %;  dotted: 3GPP InF median', 'FontSize', 7, 'FontWeight', 'normal');
grid on; box on; xlim([fGHz(1)-0.5 fGHz(end)+0.5]);
savefig_pdf(fig, cfg, 'fig_ds_vs_freq');

% ---- Fig: CDF of RMS delay spread at representative bands vs 3GPP log-normal ----
fig = newfig(W, H); hold on;
hh = []; leg = {};
for j = 1:numel(kk)
    k = kk(j); sh = (j-1)/(numel(kk)-1+eps);
    col = (1-sh)*[0.2 0.4 0.9] + sh*[0.9 0.2 0.1];
    x = sort(DS(ok & P.los==1, k) * 1e9); x = x(~isnan(x));
    if ~isempty(x)
        hh(end+1) = stairs(x, (1:numel(x))'/numel(x), '-', 'Color', col, 'LineWidth', 1); %#ok<AGROW>
        leg{end+1} = sprintf('LOS %.1f GHz', fGHz(k)); %#ok<AGROW>
    end
    x = sort(DS(ok & P.los==0, k) * 1e9); x = x(~isnan(x));
    if ~isempty(x)
        hh(end+1) = stairs(x, (1:numel(x))'/numel(x), '--', 'Color', col, 'LineWidth', 1); %#ok<AGROW>
        leg{end+1} = sprintf('NLOS %.1f GHz', fGHz(k)); %#ok<AGROW>
    end
end
xg = logspace(-1, log10(200), 400);                    % [ns]
Fl = 0.5*erfc(-(log10(xg*1e-9) - log10(ds3gpp.LOS*1e-9)) / (sqrt(2)*0.15));
hh(end+1) = plot(xg, Fl, 'k:', 'LineWidth', 1.2); leg{end+1} = '3GPP InF LOS';
if any(P.los==0)
    Fn = 0.5*erfc(-(log10(xg*1e-9) - log10(ds3gpp.NLOS*1e-9)) / (sqrt(2)*0.19));
    hh(end+1) = plot(xg, Fn, 'k-.', 'LineWidth', 1); leg{end+1} = '3GPP InF NLOS';
end
set(gca, 'XScale', 'log'); xlim([0.5 200]); ylim([0 1]);
xlabel('RMS delay spread [ns]'); ylabel('CDF');
legend(hh, leg, 'Location', 'southeast', 'FontSize', 6);
grid on; box on;
savefig_pdf(fig, cfg, 'fig_cdf_ds');

% ---- Fig: box plot of RMS delay spread per sub-band (no toolbox needed) ----
% box = 25-75 %, line = median, x = mean, whiskers = 10-90 %, dots = points outside
fig = newfig(W, H); hold on;
grp = {1, 'LOS', cL; 2, 'OLOS', [0.47 0.67 0.19]; 0, 'NLOS', cN};
present = find(arrayfun(@(g) any(ok & P.los == grp{g,1}), 1:size(grp,1)));
nG = numel(present);  bw = 0.6 / max(nG,1);           % box width per group [GHz]
hh = []; leg = {};
for gi = 1:nG
    g = present(gi); col = grp{g,3};
    off = (gi - (nG+1)/2) * bw * 1.15;                 % side-by-side offset
    for k = 1:nB
        x = DS(ok & P.los == grp{g,1}, k) * 1e9;  x = x(~isnan(x));
        if isempty(x), continue; end
        q  = [prctile_or_nan(x,10) prctile_or_nan(x,25) median(x) prctile_or_nan(x,75) prctile_or_nan(x,90)];
        xc = fGHz(k) + off;
        fill(xc + bw/2*[-1 1 1 -1], [q(2) q(2) q(4) q(4)], col, 'FaceAlpha', 0.25, 'EdgeColor', col);
        plot(xc + bw/2*[-1 1], [q(3) q(3)], '-', 'Color', col, 'LineWidth', 1.5);
        plot([xc xc], [q(1) q(2)], '-', 'Color', col);  plot([xc xc], [q(4) q(5)], '-', 'Color', col);
        plot(xc + bw/4*[-1 1], [q(1) q(1)], '-', 'Color', col);
        plot(xc + bw/4*[-1 1], [q(5) q(5)], '-', 'Color', col);
        plot(xc, mean(x), 'x', 'Color', 'k', 'MarkerSize', 5, 'LineWidth', 1);   % mean
        out = x(x < q(1) | x > q(5));
        if ~isempty(out), plot(xc*ones(size(out)), out, '.', 'Color', col, 'MarkerSize', 4); end
    end
    hh(end+1) = plot(NaN, NaN, 's', 'MarkerSize', 6, 'MarkerFaceColor', 0.75*[1 1 1] + 0.25*col, ...
                     'MarkerEdgeColor', col); leg{end+1} = grp{g,2}; %#ok<AGROW>
end
set(gca, 'YScale', 'log'); ylim([0.3 100]); xlim([fGHz(1)-0.6 fGHz(end)+0.6]);
set(gca, 'XTick', fGHz);
xlabel('Sub-band centre frequency [GHz]'); ylabel('RMS delay spread [ns]');
if ~isempty(hh), legend(hh, leg, 'Location', 'northwest', 'FontSize', 7); end
grid on; box on;
savefig_pdf(fig, cfg, 'fig_ds_boxplot');

% ---- Fig: RMS delay spread vs distance (where does multipath occur?) ----
fig = newfig(W, H); hold on;
hh = []; leg = {};
mk = {'o','s','^','d','v','p'};
for j = 1:numel(kk)
    k = kk(j); sh = (j-1)/(numel(kk)-1+eps);
    col = (1-sh)*[0.2 0.4 0.9] + sh*[0.9 0.2 0.1];
    m = ok & P.los == 1;
    hh(end+1) = plot(P.dist(m), DS(m,k)*1e9, mk{j}, 'Color', col, 'MarkerSize', 3); %#ok<AGROW>
    leg{end+1} = sprintf('%.1f GHz', fGHz(k)); %#ok<AGROW>
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', [1 2 5 10 20]);
xlim([0.9 max(P.dist)*1.05]); ylim([0.3 100]);
xlabel('Tx-Rx distance d [m]'); ylabel('RMS delay spread [ns]');
title('LOS positions', 'FontSize', 7, 'FontWeight', 'normal');
legend(hh, leg, 'Location', 'northwest', 'FontSize', 7); grid on; box on;
savefig_pdf(fig, cfg, 'fig_ds_vs_distance');

% ---- Fig: CDF of deviation from 3GPP InF-LOS path loss (LOS points) ----
fig = newfig(W, H); hold on;
hh = []; leg = {};
for j = 1:numel(kk)
    k = kk(j); sh = (j-1)/(numel(kk)-1+eps);
    col = (1-sh)*[0.2 0.4 0.9] + sh*[0.9 0.2 0.1];
    m = ok & P.los==1 & P.dist >= 1;
    dev = PL(m, k) - pl_3gpp_inf_los(P.dist(m), fGHz(k));
    x = sort(dev(~isnan(dev)));
    if ~isempty(x)
        hh(end+1) = stairs(x, (1:numel(x))'/numel(x), '-', 'Color', col, 'LineWidth', 1); %#ok<AGROW>
        leg{end+1} = sprintf('%.1f GHz (median %+.1f dB)', fGHz(k), median(x)); %#ok<AGROW>
    end
end
xg = linspace(-15, 15, 300);
hh(end+1) = plot(xg, 0.5*erfc(-xg/(sqrt(2)*4.3)), 'k:', 'LineWidth', 1.2);
leg{end+1} = '3GPP InF-LOS (\sigma = 4.3 dB)';
xlim([-15 15]); ylim([0 1]);
xlabel('PL_{meas} - PL_{3GPP InF-LOS} [dB]'); ylabel('CDF');
legend(hh, leg, 'Location', 'southeast', 'FontSize', 6);
grid on; box on;
savefig_pdf(fig, cfg, 'fig_cdf_pl_vs_3gpp');

% ---- Fig: floor plan (only if x/y available) ----
if isfield(P, 'x') && any(~isnan(P.x))
    fig = newfig(W, W * cfg.HALL_LWH(2)/cfg.HALL_LWH(1) + 0.4); hold on;
    Lh = cfg.HALL_LWH;
    plot([0 Lh(1) Lh(1) 0 0], [0 0 Lh(2) Lh(2) 0], 'k-', 'LineWidth', 1);
    hl = []; ll = {};
    if any(P.los==1)
        hl(end+1) = plot(P.x(P.los==1), P.y(P.los==1), 'o', 'Color', cL, 'MarkerSize', 2.5, 'MarkerFaceColor', cL);
        ll{end+1} = 'LOS Rx';
    end
    if any(P.los==0)
        hl(end+1) = plot(P.x(P.los==0), P.y(P.los==0), 's', 'Color', cN, 'MarkerSize', 2.5, 'MarkerFaceColor', cN);
        ll{end+1} = 'NLOS Rx';
    end
    if any(P.los==2)
        hl(end+1) = plot(P.x(P.los==2), P.y(P.los==2), '^', 'Color', [0.47 0.67 0.19], ...
                         'MarkerSize', 2.5, 'MarkerFaceColor', [0.47 0.67 0.19]);
        ll{end+1} = 'OLOS Rx';
    end
    if ~any(isnan(cfg.TX_XY))
        hl(end+1) = plot(cfg.TX_XY(1), cfg.TX_XY(2), 'p', 'Color', 'k', 'MarkerSize', 9, 'MarkerFaceColor', 'y');
        ll{end+1} = 'Tx';
    end
    axis equal; xlim([-0.5 Lh(1)+0.5]); ylim([-0.5 Lh(2)+0.5]);
    xlabel('x [m]'); ylabel('y [m]'); legend(hl, ll, 'FontSize', 7, 'Location', 'southwest');
    box on;
    savefig_pdf(fig, cfg, 'fig_floorplan');
end
end

function pl = pl_3gpp_inf_los(d, fGHz)
% 3GPP TR 38.901 Table 7.4.1-1, InF-LOS, d = d_3D [m], f in GHz, valid 1-600 m
pl = 31.84 + 21.50*log10(d) + 19.00*log10(fGHz);
end

function errbar(x, y, lo, hi, col, mk)
for k = 1:numel(x)
    if ~isnan(y(k)), plot([x(k) x(k)], [lo(k) hi(k)], '-', 'Color', col); end
end
plot(gx(x), gy(x, y), ['-' mk], 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 3);
end

function xo = gx(x)
% x values with NaN inserted where sub-bands are missing (breaks the line)
[xo, ~] = gap_insert(x, x);
end

function yo = gy(x, y)
[~, yo] = gap_insert(x, y);
end

function [xo, yo] = gap_insert(x, y)
x = x(:)'; y = y(:)'; xo = x(1); yo = y(1);
step = min(diff(x)); if isempty(step), step = 1; end
for i = 2:numel(x)
    if x(i) - x(i-1) > 1.5*step, xo(end+1) = NaN; yo(end+1) = NaN; end %#ok<AGROW>
    xo(end+1) = x(i); yo(end+1) = y(i); %#ok<AGROW>
end
end

function fig = newfig(w, h)
fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 w h], ...
             'PaperUnits', 'inches', 'PaperSize', [w h], 'PaperPosition', [0 0 w h]);
set(fig, 'DefaultAxesFontSize', 8);
end

function savefig_pdf(fig, cfg, name)
print(fig, fullfile(cfg.OUT_DIR, [name '.pdf']), '-dpdf');
print(fig, fullfile(cfg.OUT_DIR, [name '.png']), '-dpng', '-r200');
close(fig);
end