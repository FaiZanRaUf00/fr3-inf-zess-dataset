function check_1m(s2pFile, d)
% CHECK_1M  Quick in-lab check: is a measurement OK in ALL 1-GHz sub-bands?
%   check_1m('pt69_new.s2p', 1.0)
% Compares the gain-corrected path loss with free space at distance d [m]
% (use a short LOS distance, e.g. 1-2 m, horns aligned).
% PASS = within +-3 dB in every sub-band.

if nargin < 2, d = 1.0; end
gainFile = fullfile(fileparts(which(mfilename)), '..', 'data', 'antenna', 'gain.csv');
if ~exist(gainFile, 'file'), gainFile = 'gain.csv'; end

[f, S21] = read_s21(s2pFile);
G   = dlmread(gainFile, ',', 1, 0);
Gdb = 2 * interp1(G(:,1)*1e9, G(:,2), f, 'linear', 'extrap');       % Gt + Gr [dB]
H2  = abs(S21).^2 ./ 10.^(Gdb/10);                                    % antenna-free |H|^2

edges = 6e9:1e9:floor(max(f)/1e9)*1e9;
fprintf('\n%s, d = %.2f m: measured PL - free-space PL per sub-band\n', s2pFile, d);
allOK = true;
for k = 1:numel(edges)-1
    idx = f >= edges(k) & f < edges(k+1);
    if ~any(idx), continue; end
    fc  = mean(f(idx));
    pl  = -10*log10(mean(H2(idx)));
    dev = pl - 20*log10(4*pi*d*fc/299792458);
    ok  = abs(dev) <= 3;  allOK = allOK && ok;
    fprintf('  %5.1f GHz: %+6.1f dB   %s\n', fc/1e9, dev, ternary(ok, 'ok', '<-- FAIL'));
end
if allOK
    fprintf('==> PASS: all sub-bands within +-3 dB. Go ahead with the measurement.\n\n');
else
    fprintf('==> FAIL: fix the setup (cable / adapter / calibration) before measuring.\n\n');
end
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end

function [f, S21] = read_s21(fn)
fid = fopen(fn, 'r'); if fid < 0, error('Cannot open %s', fn); end
unit = 1e9; fmt = 'MA'; data = [];
while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    ln = strtrim(ln);
    if isempty(ln) || ln(1) == '!', continue; end
    if ln(1) == '#'
        tok = upper(strsplit(strtrim(ln(2:end))));
        for t = tok
            switch t{1}
                case 'HZ', unit = 1; case 'KHZ', unit = 1e3; case 'MHZ', unit = 1e6; case 'GHZ', unit = 1e9;
                case {'RI','MA','DB'}, fmt = t{1};
            end
        end
        continue
    end
    c = strfind(ln, '!'); if ~isempty(c), ln = ln(1:c(1)-1); end
    v = sscanf(ln, '%f')';
    if numel(v) >= 5, data = [data; v(1:5)]; end %#ok<AGROW>
end
fclose(fid);
f = data(:,1) * unit; a = data(:,4); b = data(:,5);
switch fmt
    case 'RI', S21 = a + 1j*b;
    case 'MA', S21 = a .* exp(1j*b*pi/180);
    case 'DB', S21 = 10.^(a/20) .* exp(1j*b*pi/180);
end
end
