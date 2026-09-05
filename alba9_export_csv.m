function alba9_export_csv(matfile,csvfile,ixfun,maxrows)
%ALBA9_EXPORT_CSV Decimated CSV view of one case, for reading without MATLAB.
%
%   The .mat file is the FULL record. This CSV is a decimated VIEW of it, so
%   that a reviewer without MATLAB can still read the trajectory.
%
%   Decimation rules, chosen so the closure is never blurred:
%     * duplicate timestamps are dropped, keeping the last row on each side of
%       the valve discontinuity, so both one-sided limits survive;
%     * every row within +-0.05 s of the closing instant is kept in full;
%     * the remainder is thinned uniformly in INDEX to at most maxrows.
%   No interpolation and no smoothing is applied anywhere.
if nargin < 4 || isempty(maxrows), maxrows = 6000; end
S = load(matfile);
ix = ixfun();
t = S.t; y = S.y;

vc = y(:,ix.valve_closed) > 0.5;
keep = false(numel(t),1);
for side = [false true]
    idx = find(vc == side);
    if isempty(idx), continue; end
    sub = idx([diff(t(idx)) > 0; true]);
    keep(sub) = true;
end

tc = NaN;
if any(vc), tc = t(find(vc,1)); end
near = ~isnan(tc) & abs(t-tc) <= 0.05;

sel = find(keep);
if numel(sel) > maxrows
    step = ceil(numel(sel)/maxrows);
    thin = false(numel(t),1);
    thin(sel(1:step:end)) = true;
    thin(sel(end)) = true;
    keep = thin;
end
keep = keep | (near & [diff(t) > 0; true]);

f = fieldnames(ix);
cols = zeros(numel(f),1);
for i = 1:numel(f), cols(i) = ix.(f{i}); end
T = array2table([t(keep), y(keep,cols)],'VariableNames',[{'time_s'}; f]);
writetable(T,csvfile);
fprintf('  %s : %d of %d rows\n',csvfile,sum(keep),numel(t));
end
