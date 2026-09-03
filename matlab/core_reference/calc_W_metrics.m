%#ok<*AGROW,*SAGROW>
function metrics = calc_W_metrics(W)
%CALC_W_METRICS Compute DMWA role-reciprocity metrics for one window.
%
% Input
%   W : Directed multi-order weighted adjacency matrix (DMWA), with one
%       row/column per hyperedge.
%
% Output
%   metrics : Struct containing global reciprocal-balance S, source/sink
%       strength summaries, entropy/participation summaries, and the
%       hyperedge-level role score q = (out - in) / (out + in).

W = double(W);
W(~isfinite(W)) = 0;
W(W < 0) = 0;
n = size(W, 1);

if n == 0
    metrics = empty_metrics();
    return;
end
if size(W, 2) ~= n
    error('W must be a square matrix.');
end
W(1:n+1:end) = 0;

if n >= 2
    mask = triu(true(n), 1);
    Wt = W.';
    wij = W(mask);
    wji = Wt(mask);
    denominator = sum(max(wij, wji));
    if denominator > 0
        metrics.S = sum(min(wij, wji)) / denominator;
    else
        metrics.S = 0;
    end
else
    metrics.S = 0;
end

outStrength = sum(W, 2);
inStrength = sum(W, 1).';
q = (outStrength - inStrength) ./ (outStrength + inStrength + eps);
sourceLoad = max(q, 0);
sinkLoad = max(-q, 0);

metrics.n_hyperedges = n;
metrics.total_weight = sum(W, 'all');
metrics.mean_nonzero_weight = mean(W(W > 0), 'omitnan');
if isnan(metrics.mean_nonzero_weight)
    metrics.mean_nonzero_weight = 0;
end
metrics.out_strength = outStrength;
metrics.in_strength = inStrength;
metrics.q = q;
metrics.source_fraction = mean(q > 0);
metrics.sink_fraction = mean(q < 0);
metrics.mean_abs_q = mean(abs(q), 'omitnan');
metrics.source_gini = gini_nonnegative(sourceLoad);
metrics.sink_gini = gini_nonnegative(sinkLoad);
metrics.source_entropy = entropy_norm(sourceLoad);
metrics.sink_entropy = entropy_norm(sinkLoad);
metrics.source_participation = participation_norm(sourceLoad);
metrics.sink_participation = participation_norm(sinkLoad);
metrics.source_sink_gini_gap = metrics.source_gini - metrics.sink_gini;
metrics.sink_minus_source_entropy_gap = metrics.sink_entropy - metrics.source_entropy;
metrics.sink_minus_source_pr_gap = metrics.sink_participation - metrics.source_participation;
end

function metrics = empty_metrics()
metrics = struct( ...
    'n_hyperedges', 0, ...
    'S', 0, ...
    'total_weight', 0, ...
    'mean_nonzero_weight', 0, ...
    'out_strength', [], ...
    'in_strength', [], ...
    'q', [], ...
    'source_fraction', 0, ...
    'sink_fraction', 0, ...
    'mean_abs_q', 0, ...
    'source_gini', 0, ...
    'sink_gini', 0, ...
    'source_entropy', 0, ...
    'sink_entropy', 0, ...
    'source_participation', 0, ...
    'sink_participation', 0, ...
    'source_sink_gini_gap', 0, ...
    'sink_minus_source_entropy_gap', 0, ...
    'sink_minus_source_pr_gap', 0);
end

function g = gini_nonnegative(x)
x = x(:);
x(~isfinite(x)) = 0;
x(x < 0) = 0;
n = numel(x);
if n == 0 || sum(x) <= 0
    g = 0;
    return;
end
x = sort(x, 'ascend');
g = (2 * sum((1:n)' .* x) / (n * sum(x))) - (n + 1) / n;
end

function h = entropy_norm(x)
x = x(:);
x(~isfinite(x)) = 0;
x(x < 0) = 0;
s = sum(x);
n = numel(x);
if n <= 1 || s <= 0
    h = 0;
    return;
end
p = x / s;
p = p(p > 0);
h = -sum(p .* log(p)) / log(n);
end

function pr = participation_norm(x)
x = x(:);
x(~isfinite(x)) = 0;
x(x < 0) = 0;
n = numel(x);
if n == 0 || sum(x) <= 0 || sum(x.^2) <= 0
    pr = 0;
    return;
end
pr = (sum(x)^2 / sum(x.^2)) / n;
end
