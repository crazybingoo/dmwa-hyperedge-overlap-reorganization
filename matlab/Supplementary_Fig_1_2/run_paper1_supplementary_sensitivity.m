%% Paper 1 supplementary sensitivity analyses (2026-07-14)
% Full-window reconstruction checks for Supplementary Figs. 1 and 2.
% Analyses are performed at seizure-window level and later aggregated first
% within seizure and then within patient, matching the manuscript hierarchy.

clear; clc;
rng(20260714, 'twister');

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot), projectRoot = pwd; end
analysisRoot = getenv('PAPER1_ANALYSIS_CACHE_ROOT');
if isempty(analysisRoot)
    analysisRoot = fullfile(projectRoot, 'private_cache', 'supplementary_sensitivity');
end
cacheDir = fullfile(analysisRoot, 'cache');
outDir = fullfile(projectRoot, 'stats', 'supplementary_sensitivity');
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
if ~exist(outDir, 'dir'), mkdir(outDir); end

meta = readtable(fullfile(projectRoot, 'data_tables', 'paper1_case_metadata.csv'), 'TextType', 'string');
thr = readtable(fullfile(projectRoot, 'data_tables', 'paper2_threshold_elbow_case_summary.csv'), 'TextType', 'string');

fs = 1024;
windowSamples = 3 * fs;
fixedDensity = 0.216;
ratioGrid = (0.70:0.05:0.95)';
variantNames = ["PLV baseline", "PLV pre-ictal threshold", "PLV threshold -10%", ...
    "PLV threshold +10%", "PLV fixed density", "PLI fixed density", "wPLI fixed density"];

metricRows = cell(120000, 14);
ratioRows = cell(120000, 12);
caseRows = cell(height(meta), 8);
metricCursor = 0;
ratioCursor = 0;

for c = 1:height(meta)
    caseId = meta.case_id(c);
    patientId = infer_patient_id(caseId);
    match = thr(thr.case_id == caseId, :);
    if height(match) ~= 1
        error('Expected one baseline threshold for %s.', caseId);
    end
    baselineThreshold = match.selected_threshold(1);
    fprintf('[%02d/%02d] %s\n', c, height(meta), caseId);

    cachePath = fullfile(cacheDir, caseId + "_sync_edges.mat");
    if exist(cachePath, 'file')
        C = load(cachePath, 'plvEdges', 'pliEdges', 'wpliEdges', 'nChannel', 'nWindow');
        plvEdges = C.plvEdges;
        pliEdges = C.pliEdges;
        wpliEdges = C.wpliEdges;
        nChannel = C.nChannel;
        nWindow = C.nWindow;
    else
        D = load(meta.file_path(c));
        if ~isfield(D, 'X1')
            error('Missing X1 in %s.', meta.file_path(c));
        end
        X = double(D.X1);
        nChannel = size(X, 1);
        nWindow = floor(size(X, 2) / fs) - 3 + 1;
        [pairI, pairJ] = find(triu(true(nChannel), 1));
        nPairs = numel(pairI);
        plvEdges = zeros(nWindow, nPairs, 'single');
        pliEdges = zeros(nWindow, nPairs, 'single');
        wpliEdges = zeros(nWindow, nPairs, 'single');
        for w = 1:nWindow
            idx = fs * (w - 1) + (1:windowSamples);
            phase = angle(hilbert(X(:, idx).')).';
            dphi = phase(pairI, :) - phase(pairJ, :);
            sinDiff = sin(dphi);
            plvEdges(w, :) = single(abs(mean(exp(1i * dphi), 2))).';
            pliEdges(w, :) = single(abs(mean(sign(sinDiff), 2))).';
            wpliEdges(w, :) = single((abs(mean(sinDiff, 2)) ./ (mean(abs(sinDiff), 2) + eps))).';
            if mod(w, 100) == 0
                fprintf('  synchrony %d/%d\n', w, nWindow);
            end
        end
        save(cachePath, 'plvEdges', 'pliEdges', 'wpliEdges', 'nChannel', 'nWindow', '-v7.3');
        clear X D phase dphi sinDiff;
    end

    stage = assign_stage_labels(nWindow, meta.pre_start(c):meta.pre_end(c), ...
        meta.ictal_start(c):meta.ictal_end(c), meta.post_start(c):meta.post_end(c));
    preMask = stage == "pre-ictal";
    preThreshold = threshold_density_elbow(double(plvEdges(preMask, :)));
    caseRows(c, :) = {char(caseId), char(patientId), nChannel, nWindow, baselineThreshold, ...
        preThreshold, fixedDensity, mean(preMask)};

    [pairI, pairJ] = find(triu(true(nChannel), 1));
    for w = 1:nWindow
        values = {double(plvEdges(w, :)), double(plvEdges(w, :)), double(plvEdges(w, :)), ...
            double(plvEdges(w, :)), double(plvEdges(w, :)), double(pliEdges(w, :)), double(wpliEdges(w, :))};
        rules = [baselineThreshold, preThreshold, max(0, baselineThreshold * 0.90), ...
            min(1, baselineThreshold * 1.10), NaN, NaN, NaN];
        baselineW = [];
        for v = 1:numel(variantNames)
            if v <= 4
                A = adjacency_threshold(values{v}, pairI, pairJ, nChannel, rules(v));
                ruleValue = rules(v);
            else
                A = adjacency_fixed_density(values{v}, pairI, pairJ, nChannel, fixedDensity);
                ruleValue = fixedDensity;
            end
            HE = hyperedges_2_3(A);
            W = build_dmwa_from_hyperedges(HE);
            M = directed_dmwa_metrics(W);
            metricCursor = metricCursor + 1;
            metricRows(metricCursor, :) = {char(variantNames(v)), char(caseId), char(patientId), w, ...
                char(stage(w)), nChannel, ruleValue, nnz(triu(A, 1)) / max(1, nchoosek(nChannel, 2)), ...
                numel(HE), M.S, M.source_sink_gini_gap, M.sink_minus_source_entropy_gap, ...
                M.sink_minus_source_pr_gap, mean(values{v}, 'omitnan')};
            if v == 1
                baselineW = W;
            end
        end

        if ~isempty(baselineW)
            triads = common_triads(size(baselineW, 1), 5000, 20260714 + c * 100000 + w);
            lut = triad13_lut();
            for rr = 1:numel(ratioGrid)
                A_dir = adjacency_from_W(baselineW, ratioGrid(rr));
                G = motif_group_fractions(A_dir, triads, lut);
                ratioCursor = ratioCursor + 1;
                ratioRows(ratioCursor, :) = {char(caseId), char(patientId), w, char(stage(w)), ...
                    ratioGrid(rr), size(baselineW, 1), size(triads, 1), G.connected_count, ...
                    G.simple_fraction, G.reciprocal_fraction, G.dense_fraction, G.complex_fraction};
            end
        end
        if mod(w, 50) == 0
            fprintf('  reconstruction %d/%d\n', w, nWindow);
        end
    end

    metricT = cell2table(metricRows(1:metricCursor, :), 'VariableNames', ...
        {'variant','case_id','patient_id','window_idx','stage','n_channels','rule_value', ...
        'pair_density','n_hyperedges','S','source_sink_gini_gap', ...
        'sink_minus_source_entropy_gap','sink_minus_source_pr_gap','mean_synchrony'});
    ratioT = cell2table(ratioRows(1:ratioCursor, :), 'VariableNames', ...
        {'case_id','patient_id','window_idx','stage','reciprocity_ratio','n_dmwa_nodes', ...
        'triads_evaluated','connected_count','simple_fraction','reciprocal_fraction', ...
        'dense_fraction','complex_fraction'});
    metricT = normalize_numeric_columns(metricT, [4, 6:14]);
    ratioT = normalize_numeric_columns(ratioT, [3, 5:12]);
    writetable(metricT, fullfile(outDir, 'supp_fig1_synchrony_sensitivity_window.csv'));
    writetable(ratioT, fullfile(outDir, 'supp_fig2_reciprocity_ratio_window.csv'));
end

caseT = cell2table(caseRows, 'VariableNames', {'case_id','patient_id','n_channels','n_windows', ...
    'baseline_threshold','preictal_only_threshold','fixed_density','preictal_window_fraction'});
caseT = normalize_numeric_columns(caseT, 3:8);
writetable(caseT, fullfile(outDir, 'supp_fig1_case_rules.csv'));

fprintf('Completed supplementary sensitivity extraction.\n');

function T = normalize_numeric_columns(T, idx)
for i = idx
    name = T.Properties.VariableNames{i};
    if iscell(T.(name))
        T.(name) = cell2mat(T.(name));
    end
end
end

function patientId = infer_patient_id(caseId)
s = char(caseId);
k = strfind(s, '_cut');
if isempty(k), k = strfind(s, '_process'); end
if isempty(k), patientId = string(s); else, patientId = string(s(1:k(1)-1)); end
end

function stage = assign_stage_labels(nWindow, preIdx, ictalIdx, postIdx)
stage = repmat("unassigned", nWindow, 1);
preIdx = sanitize_idx(preIdx, nWindow);
ictalIdx = setdiff(sanitize_idx(ictalIdx, nWindow), preIdx, 'stable');
postIdx = setdiff(sanitize_idx(postIdx, nWindow), union(preIdx, ictalIdx), 'stable');
splitPts = round(linspace(0, numel(ictalIdx), 4));
stage(preIdx) = "pre-ictal";
stage(ictalIdx(1:splitPts(2))) = "early";
stage(ictalIdx(splitPts(2)+1:splitPts(3))) = "mid";
stage(ictalIdx(splitPts(3)+1:end)) = "late";
stage(postIdx) = "post-ictal";
end

function idx = sanitize_idx(idx, nWindow)
idx = unique(idx(idx >= 1 & idx <= nWindow), 'stable');
end

function threshold = threshold_density_elbow(edgeValues)
grid = (0:0.0025:1)';
density = zeros(size(grid));
for i = 1:numel(grid)
    density(i) = mean(edgeValues(:) >= grid(i));
end
x = (grid - grid(1)) / (grid(end) - grid(1) + eps);
y = (density - density(end)) / (density(1) - density(end) + eps);
p1 = [x(1), y(1)]; p2 = [x(end), y(end)];
dist = abs((p2(2)-p1(2))*x - (p2(1)-p1(1))*y + p2(1)*p1(2)-p2(2)*p1(1)) ...
    / (hypot(p2(2)-p1(2), p2(1)-p1(1)) + eps);
[~, k] = max(dist);
threshold = grid(k);
end

function A = adjacency_threshold(values, pairI, pairJ, n, threshold)
A = false(n);
keep = values(:) >= threshold;
A(sub2ind([n n], pairI(keep), pairJ(keep))) = true;
A = A | A.';
A(1:n+1:end) = false;
end

function A = adjacency_fixed_density(values, pairI, pairJ, n, density)
A = false(n);
k = max(0, min(numel(values), round(density * numel(values))));
if k > 0
    [~, ord] = sort(values(:), 'descend');
    keep = ord(1:k);
    A(sub2ind([n n], pairI(keep), pairJ(keep))) = true;
    A = A | A.';
end
A(1:n+1:end) = false;
end

function HE = hyperedges_2_3(A)
n = size(A, 1);
[row, col] = find(triu(A, 1));
HE2 = arrayfun(@(i) [row(i), col(i)], 1:numel(row), 'UniformOutput', false)';
tri = zeros(0, 3);
for i = 1:n-2
    neigh = find(A(i, i+1:end)) + i;
    for jj = 1:numel(neigh)-1
        j = neigh(jj);
        ks = neigh(jj+1:end);
        ks = ks(A(j, ks));
        if ~isempty(ks)
            tri = [tri; [repmat(i, numel(ks), 1), repmat(j, numel(ks), 1), ks(:)]]; %#ok<AGROW>
        end
    end
end
if isempty(tri)
    HE3 = cell(0, 1);
else
    tri = unique(tri, 'rows', 'stable');
    HE3 = mat2cell(tri, ones(size(tri, 1), 1), 3);
end
HE = [HE2(:); HE3(:)];
end

function W = build_dmwa_from_hyperedges(HE)
m = numel(HE);
W = zeros(m);
if m < 2, return; end
sizes = cellfun(@numel, HE);
maxNode = max(cellfun(@max, HE));
rowIdx = zeros(sum(sizes), 1); colIdx = zeros(sum(sizes), 1); cursor = 0;
for i = 1:m
    nodes = HE{i}(:); k = numel(nodes);
    rowIdx(cursor + (1:k)) = i; colIdx(cursor + (1:k)) = nodes; cursor = cursor + k;
end
H = sparse(rowIdx, colIdx, 1, m, maxNode);
O = full(H * H.'); O(1:m+1:end) = 0;
D = sum(O > 0, 2);
Si = repmat(sizes(:), 1, m); Sj = repmat(sizes(:)', m, 1);
Di = repmat(D, 1, m); Dj = repmat(D', m, 1);
W = O ./ Si;
degreeSum = Di + Dj;
degreeBias = Dj ./ degreeSum; degreeBias(degreeSum == 0) = 0.5;
sameOrder = (Si == Sj) & (O > 0);
W(sameOrder) = W(sameOrder) .* degreeBias(sameOrder);
W(O == 0) = 0; W(1:m+1:end) = 0;
end

function M = directed_dmwa_metrics(W)
n = size(W, 1);
if n < 2
    M.S = 0; M.source_sink_gini_gap = 0; M.sink_minus_source_entropy_gap = 0; M.sink_minus_source_pr_gap = 0; return;
end
mask = triu(true(n), 1); Wij = W(mask); WT = W.'; Wji = WT(mask);
den = sum(max(Wij, Wji));
if den > 0, M.S = sum(min(Wij, Wji)) / den; else, M.S = 0; end
sout = sum(W, 2); sinv = sum(W, 1)'; q = (sout - sinv) ./ (sout + sinv + eps);
source = max(q, 0); sink = max(-q, 0);
M.source_sink_gini_gap = gini_nonnegative(source) - gini_nonnegative(sink);
M.sink_minus_source_entropy_gap = entropy_norm(sink) - entropy_norm(source);
M.sink_minus_source_pr_gap = participation_norm(sink) - participation_norm(source);
end

function g = gini_nonnegative(x)
x = max(0, x(:)); n = numel(x);
if n == 0 || sum(x) <= 0, g = 0; return; end
x = sort(x); g = 2 * sum((1:n)' .* x) / (n * sum(x)) - (n + 1) / n;
end

function h = entropy_norm(x)
x = max(0, x(:)); n = numel(x); s = sum(x);
if n <= 1 || s <= 0, h = 0; return; end
p = x / s; p = p(p > 0); h = -sum(p .* log(p)) / log(n);
end

function pr = participation_norm(x)
x = max(0, x(:)); n = numel(x);
if n == 0 || sum(x) <= 0 || sum(x.^2) <= 0, pr = 0; return; end
pr = (sum(x)^2 / sum(x.^2)) / n;
end

function A = adjacency_from_W(W, reciprocityRatio)
n = size(W, 1); A = false(n);
for i = 1:n-1
    for j = i+1:n
        wij = W(i,j); wji = W(j,i); vmax = max(wij,wji);
        if vmax <= 1e-12, continue; end
        if min(wij,wji) / (vmax + eps) >= reciprocityRatio
            A(i,j) = true; A(j,i) = true;
        elseif wij > wji
            A(i,j) = true;
        elseif wji > wij
            A(j,i) = true;
        else
            A(i,j) = true; A(j,i) = true;
        end
    end
end
end

function triads = common_triads(n, maxTriads, seed)
if n < 3, triads = zeros(0,3); return; end
total = nchoosek(n,3);
if total <= maxTriads
    triads = nchoosek(1:n,3);
    return;
end
rng(seed, 'twister'); triads = zeros(0,3);
while size(triads,1) < maxTriads
    need = maxTriads - size(triads,1);
    cand = sort(randi(n, max(need*2,1000), 3), 2);
    cand = cand(cand(:,1) ~= cand(:,2) & cand(:,2) ~= cand(:,3), :);
    triads = unique([triads; cand], 'rows', 'stable'); %#ok<AGROW>
end
triads = triads(1:maxTriads,:);
end

function G = motif_group_fractions(A, triads, lut)
G.connected_count = 0; G.simple_fraction = NaN; G.reciprocal_fraction = NaN; G.dense_fraction = NaN; G.complex_fraction = NaN;
if isempty(triads), return; end
n = size(A,1); i=triads(:,1); j=triads(:,2); k=triads(:,3);
codes = A(sub2ind([n n],i,j)) + 2*A(sub2ind([n n],i,k)) + ...
    4*A(sub2ind([n n],j,i)) + 8*A(sub2ind([n n],j,k)) + ...
    16*A(sub2ind([n n],k,i)) + 32*A(sub2ind([n n],k,j));
idx = lut(double(codes)+1); idx = idx(idx>0); G.connected_count = numel(idx);
if isempty(idx), return; end
G.simple_fraction = mean(ismember(idx,[1 2 3]));
G.reciprocal_fraction = mean(ismember(idx,[4 5 8 9 10 11 12 13]));
G.dense_fraction = mean(ismember(idx,[9 10 11 12 13]));
G.complex_fraction = mean(ismember(idx,4:13));
end

function lut = triad13_lut()
T = cell(13,1);
T{1}=logical([0 1 1;0 0 0;0 0 0]); T{2}=logical([0 0 0;1 0 0;1 0 0]);
T{3}=logical([0 1 0;0 0 1;0 0 0]); T{4}=logical([0 1 0;1 0 0;1 0 0]);
T{5}=logical([0 1 1;1 0 0;0 0 0]); T{6}=logical([0 1 1;0 0 1;0 0 0]);
T{7}=logical([0 1 0;0 0 1;1 0 0]); T{8}=logical([0 1 1;1 0 0;1 0 0]);
T{9}=logical([0 1 1;1 0 1;0 0 0]); T{10}=logical([0 1 0;1 0 0;1 1 0]);
T{11}=logical([0 1 1;1 0 0;0 1 0]); T{12}=logical([0 1 1;1 0 1;1 0 0]);
T{13}=logical([0 1 1;1 0 1;1 1 0]);
lut = zeros(64,1); P = perms(1:3); weights=[1;2;4;8;16;32];
for m=1:13
    for p=1:size(P,1)
        B=T{m}(P(p,:),P(p,:)); bits=[B(1,2),B(1,3),B(2,1),B(2,3),B(3,1),B(3,2)];
        lut(bits*weights+1)=m;
    end
end
end
