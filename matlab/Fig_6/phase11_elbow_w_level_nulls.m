%#ok<*AGROW,*SAGROW>
%% Phase 11: full-window elbow-threshold W-level nulls for Paper 1
%
% This script upgrades the previous sampled W-level diagnostics by applying
% W-level null models to every eligible elbow-threshold window.
%
% Null models:
%   1. weight_shuffle: preserves the nonzero directed support of W and
%      shuffles the marginal nonzero weight distribution.
%   2. orientation_pair_swap: randomly swaps W_ij and W_ji for each dyadic
%      hyperedge pair, preserving pair-level reciprocal weights but
%      disrupting consistent source/sink orientation.
%
% Optional caller variables:
%   PHASE11_NITER      number of null iterations, default 100.
%   PHASE11_CASE_LIMIT number of cases for smoke testing, default all cases.
%   PHASE11_CASE_START first case row in metadata, default 1.
%   PHASE11_CASE_END   last case row in metadata, default all cases.
%   PHASE11_OUTPUT_SUFFIX suffix inserted before .csv output extension.
%   PHASE11_ITER_OFFSET value added to iteration IDs, default 0.

clearvars -except PHASE11_NITER PHASE11_CASE_LIMIT PHASE11_CASE_START PHASE11_CASE_END PHASE11_OUTPUT_SUFFIX PHASE11_ITER_OFFSET; clc;
rng(20260617);

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot)
    projectRoot = pwd;
end
dataDir = fullfile(projectRoot, 'data_tables');
statsDir = fullfile(projectRoot, 'stats');
noteDir = fullfile(projectRoot, 'notes');
if ~exist(statsDir, 'dir'), mkdir(statsDir); end
if ~exist(noteDir, 'dir'), mkdir(noteDir); end

if ~exist('PHASE11_NITER', 'var') || isempty(PHASE11_NITER)
    nIter = 100;
else
    nIter = PHASE11_NITER;
end
if ~exist('PHASE11_CASE_LIMIT', 'var') || isempty(PHASE11_CASE_LIMIT)
    caseLimit = Inf;
else
    caseLimit = PHASE11_CASE_LIMIT;
end
if ~exist('PHASE11_CASE_START', 'var') || isempty(PHASE11_CASE_START)
    caseStart = 1;
else
    caseStart = PHASE11_CASE_START;
end
if ~exist('PHASE11_CASE_END', 'var') || isempty(PHASE11_CASE_END)
    caseEnd = Inf;
else
    caseEnd = PHASE11_CASE_END;
end
if ~exist('PHASE11_OUTPUT_SUFFIX', 'var') || isempty(PHASE11_OUTPUT_SUFFIX)
    outputSuffix = '';
else
    outputSuffix = char(PHASE11_OUTPUT_SUFFIX);
end
if ~exist('PHASE11_ITER_OFFSET', 'var') || isempty(PHASE11_ITER_OFFSET)
    iterOffset = 0;
else
    iterOffset = PHASE11_ITER_OFFSET;
end

metadataPath = fullfile(dataDir, 'paper1_case_metadata.csv');
thresholdPath = fullfile(dataDir, 'paper2_threshold_elbow_case_summary.csv');
empOut = fullfile(statsDir, ['phase11_elbow_w_null_empirical_case_stage' outputSuffix '.csv']);
nullOut = fullfile(statsDir, ['phase11_elbow_w_null_case_stage' outputSuffix '.csv']);
noteOut = fullfile(noteDir, ['phase11_elbow_w_level_null_extraction' outputSuffix '.md']);

caseMeta = readtable(metadataPath, 'TextType', 'string');
thresholdSummary = readtable(thresholdPath, 'TextType', 'string');
caseStart = max(1, caseStart);
caseEnd = min(height(caseMeta), min(caseEnd, caseStart + caseLimit - 1));
caseMeta = caseMeta(caseStart:caseEnd, :);

stageNames = ["pre-ictal","early","mid","late","post-ictal"];
metricNames = ["S","source_sink_gini_gap","sink_minus_source_entropy_gap","sink_minus_source_pr_gap"];
nullModels = ["weight_shuffle","orientation_pair_swap"];

fs = 1024;
winLenSec = 3;
nCase = height(caseMeta);
nStage = numel(stageNames);
nMetric = numel(metricNames);
nModel = numel(nullModels);

empSums = zeros(nCase, nStage, nMetric);
empCounts = zeros(nCase, nStage);
nullSums = zeros(nModel, nIter, nCase, nStage, nMetric);
nullCounts = zeros(nModel, nIter, nCase, nStage);

fprintf('Phase 11 full-window W-level nulls: cases %d-%d, %d cases, %d iterations.\n', ...
    caseStart, caseEnd, nCase, nIter);

for c = 1:nCase
    caseId = string(caseMeta.case_id(c));
    match = thresholdSummary(thresholdSummary.case_id == caseId, :);
    if height(match) ~= 1
        error('Cannot find one selected threshold for %s', caseId);
    end
    threshold = match.selected_threshold(1);
    filePath = resolve_local_path(char(caseMeta.file_path(c)), projectRoot);
    fprintf('[%02d/%02d] %s, threshold %.4f\n', c, nCase, caseId, threshold);

    S0 = load(filePath, 'X1');
    X1 = S0.X1;
    nSample = size(X1, 2);
    nWindow = floor(nSample / fs) - winLenSec + 1;
    stageLabels = assign_stage_labels(nWindow, ...
        caseMeta.pre_start(c):caseMeta.pre_end(c), ...
        caseMeta.ictal_start(c):caseMeta.ictal_end(c), ...
        caseMeta.post_start(c):caseMeta.post_end(c));

    for w = 1:nWindow
        stage = string(stageLabels{w});
        stageIdx = find(stageNames == stage, 1);
        if isempty(stageIdx)
            continue;
        end

        startPoint = fs * (w - 1) + 1;
        endPoint = fs * (w - 1 + winLenSec);
        datanew = X1(:, startPoint:endPoint);

        hyperEdges = extract_hyperedges_at_threshold(datanew, threshold);
        W = build_dmwa_from_hyperedges(hyperEdges);

        pairRep = w_pair_representation(W);
        empiricalMetrics = pair_directed_metrics(pairRep.n, pairRep.i, pairRep.j, pairRep.wU, pairRep.wL);
        empSums(c, stageIdx, :) = squeeze(empSums(c, stageIdx, :)) + empiricalMetrics(:);
        empCounts(c, stageIdx) = empCounts(c, stageIdx) + 1;

        for mi = 1:nModel
            model = nullModels(mi);
            for it = 1:nIter
                if model == "weight_shuffle"
                    nullMetrics = weight_shuffle_metrics(pairRep);
                else
                    nullMetrics = orientation_pair_swap_metrics(pairRep);
                end
                nullSums(mi, it, c, stageIdx, :) = squeeze(nullSums(mi, it, c, stageIdx, :)) + nullMetrics(:);
                nullCounts(mi, it, c, stageIdx) = nullCounts(mi, it, c, stageIdx) + 1;
            end
        end
    end
end

empTable = build_empirical_table(caseMeta, stageNames, metricNames, empSums, empCounts);
nullTable = build_null_table(caseMeta, stageNames, metricNames, nullModels, nullSums, nullCounts, iterOffset);
writetable(empTable, empOut);
writetable(nullTable, nullOut);

fid = fopen(noteOut, 'w');
fprintf(fid, '# Phase 11 full-window elbow-threshold W-level null extraction\n\n');
fprintf(fid, 'Date: 2026-06-17\n\n');
fprintf(fid, '- Case rows processed: %d-%d\n', caseStart, caseEnd);
fprintf(fid, '- Cases processed: %d\n', nCase);
fprintf(fid, '- Null iterations per model: %d\n', nIter);
fprintf(fid, '- Iteration offset: %d\n', iterOffset);
fprintf(fid, '- Threshold source: `%s`\n', thresholdPath);
fprintf(fid, '- Case metadata: `%s`\n', metadataPath);
fprintf(fid, '- Empirical case-stage output: `%s`\n', empOut);
fprintf(fid, '- Null case-stage output: `%s`\n', nullOut);
fprintf(fid, '- Rule: one seizure-specific elbow threshold applied to all windows from that seizure.\n');
fprintf(fid, '- Null model `weight_shuffle`: preserves W support and shuffles nonzero weights.\n');
fprintf(fid, '- Null model `orientation_pair_swap`: swaps W_ij/W_ji within hyperedge pairs.\n');
fclose(fid);

fprintf('Wrote empirical case-stage table: %s\n', empOut);
fprintf('Wrote null case-stage table: %s\n', nullOut);

%% Local functions
function outPath = resolve_local_path(rawPath, projectRoot)
rawPath = char(rawPath);
if isempty(rawPath)
    error('Encountered an empty path in case metadata.');
end
isAbs = startsWith(rawPath, filesep) || ~isempty(regexp(rawPath, '^[A-Za-z]:[\\/]', 'once'));
if isAbs
    outPath = rawPath;
else
    outPath = fullfile(projectRoot, rawPath);
end
end

function T = build_empirical_table(caseMeta, stageNames, metricNames, sums, counts)
rows = {};
rowIdx = 0;
for c = 1:height(caseMeta)
    caseId = char(caseMeta.case_id(c));
    patientId = infer_patient_id(caseId);
    for si = 1:numel(stageNames)
        nWin = counts(c, si);
        vals = nan(1, numel(metricNames));
        if nWin > 0
            vals = squeeze(sums(c, si, :))' ./ nWin;
        end
        rowIdx = rowIdx + 1;
        rows(rowIdx, :) = [{caseId, patientId, char(stageNames(si)), nWin}, num2cell(vals)];
    end
end
T = cell2table(rows, 'VariableNames', [{'case_id','patient_id','phase5','n_windows'}, cellstr(metricNames)]);
end

function T = build_null_table(caseMeta, stageNames, metricNames, nullModels, sums, counts, iterOffset)
rows = {};
rowIdx = 0;
for mi = 1:numel(nullModels)
    for it = 1:size(sums, 2)
        for c = 1:height(caseMeta)
            caseId = char(caseMeta.case_id(c));
            patientId = infer_patient_id(caseId);
            for si = 1:numel(stageNames)
                nWin = counts(mi, it, c, si);
                vals = nan(1, numel(metricNames));
                if nWin > 0
                    vals = squeeze(sums(mi, it, c, si, :))' ./ nWin;
                end
                rowIdx = rowIdx + 1;
                rows(rowIdx, :) = [{char(nullModels(mi)), it + iterOffset, caseId, patientId, char(stageNames(si)), nWin}, num2cell(vals)];
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', [{'null_model','iter','case_id','patient_id','phase5','n_windows'}, cellstr(metricNames)]);
end

function patient_id = infer_patient_id(case_id)
idx = strfind(case_id, '_cut');
if ~isempty(idx)
    patient_id = case_id(1:idx(1)-1);
    return;
end
idx = strfind(case_id, '_process');
if ~isempty(idx)
    patient_id = case_id(1:idx(1)-1);
    return;
end
patient_id = case_id;
end

function stage = assign_stage_labels(nWindow, pre_idx, ictal_idx, post_idx)
stage = repmat({'unassigned'}, nWindow, 1);
pre_idx = sanitize_idx(pre_idx, nWindow);
ictal_idx = sanitize_idx(ictal_idx, nWindow);
post_idx = sanitize_idx(post_idx, nWindow);
ictal_idx = setdiff(ictal_idx, pre_idx, 'stable');
post_idx = setdiff(post_idx, union(pre_idx, ictal_idx), 'stable');
splitPts = round(linspace(0, numel(ictal_idx), 4));
early_idx = ictal_idx(1:splitPts(2));
mid_idx = ictal_idx(splitPts(2)+1:splitPts(3));
late_idx = ictal_idx(splitPts(3)+1:end);
stage(pre_idx) = {'pre-ictal'};
stage(early_idx) = {'early'};
stage(mid_idx) = {'mid'};
stage(late_idx) = {'late'};
stage(post_idx) = {'post-ictal'};
end

function idx = sanitize_idx(idx, nWindow)
idx = idx(:)';
idx = idx(idx >= 1 & idx <= nWindow);
idx = unique(idx, 'stable');
end

function hyperEdges = extract_hyperedges_at_threshold(datanew, threshold)
[nChannel, ~] = size(datanew);
phase = angle(hilbert(datanew.')).';
z = exp(1i * phase);
plvMatrix = abs((z * z') ./ size(z, 2));
plvMatrix(1:nChannel+1:end) = 0;

A = plvMatrix >= threshold;
A(1:nChannel+1:end) = false;
A = A | A.';

[row, col] = find(triu(A, 1));
hyperEdges2 = arrayfun(@(i) [row(i), col(i)], 1:numel(row), 'UniformOutput', false)';

tri = [];
for i = 1:nChannel-2
    neigh = find(A(i, i+1:end)) + i;
    for jj = 1:numel(neigh)-1
        j = neigh(jj);
        ks = neigh(jj+1:end);
        ks = ks(A(j, ks));
        if ~isempty(ks)
            tri = [tri; [repmat(i, numel(ks), 1), repmat(j, numel(ks), 1), ks(:)]];
        end
    end
end
if isempty(tri)
    hyperEdges3 = cell(0, 1);
else
    tri = unique(sort(tri, 2), 'rows', 'stable');
    hyperEdges3 = mat2cell(tri, ones(size(tri, 1), 1), 3);
end

hyperEdges = [hyperEdges2(:); hyperEdges3(:)];
end

function W = build_dmwa_from_hyperedges(HE)
m = numel(HE);
W = zeros(m, m);
if m < 2
    return;
end
sizes = cellfun(@numel, HE);
maxNode = max(cellfun(@max, HE));
rowIdx = zeros(sum(sizes), 1);
colIdx = zeros(sum(sizes), 1);
cursor = 0;
for i = 1:m
    nodes = HE{i}(:);
    k = numel(nodes);
    rowIdx(cursor + (1:k)) = i;
    colIdx(cursor + (1:k)) = nodes;
    cursor = cursor + k;
end
H = sparse(rowIdx, colIdx, 1, m, maxNode);
O = full(H * H.');
O(1:m+1:end) = 0;
D = sum(O > 0, 2);
Sizes_i = repmat(sizes(:), 1, m);
Sizes_j = repmat(sizes(:)', m, 1);
D_i = repmat(D, 1, m);
D_j = repmat(D', m, 1);
Coverage = O ./ Sizes_i;
DegreeSum = D_i + D_j;
DegreeBias = D_j ./ DegreeSum;
DegreeBias(DegreeSum == 0) = 0.5;
W = Coverage;
sameOrder = (Sizes_i == Sizes_j) & (O > 0);
W(sameOrder) = W(sameOrder) .* DegreeBias(sameOrder);
W(O == 0) = 0;
W(1:m+1:end) = 0;
end

function pairRep = w_pair_representation(W)
n = size(W, 1);
if n < 2
    pairRep = struct('n', n, 'i', [], 'j', [], 'wU', [], 'wL', [], 'vals', []);
    return;
end
A = (W > 0) | (W.' > 0);
[i, j] = find(triu(A, 1));
idxU = sub2ind([n n], i, j);
idxL = sub2ind([n n], j, i);
wU = W(idxU);
wL = W(idxL);
pairRep = struct('n', n, 'i', i, 'j', j, 'wU', wU, 'wL', wL, 'vals', [wU; wL]);
end

function metrics = weight_shuffle_metrics(pairRep)
nPairs = numel(pairRep.wU);
if nPairs == 0
    metrics = zeros(4, 1);
    return;
end
vals = pairRep.vals(randperm(numel(pairRep.vals)));
wU = vals(1:nPairs);
wL = vals(nPairs+1:end);
metrics = pair_directed_metrics(pairRep.n, pairRep.i, pairRep.j, wU, wL);
end

function metrics = orientation_pair_swap_metrics(pairRep)
nPairs = numel(pairRep.wU);
if nPairs == 0
    metrics = zeros(4, 1);
    return;
end
swap = rand(nPairs, 1) < 0.5;
wU = pairRep.wU;
wL = pairRep.wL;
tmp = wU(swap);
wU(swap) = wL(swap);
wL(swap) = tmp;
metrics = pair_directed_metrics(pairRep.n, pairRep.i, pairRep.j, wU, wL);
end

function metrics = pair_directed_metrics(n, i, j, wU, wL)
if n < 2 || isempty(i)
    metrics = zeros(4, 1);
    return;
end
wU(~isfinite(wU)) = 0;
wL(~isfinite(wL)) = 0;
wU(wU < 0) = 0;
wL(wL < 0) = 0;
denom = sum(max(wU, wL));
if denom > 0
    Sval = sum(min(wU, wL)) / denom;
else
    Sval = 0;
end
sout = accumarray([i; j], [wU; wL], [n 1], @sum, 0);
sin = accumarray([j; i], [wU; wL], [n 1], @sum, 0);
q = (sout - sin) ./ (sout + sin + eps);
sourceLoad = max(q, 0);
sinkLoad = max(-q, 0);
giniGap = gini_nonnegative(sourceLoad) - gini_nonnegative(sinkLoad);
entropyGap = entropy_norm(sinkLoad) - entropy_norm(sourceLoad);
prGap = participation_norm(sinkLoad) - participation_norm(sourceLoad);
metrics = [Sval; giniGap; entropyGap; prGap];
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

