%#ok<*AGROW,*SAGROW>
%% Phase 9: elbow-threshold DMWA and comparator metrics for Paper 1
%
% This script follows the Paper 2 threshold-selection rule:
%   1. One PLV threshold per seizure.
%   2. The threshold is selected from the threshold-density curve elbow.
%   3. The same seizure-specific threshold is applied to every window.
%
% The selected thresholds are read from Paper 2 Supplementary Fig. 1
% source tables copied into the Paper 1 data_tables directory.

clear; clc;

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot)
    projectRoot = pwd;
end
dataDir = fullfile(projectRoot, 'data_tables');
outDir = fullfile(projectRoot, 'stats');
noteDir = fullfile(projectRoot, 'notes');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(noteDir, 'dir'), mkdir(noteDir); end

caseMetaPath = fullfile(dataDir, 'paper1_case_metadata.csv');
thresholdSummaryPath = fullfile(dataDir, 'paper2_threshold_elbow_case_summary.csv');
windowDensityPath = fullfile(dataDir, 'paper2_threshold_elbow_window_density.csv');
outCsv = fullfile(outDir, 'phase9_elbow_threshold_window_metrics.csv');
outNote = fullfile(noteDir, 'phase9_elbow_threshold_dmwa_comparator_extraction.md');

if ~exist(caseMetaPath, 'file')
    error('Missing case metadata: %s', caseMetaPath);
end
if ~exist(thresholdSummaryPath, 'file')
    error('Missing threshold summary: %s', thresholdSummaryPath);
end
if ~exist(windowDensityPath, 'file')
    error('Missing window-density table: %s', windowDensityPath);
end

thresholdSummary = readtable(thresholdSummaryPath, 'TextType', 'string');
windowDensity = readtable(windowDensityPath, 'TextType', 'string');

caseInfo = load_case_metadata(caseMetaPath);
fs = 1024;
winLenSec = 3;

rows = {};
rowIdx = 0;

fprintf('Phase 9 elbow-threshold extraction for %d seizures.\n', numel(caseInfo));

for c = 1:numel(caseInfo)
    case_id = string(caseInfo(c).case_id);
    file_path = caseInfo(c).file_path;
    match = thresholdSummary(thresholdSummary.case_id == case_id, :);
    if height(match) ~= 1
        error('Cannot find one selected threshold for %s', case_id);
    end
    selectedThreshold = match.selected_threshold(1);

    fprintf('[%02d/%02d] %s, threshold %.4f\n', c, numel(caseInfo), case_id, selectedThreshold);

    S0 = load(file_path, 'X1');
    X1 = S0.X1;
    nChannel = size(X1, 1);
    nSample = size(X1, 2);
    nWindow = floor(nSample / fs) - winLenSec + 1;

    stageLabels = assign_stage_labels(nWindow, caseInfo(c).pre_idx, caseInfo(c).ictal_idx, caseInfo(c).post_idx);
    densityRows = windowDensity(windowDensity.case_id == case_id, :);

    for w = 1:nWindow
        phase5 = stageLabels{w};
        if strcmp(phase5, 'unassigned')
            continue;
        end

        startPoint = fs * (w - 1) + 1;
        endPoint = fs * (w - 1 + winLenSec);
        datanew = X1(:, startPoint:endPoint);

        [hyperEdges, pairDensity, meanPLV] = extract_hyperedges_at_threshold(datanew, selectedThreshold);
        W = build_dmwa_from_hyperedges(hyperEdges);

        dMetrics = directed_dmwa_metrics(W);
        uMetrics = undirected_line_metrics(W);
        oMetrics = order_metrics(hyperEdges);

        densityMatch = densityRows(densityRows.window_index == w, :);
        if height(densityMatch) == 1
            pairDensityFromPaper2 = densityMatch.pair_density_at_elbow(1);
            meanPLVFromPaper2 = densityMatch.mean_plv(1);
        else
            pairDensityFromPaper2 = NaN;
            meanPLVFromPaper2 = NaN;
        end

        rowIdx = rowIdx + 1;
        rows(rowIdx, :) = { ...
            char(case_id), get_case_patient_id(caseInfo(c), char(case_id)), w, phase5, selectedThreshold, ...
            nChannel, pairDensity, meanPLV, pairDensityFromPaper2, meanPLVFromPaper2, ...
            oMetrics.n_hyperedges, oMetrics.mean_order, oMetrics.frac_order2, oMetrics.frac_order3, oMetrics.order_entropy_norm, ...
            dMetrics.S, dMetrics.source_sink_gini_gap, dMetrics.sink_minus_source_entropy_gap, dMetrics.sink_minus_source_pr_gap, ...
            uMetrics.n_edges_undirected, uMetrics.edge_density_undirected, uMetrics.mean_weight_nonzero, ...
            uMetrics.weight_sum_norm, uMetrics.degree_gini, uMetrics.degree_entropy_norm, ...
            uMetrics.degree_pr_norm, uMetrics.two_path_triple_density ...
            };
    end
end

varNames = {'case_id','patient_id','window_idx','phase5','selected_threshold', ...
    'n_channels','pair_density_elbow','mean_plv','pair_density_from_paper2','mean_plv_from_paper2', ...
    'n_hyperedges','mean_order','frac_order2','frac_order3','order_entropy_norm', ...
    'S','source_sink_gini_gap','sink_minus_source_entropy_gap','sink_minus_source_pr_gap', ...
    'n_edges_undirected','edge_density_undirected','mean_weight_nonzero', ...
    'weight_sum_norm','degree_gini','degree_entropy_norm','degree_pr_norm','two_path_triple_density'};

T = cell2table(rows, 'VariableNames', varNames);
writetable(T, outCsv);

fid = fopen(outNote, 'w');
fprintf(fid, '# Phase 9 elbow-threshold DMWA comparator extraction\n\n');
fprintf(fid, '- Threshold source: `%s`\n', thresholdSummaryPath);
fprintf(fid, '- Output CSV: `%s`\n', outCsv);
fprintf(fid, '- Rows written: %d\n', height(T));
fprintf(fid, '- Rule: one seizure-specific PLV threshold selected at the threshold-density elbow, then applied to all windows from that seizure.\n');
fprintf(fid, '- This supersedes older per-window 55%% quantile thresholding for new Paper 1 comparator analyses.\n');
fclose(fid);

fprintf('Wrote %d rows to %s\n', height(T), outCsv);

function [hyperEdges, pairDensity, meanPLV] = extract_hyperedges_at_threshold(datanew, threshold)
    [nChannel, ~] = size(datanew);
    phase = angle(hilbert(datanew.')).';
    z = exp(1i * phase);
    plvMatrix = abs((z * z') ./ size(z, 2));
    plvMatrix(1:nChannel+1:end) = 0;

    upperMask = triu(true(nChannel), 1);
    edgePLV = plvMatrix(upperMask);
    meanPLV = mean(edgePLV, 'omitnan');

    A = plvMatrix >= threshold;
    A(1:nChannel+1:end) = false;
    A = A | A.';
    pairDensity = nnz(triu(A, 1)) / nnz(upperMask);

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

function metrics = directed_dmwa_metrics(W)
    n = size(W, 1);
    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    if n < 2
        metrics.S = 0;
        metrics.source_sink_gini_gap = 0;
        metrics.sink_minus_source_entropy_gap = 0;
        metrics.sink_minus_source_pr_gap = 0;
        return;
    end
    W(1:n+1:end) = 0;
    mask = triu(true(n), 1);
    Wij = W(mask);
    Wji = W.'; Wji = Wji(mask);
    denom = sum(max(Wij, Wji));
    if denom > 0
        metrics.S = sum(min(Wij, Wji)) / denom;
    else
        metrics.S = 0;
    end
    sout = sum(W, 2);
    sin = sum(W, 1)';
    q = (sout - sin) ./ (sout + sin + eps);
    sourceLoad = max(q, 0);
    sinkLoad = max(-q, 0);
    metrics.source_sink_gini_gap = gini_nonnegative(sourceLoad) - gini_nonnegative(sinkLoad);
    metrics.sink_minus_source_entropy_gap = entropy_norm(sinkLoad) - entropy_norm(sourceLoad);
    metrics.sink_minus_source_pr_gap = participation_norm(sinkLoad) - participation_norm(sourceLoad);
end

function metrics = undirected_line_metrics(W)
    n = size(W, 1);
    if n < 2
        metrics.n_edges_undirected = 0;
        metrics.edge_density_undirected = 0;
        metrics.mean_weight_nonzero = 0;
        metrics.weight_sum_norm = 0;
        metrics.degree_gini = 0;
        metrics.degree_entropy_norm = 0;
        metrics.degree_pr_norm = 0;
        metrics.two_path_triple_density = 0;
        return;
    end
    U = max(W, W.');
    U(1:n+1:end) = 0;
    A = U > 0;
    upperMask = triu(true(n), 1);
    weights = U(upperMask);
    weights = weights(weights > 0);
    metrics.n_edges_undirected = numel(weights);
    metrics.edge_density_undirected = metrics.n_edges_undirected / nchoosek(n, 2);
    if isempty(weights)
        metrics.mean_weight_nonzero = 0;
    else
        metrics.mean_weight_nonzero = mean(weights);
    end
    metrics.weight_sum_norm = sum(weights) / max(1, n);
    deg = sum(U, 2);
    metrics.degree_gini = gini_nonnegative(deg);
    metrics.degree_entropy_norm = entropy_norm(deg);
    metrics.degree_pr_norm = participation_norm(deg);
    binaryDeg = sum(A, 2);
    twoPathTriples = sum(binaryDeg .* max(binaryDeg - 1, 0) / 2);
    if n >= 3
        metrics.two_path_triple_density = twoPathTriples / nchoosek(n, 3);
    else
        metrics.two_path_triple_density = 0;
    end
end

function metrics = order_metrics(HE)
    n = numel(HE);
    metrics.n_hyperedges = n;
    if n == 0
        metrics.mean_order = NaN;
        metrics.frac_order2 = NaN;
        metrics.frac_order3 = NaN;
        metrics.order_entropy_norm = NaN;
        return;
    end
    orders = cellfun(@numel, HE);
    metrics.mean_order = mean(orders);
    metrics.frac_order2 = mean(orders == 2);
    metrics.frac_order3 = mean(orders == 3);
    metrics.order_entropy_norm = entropy_norm([sum(orders == 2); sum(orders == 3)]);
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

function caseInfo = load_case_metadata(caseMetaPath)
    meta = readtable(caseMetaPath, 'TextType', 'string');
    required = ["case_id","file_path","pre_start","pre_end", ...
        "ictal_start","ictal_end","post_start","post_end"];
    missing = setdiff(required, string(meta.Properties.VariableNames));
    if ~isempty(missing)
        error('Case metadata is missing required columns: %s', strjoin(missing, ', '));
    end

    projectRoot = fileparts(fileparts(caseMetaPath));
    hasPatient = ismember("patient_id", string(meta.Properties.VariableNames));
    caseInfo = struct([]);
    for i = 1:height(meta)
        caseInfo(i).case_id = char(meta.case_id(i));
        caseInfo(i).file_path = resolve_local_path(char(meta.file_path(i)), projectRoot);
        if hasPatient && strlength(meta.patient_id(i)) > 0
            caseInfo(i).patient_id = char(meta.patient_id(i));
        else
            caseInfo(i).patient_id = '';
        end
        caseInfo(i).pre_idx = meta.pre_start(i):meta.pre_end(i);
        caseInfo(i).ictal_idx = meta.ictal_start(i):meta.ictal_end(i);
        caseInfo(i).post_idx = meta.post_start(i):meta.post_end(i);
    end
end

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

function patient_id = get_case_patient_id(info, case_id)
    if isfield(info, 'patient_id') && ~isempty(info.patient_id)
        patient_id = info.patient_id;
    else
        patient_id = infer_patient_id(case_id);
    end
end

