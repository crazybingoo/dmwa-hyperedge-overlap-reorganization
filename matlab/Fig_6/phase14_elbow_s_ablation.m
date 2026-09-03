%#ok<*AGROW,*SAGROW>
%% Phase 14: elbow-threshold S ablation for Paper 1 Fig. 6
% Recomputes Original, Coverage-only and DegreeBias-only S variants under
% the Paper 2 threshold-density elbow rule. This replaces legacy Fig. 6
% ablation panels that used the older 55% PLV quantile hyperedge rule.

clear; clc;

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot)
    projectRoot = pwd;
end
dataDir = fullfile(projectRoot, 'data_tables');
statsDir = fullfile(projectRoot, 'stats');
noteDir = fullfile(projectRoot, 'notes');
if ~exist(statsDir, 'dir'), mkdir(statsDir); end
if ~exist(noteDir, 'dir'), mkdir(noteDir); end

caseMetaPath = fullfile(dataDir, 'paper1_case_metadata.csv');
thresholdPath = fullfile(dataDir, 'paper2_threshold_elbow_case_summary.csv');
caseMeta = readtable(caseMetaPath, 'TextType', 'string');
thresholdSummary = readtable(thresholdPath, 'TextType', 'string');

caseStart = read_env_int('PHASE14_CASE_START', 1);
caseLimit = read_env_int('PHASE14_CASE_LIMIT', height(caseMeta));
caseEnd = read_env_int('PHASE14_CASE_END', height(caseMeta));
caseStart = max(1, caseStart);
caseEnd = min(height(caseMeta), min(caseEnd, caseStart + caseLimit - 1));

outputSuffix = getenv('PHASE14_OUTPUT_SUFFIX');
if strlength(string(outputSuffix)) > 0
    outputSuffix = ['_' char(outputSuffix)];
else
    outputSuffix = '';
end

winLenSec = 3;
windowOut = fullfile(statsDir, ['phase14_elbow_s_ablation_window_summary' outputSuffix '.csv']);
stageOut = fullfile(statsDir, ['phase14_elbow_s_ablation_stage_summary' outputSuffix '.csv']);
noteOut = fullfile(noteDir, ['phase14_elbow_s_ablation' outputSuffix '.md']);

windowRows = {};
stageTables = cell(caseEnd - caseStart + 1, 1);
rowIdx = 0;

fprintf('Phase 14 elbow S ablation: cases %d-%d (%d cases).\n', ...
    caseStart, caseEnd, caseEnd - caseStart + 1);

for c = caseStart:caseEnd
    caseId = string(caseMeta.case_id(c));
    patientId = get_patient_id(caseId);
    filePath = resolve_local_path(char(caseMeta.file_path(c)), projectRoot);

    match = thresholdSummary(thresholdSummary.case_id == caseId, :);
    if height(match) ~= 1
        error('Cannot find exactly one elbow threshold for %s', caseId);
    end
    threshold = match.selected_threshold(1);
    fprintf('[%02d/%02d] %s, threshold %.4f\n', c, height(caseMeta), caseId, threshold);

    matVars = whos('-file', filePath);
    hasFs = any(strcmp({matVars.name}, 'fs'));
    if hasFs
        S = load(filePath, 'X1', 'fs');
    else
        S = load(filePath, 'X1');
    end
    if ~isfield(S, 'X1')
        error('Missing X1 in %s', filePath);
    end
    X1 = double(S.X1);
    if isfield(S, 'fs') && ~isempty(S.fs)
        fs = double(S.fs);
    else
        fs = 1024;
    end

    nSample = size(X1, 2);
    nWindow = floor(nSample / fs) - winLenSec + 1;
    if nWindow < 1
        warning('%s has insufficient samples for 3 s windows.', caseId);
        continue;
    end

    stageLabels = assign_stage_labels(nWindow, ...
        caseMeta.pre_start(c):caseMeta.pre_end(c), ...
        caseMeta.ictal_start(c):caseMeta.ictal_end(c), ...
        caseMeta.post_start(c):caseMeta.post_end(c));

    caseRows = {};
    caseRowIdx = 0;
    for w = 1:nWindow
        phase5 = string(stageLabels{w});
        if phase5 == "unassigned"
            continue;
        end

        startPoint = fs * (w - 1) + 1;
        endPoint = fs * (w - 1 + winLenSec);
        datanew = X1(:, startPoint:endPoint);

        hyperEdges = extract_hyperedges_at_threshold(datanew, threshold);
        [Worig, Wcov, Wdeg] = build_ablation_matrices(hyperEdges);

        row = {char(caseId), char(patientId), w, char(phase5), threshold, numel(hyperEdges), ...
            calculate_S_from_W(Worig), calculate_S_from_W(Wcov), calculate_S_from_W(Wdeg)};
        rowIdx = rowIdx + 1;
        windowRows(rowIdx, :) = row;
        caseRowIdx = caseRowIdx + 1;
        caseRows(caseRowIdx, :) = row;
    end

    if ~isempty(caseRows)
        caseTable = rows_to_table(caseRows);
        stageTables{c - caseStart + 1} = aggregate_stage_metrics(caseTable);
    end
end

windowTable = rows_to_table(windowRows);
stageTables = stageTables(~cellfun(@isempty, stageTables));
if isempty(stageTables)
    stageTable = table();
else
    stageTable = vertcat(stageTables{:});
end

writetable(windowTable, windowOut);
writetable(stageTable, stageOut);

fid = fopen(noteOut, 'w');
fprintf(fid, '# Phase 14 elbow-threshold S ablation\n\n');
fprintf(fid, '- Case rows processed: %d-%d\n', caseStart, caseEnd);
fprintf(fid, '- Window length: %d s\n', winLenSec);
fprintf(fid, '- Sampling rate: read `fs` from `.mat` when present; otherwise fallback to 1024 Hz to match Phase 9.\n');
fprintf(fid, '- Threshold source: `%s`\n', thresholdPath);
fprintf(fid, '- Rule: one seizure-specific threshold-density elbow applied to all windows from that seizure.\n');
fprintf(fid, '- Original: Coverage plus same-order DegreeBias.\n');
fprintf(fid, '- Coverage-only: `O / size(e_i)` without same-order DegreeBias.\n');
fprintf(fid, '- DegreeBias-only: `O * DegreeBias`, ignoring hyperedge order normalization.\n');
fprintf(fid, '- Window rows: %d\n', height(windowTable));
fprintf(fid, '- Stage rows: %d\n', height(stageTable));
fprintf(fid, '\n## Outputs\n\n');
fprintf(fid, '- `%s`\n', windowOut);
fprintf(fid, '- `%s`\n', stageOut);
fclose(fid);

fprintf('Wrote window table: %s\n', windowOut);
fprintf('Wrote stage table: %s\n', stageOut);
fprintf('Wrote note: %s\n', noteOut);

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

function value = read_env_int(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
else
    value = str2double(raw);
    if isnan(value), value = defaultValue; end
end
end

function patientId = get_patient_id(caseId)
patientId = regexprep(char(caseId), '_(cut\d+|process)$', '');
end

function stage = assign_stage_labels(nWindow, preIdx, ictalIdx, postIdx)
stage = repmat({'unassigned'}, nWindow, 1);
preIdx = sanitize_idx(preIdx, nWindow);
ictalIdx = sanitize_idx(ictalIdx, nWindow);
postIdx = sanitize_idx(postIdx, nWindow);
ictalIdx = setdiff(ictalIdx, preIdx, 'stable');
postIdx = setdiff(postIdx, union(preIdx, ictalIdx), 'stable');

if ~isempty(preIdx)
    stage(preIdx) = {'pre-ictal'};
end
if ~isempty(ictalIdx)
    splitPts = round(linspace(0, numel(ictalIdx), 4));
    earlyIdx = ictalIdx(1:splitPts(2));
    midIdx = ictalIdx(splitPts(2)+1:splitPts(3));
    lateIdx = ictalIdx(splitPts(3)+1:end);
    stage(earlyIdx) = {'early'};
    stage(midIdx) = {'mid'};
    stage(lateIdx) = {'late'};
end
if ~isempty(postIdx)
    stage(postIdx) = {'post-ictal'};
end
end

function idx = sanitize_idx(idx, nWindow)
idx = idx(:)';
idx = idx(idx >= 1 & idx <= nWindow);
idx = unique(idx, 'stable');
end

function hyperEdges = extract_hyperedges_at_threshold(datanew, threshold)
[nChannel, ~] = size(datanew);
phase = angle(hilbert(datanew.')).';
plvMatrix = eye(nChannel);
for i = 1:nChannel-1
    for j = i+1:nChannel
        phaseDiff = phase(i, :) - phase(j, :);
        plv = abs(mean(exp(1i * phaseDiff)));
        plvMatrix(i, j) = plv;
        plvMatrix(j, i) = plv;
    end
end
A = plvMatrix >= threshold;
A(1:nChannel+1:end) = false;

[row, col] = find(triu(A, 1));
if isempty(row)
    hyperEdges2 = cell(0, 1);
else
    hyperEdges2 = arrayfun(@(i) [row(i), col(i)], 1:numel(row), 'UniformOutput', false)';
end

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

function [Worig, Wcov, Wdeg] = build_ablation_matrices(HE)
m = numel(HE);
Worig = zeros(m, m);
Wcov = zeros(m, m);
Wdeg = zeros(m, m);
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
Coverage = O ./ Sizes_i;
DegreeSum = repmat(D, 1, m) + repmat(D', m, 1);
DegreeBias = repmat(0.5, m, m);
degreeBiasRaw = repmat(D', m, 1) ./ DegreeSum;
degreeMask = DegreeSum > 0;
DegreeBias(degreeMask) = degreeBiasRaw(degreeMask);

Worig = Coverage;
sameOrder = (Sizes_i == Sizes_j) & (O > 0);
Worig(sameOrder) = Worig(sameOrder) .* DegreeBias(sameOrder);
Worig(O == 0) = 0;
Worig(1:m+1:end) = 0;

Wcov = Coverage;
Wcov(O == 0) = 0;
Wcov(1:m+1:end) = 0;

Wdeg = O .* DegreeBias;
Wdeg(O == 0) = 0;
Wdeg(1:m+1:end) = 0;
end

function S = calculate_S_from_W(W)
n = size(W, 1);
if n < 2
    S = 0;
    return;
end
W(~isfinite(W)) = 0;
W(W < 0) = 0;
W(1:n+1:end) = 0;
mask = triu(true(n), 1);
Wij = W(mask);
Wji = W.';
Wji = Wji(mask);
denom = sum(max(Wij, Wji));
if denom > 0
    S = sum(min(Wij, Wji)) / denom;
else
    S = 0;
end
end

function T = rows_to_table(rows)
varNames = {'case_id','patient_id','window_idx','phase5','selected_threshold', ...
    'n_hyperedges','S_orig','S_cov','S_deg'};
if isempty(rows)
    T = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
else
    T = cell2table(rows, 'VariableNames', varNames);
end
numericVars = varNames([3, 5:numel(varNames)]);
for v = 1:numel(numericVars)
    name = numericVars{v};
    if iscell(T.(name))
        T.(name) = cell2mat(T.(name));
    end
end
end

function stageT = aggregate_stage_metrics(winT)
stageOrder = {'pre-ictal','early','mid','late','post-ictal'};
stageT = table();
if isempty(winT), return; end
caseId = winT.case_id{1};
patientId = winT.patient_id{1};
for s = 1:numel(stageOrder)
    phase5 = stageOrder{s};
    idx = strcmp(winT.phase5, phase5);
    if ~any(idx), continue; end
    row = table({caseId}, {patientId}, {phase5}, sum(idx), ...
        mean(winT.S_orig(idx), 'omitnan'), mean(winT.S_cov(idx), 'omitnan'), mean(winT.S_deg(idx), 'omitnan'), ...
        std(winT.S_orig(idx), 'omitnan'), std(winT.S_cov(idx), 'omitnan'), std(winT.S_deg(idx), 'omitnan'), ...
        'VariableNames', {'case_id','patient_id','phase5','n_windows', ...
        'S_orig_mean','S_cov_mean','S_deg_mean','S_orig_std','S_cov_std','S_deg_std'});
    stageT = [stageT; row];
end

preIdx = strcmp(stageT.phase5, 'pre-ictal');
if any(preIdx)
    preOrig = stageT.S_orig_mean(preIdx);
    preCov = stageT.S_cov_mean(preIdx);
    preDeg = stageT.S_deg_mean(preIdx);
    stageT.Delta_orig = stageT.S_orig_mean - preOrig;
    stageT.Delta_cov = stageT.S_cov_mean - preCov;
    stageT.Delta_deg = stageT.S_deg_mean - preDeg;
else
    stageT.Delta_orig = nan(height(stageT), 1);
    stageT.Delta_cov = nan(height(stageT), 1);
    stageT.Delta_deg = nan(height(stageT), 1);
end
end

