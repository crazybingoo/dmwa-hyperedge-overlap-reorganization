%#ok<*AGROW,*SAGROW>
%% Phase 13: elbow-threshold 13-motif extraction for Paper 1 Fig. 5
% Recomputes local 3-motif composition under the Paper 2 threshold-density
% elbow rule. This replaces legacy motif source tables that used the fixed
% 55% PLV quantile rule.

clear; clc;

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot)
    projectRoot = pwd;
end
dataDir = fullfile(projectRoot, 'data_tables');
statsDir = fullfile(projectRoot, 'stats');
analysisDir = fullfile(projectRoot, 'analysis_tables');
noteDir = fullfile(projectRoot, 'notes');
if ~exist(statsDir, 'dir'), mkdir(statsDir); end
if ~exist(analysisDir, 'dir'), mkdir(analysisDir); end
if ~exist(noteDir, 'dir'), mkdir(noteDir); end

caseMetaPath = fullfile(dataDir, 'paper1_case_metadata.csv');
thresholdPath = fullfile(dataDir, 'paper2_threshold_elbow_case_summary.csv');
caseMeta = readtable(caseMetaPath, 'TextType', 'string');
thresholdSummary = readtable(thresholdPath, 'TextType', 'string');

caseStart = read_env_int('PHASE13_CASE_START', 1);
caseLimit = read_env_int('PHASE13_CASE_LIMIT', height(caseMeta));
caseEnd = read_env_int('PHASE13_CASE_END', height(caseMeta));
caseStart = max(1, caseStart);
caseEnd = min(height(caseMeta), min(caseEnd, caseStart + caseLimit - 1));

outputSuffix = getenv('PHASE13_OUTPUT_SUFFIX');
if strlength(string(outputSuffix)) > 0
    outputSuffix = ['_' char(outputSuffix)];
else
    outputSuffix = '';
end

winLenSec = 3;
maxExactNodes = read_env_int('PHASE13_MOTIF_MAX_EXACT_N', 180);
sampleTriads = read_env_int('PHASE13_MOTIF_SAMPLE_TRIADS', 50000);
maxWindowsPerStage = read_env_int('PHASE13_MAX_WINDOWS_PER_STAGE', inf);
randomSeed = read_env_int('PHASE13_RANDOM_SEED', 260508);
edgeEps = read_env_double('PHASE13_EDGE_EPS', 1e-12);
reciprocityRatio = read_env_double('PHASE13_RECIPROCITY_RATIO', 0.85);
rng(randomSeed, 'twister');

motifNames = {'021D','021U','021C','111D','111U','030T','030C', ...
    '201','120D','120U','120C','210','300'};
triadLUT = build_triad13_lut(motifNames);

windowOut = fullfile(statsDir, ['phase13_elbow_motif_window_summary' outputSuffix '.csv']);
stageOut = fullfile(statsDir, ['phase13_elbow_motif_stage_summary' outputSuffix '.csv']);
analysisWindowOut = fullfile(analysisDir, ['motif_group_window_summary_elbow' outputSuffix '.csv']);
analysisStageOut = fullfile(analysisDir, ['motif_group_stage_summary_elbow' outputSuffix '.csv']);
noteOut = fullfile(noteDir, ['phase13_elbow_motif13_extraction' outputSuffix '.md']);

windowRows = {};
rowIdx = 0;
stageTables = cell(caseEnd - caseStart + 1, 1);

fprintf('Phase 13 elbow motif extraction: cases %d-%d (%d cases).\n', ...
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

    preIdx = caseMeta.pre_start(c):caseMeta.pre_end(c);
    ictalIdx = caseMeta.ictal_start(c):caseMeta.ictal_end(c);
    postIdx = caseMeta.post_start(c):caseMeta.post_end(c);
    stageLabels = assign_stage_labels(nWindow, preIdx, ictalIdx, postIdx);

    caseRows = table();
    stageQuota = containers.Map({'pre-ictal','early','mid','late','post-ictal'}, zeros(1, 5));
    for w = 1:nWindow
        stage = string(stageLabels{w});
        if stage == "unassigned"
            continue;
        end
        stageKey = char(stage);
        if isfinite(maxWindowsPerStage)
            if stageQuota(stageKey) >= maxWindowsPerStage
                continue;
            end
            stageQuota(stageKey) = stageQuota(stageKey) + 1;
        end

        startPoint = fs * (w - 1) + 1;
        endPoint = fs * (w - 1 + winLenSec);
        datanew = X1(:, startPoint:endPoint);

        hyperEdges = extract_hyperedges_at_threshold(datanew, threshold);
        W = build_dmwa_from_hyperedges(hyperEdges);
        nnzWPositive = nnz(W > edgeEps);
        A = build_adjacency_from_W(W, edgeEps, reciprocityRatio);
        nnzADirected = nnz(A);

        [counts, numTriplets, numConnected, censusMode, triadsEvaluated, samplingFraction] = ...
            triad13_counts(A, triadLUT, maxExactNodes, sampleTriads);
        metrics = motif_group_metrics(counts, numTriplets, numConnected, size(W, 1), nnzWPositive, nnzADirected);

        rowIdx = rowIdx + 1;
        windowRows(rowIdx, :) = [{char(caseId), char(patientId), w, char(stage), ...
            censusMode, triadsEvaluated, samplingFraction}, ...
            num2cell(metrics)];
    end

    if rowIdx > 0
        allWindowPartial = rows_to_table(windowRows);
        caseRows = allWindowPartial(strcmp(allWindowPartial.case_id, char(caseId)), :);
        stageTables{c - caseStart + 1} = aggregate_stage_metrics(caseRows);
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
writetable(windowTable, analysisWindowOut);
writetable(stageTable, analysisStageOut);

fid = fopen(noteOut, 'w');
fprintf(fid, '# Phase 13 elbow-threshold motif13 extraction\n\n');
fprintf(fid, '- Case rows processed: %d-%d\n', caseStart, caseEnd);
fprintf(fid, '- Window length: %d s\n', winLenSec);
fprintf(fid, '- Sampling rate: read `fs` from `.mat` when present; otherwise fallback to 1024 Hz to match Phase 9 and legacy Fig. 5 scripts.\n');
fprintf(fid, '- Threshold source: `%s`\n', thresholdPath);
fprintf(fid, '- Rule: one seizure-specific threshold-density elbow applied to all windows from that seizure.\n');
fprintf(fid, '- Directional motif adjacency: stronger dyadic DMWA direction unless min/max >= %.3f, then reciprocal.\n', reciprocityRatio);
fprintf(fid, '- Edge epsilon for motif adjacency: %.3g.\n', edgeEps);
fprintf(fid, '- Motif census: exact for DMWA node count <= %d; otherwise %d uniformly sampled triads with fixed seed %d.\n', ...
    maxExactNodes, sampleTriads, randomSeed);
if isfinite(maxWindowsPerStage)
    fprintf(fid, '- Window quota: at most %d windows per stage and case.\n', maxWindowsPerStage);
else
    fprintf(fid, '- Window quota: all assigned windows.\n');
end
fprintf(fid, '- Motif taxonomy: 13 connected directed triads grouped as simple, reciprocal, dense and complex.\n');
fprintf(fid, '- Window rows: %d\n', height(windowTable));
fprintf(fid, '- Stage rows: %d\n', height(stageTable));
fprintf(fid, '\n## Outputs\n\n');
fprintf(fid, '- `%s`\n', windowOut);
fprintf(fid, '- `%s`\n', stageOut);
fprintf(fid, '- `%s`\n', analysisWindowOut);
fprintf(fid, '- `%s`\n', analysisStageOut);
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

function value = read_env_double(name, defaultValue)
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
stage(preIdx) = {'pre-ictal'};

if ~isempty(ictalIdx)
    nIctal = numel(ictalIdx);
    cut1 = round(nIctal / 3);
    cut2 = round(2 * nIctal / 3);
    earlyIdx = ictalIdx(1:max(cut1, 1));
    midIdx = ictalIdx(max(cut1, 1)+1:max(cut2, max(cut1, 1)));
    lateIdx = ictalIdx(max(cut2, max(cut1, 1))+1:end);
    stage(earlyIdx) = {'early'};
    stage(midIdx) = {'mid'};
    stage(lateIdx) = {'late'};
end

stage(postIdx) = {'post-ictal'};
end

function idx = sanitize_idx(idx, nWindow)
idx = idx(:)';
idx = idx(idx >= 1 & idx <= nWindow);
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
hyperEdges2 = arrayfun(@(i) [row(i), col(i)], 1:numel(row), 'UniformOutput', false)';

tri = [];
for i = 1:nChannel-2
    for j = i+1:nChannel-1
        if ~A(i, j), continue; end
        ks = find(A(i, j+1:end) & A(j, j+1:end)) + j;
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
O = H * H.';
O = full(O);
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
W = Coverage;
sameOrder = (Sizes_i == Sizes_j) & (O > 0);
W(sameOrder) = W(sameOrder) .* DegreeBias(sameOrder);
W(O == 0) = 0;
W(1:m+1:end) = 0;
end

function A = build_adjacency_from_W(W, edgeEps, reciprocityRatio)
n = size(W, 1);
A = false(n, n);
for i = 1:n-1
    for j = i+1:n
        wij = W(i, j);
        wji = W(j, i);
        vmax = max(wij, wji);
        if vmax <= edgeEps
            continue;
        end
        vmin = min(wij, wji);
        ratio = vmin / (vmax + eps);
        if ratio >= reciprocityRatio
            A(i, j) = true;
            A(j, i) = true;
        elseif wij > wji
            A(i, j) = true;
        elseif wji > wij
            A(j, i) = true;
        else
            A(i, j) = true;
            A(j, i) = true;
        end
    end
end
A(1:n+1:end) = false;
end

function lut = build_triad13_lut(motifNames)
templates = get_triad_templates(motifNames);
permList = perms(1:3);
lut = zeros(64, 1);
for m = 1:numel(motifNames)
    T = templates{m};
    for p = 1:size(permList, 1)
        ord = permList(p, :);
        Tp = T(ord, ord);
        code = triad_code(Tp);
        lut(code + 1) = m;
    end
end
end

function templates = get_triad_templates(motifNames)
templateMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
templateMap('021D') = logical([0 1 1; 0 0 0; 0 0 0]);
templateMap('021U') = logical([0 0 0; 1 0 0; 1 0 0]);
templateMap('021C') = logical([0 1 0; 0 0 1; 0 0 0]);
templateMap('111D') = logical([0 1 0; 1 0 0; 1 0 0]);
templateMap('111U') = logical([0 1 1; 1 0 0; 0 0 0]);
templateMap('030T') = logical([0 1 1; 0 0 1; 0 0 0]);
templateMap('030C') = logical([0 1 0; 0 0 1; 1 0 0]);
templateMap('201') = logical([0 1 1; 1 0 0; 1 0 0]);
templateMap('120D') = logical([0 1 1; 1 0 1; 0 0 0]);
templateMap('120U') = logical([0 1 0; 1 0 0; 1 1 0]);
templateMap('120C') = logical([0 1 1; 1 0 0; 0 1 0]);
templateMap('210') = logical([0 1 1; 1 0 1; 1 0 0]);
templateMap('300') = logical([0 1 1; 1 0 1; 1 1 0]);
templates = cell(size(motifNames));
for m = 1:numel(motifNames)
    templates{m} = templateMap(motifNames{m});
end
end

function code = triad_code(A3)
bits = [A3(1,2), A3(1,3), A3(2,1), A3(2,3), A3(3,1), A3(3,2)];
code = bits * [1; 2; 4; 8; 16; 32];
end

function [counts, numTriplets, numConnected, mode, triadsEvaluated, samplingFraction] = ...
    triad13_counts(B, triadLUT, maxExactNodes, sampleTriads)
n = size(B, 1);
counts = zeros(1, 13);
mode = 'exact';
triadsEvaluated = 0;
samplingFraction = NaN;
if n < 3
    numTriplets = 0;
    numConnected = 0;
    return;
end
numTriplets = nchoosek(n, 3);
numConnected = 0;

if n <= maxExactNodes || sampleTriads >= numTriplets
    mode = 'exact';
    triadsEvaluated = numTriplets;
    samplingFraction = 1;
    for i = 1:n-2
        for j = i+1:n-1
            k = (j+1:n)';
            ij = repmat(B(i, j), numel(k), 1);
            ji = repmat(B(j, i), numel(k), 1);
            codes = ...
                ij + ...
                2 * B(i, k).' + ...
                4 * ji + ...
                8 * B(j, k).' + ...
                16 * B(k, i) + ...
                32 * B(k, j);
            motifIdx = triadLUT(codes + 1);
            keep = motifIdx > 0;
            numConnected = numConnected + nnz(keep);
            if any(keep)
                counts = counts + accumarray(motifIdx(keep), 1, [13 1]).';
            end
        end
    end
else
    mode = 'sampled';
    triadsEvaluated = min(sampleTriads, numTriplets);
    samplingFraction = triadsEvaluated / numTriplets;
    rawCounts = zeros(1, 13);
    filled = 0;
    while filled < triadsEvaluated
        batchSize = max(ceil((triadsEvaluated - filled) * 1.4), 1000);
        cand = randi(n, batchSize, 3);
        keepDistinct = cand(:, 1) ~= cand(:, 2) & cand(:, 1) ~= cand(:, 3) & cand(:, 2) ~= cand(:, 3);
        cand = sort(cand(keepDistinct, :), 2);
        if isempty(cand)
            continue;
        end
        take = min(size(cand, 1), triadsEvaluated - filled);
        cand = cand(1:take, :);
        i = cand(:, 1); j = cand(:, 2); k = cand(:, 3);
        codes = ...
            B(sub2ind([n n], i, j)) + ...
            2 * B(sub2ind([n n], i, k)) + ...
            4 * B(sub2ind([n n], j, i)) + ...
            8 * B(sub2ind([n n], j, k)) + ...
            16 * B(sub2ind([n n], k, i)) + ...
            32 * B(sub2ind([n n], k, j));
        motifIdx = triadLUT(codes + 1);
        keep = motifIdx > 0;
        if any(keep)
            rawCounts = rawCounts + accumarray(motifIdx(keep), 1, [13 1]).';
        end
        filled = filled + take;
    end
    counts = rawCounts ./ triadsEvaluated .* numTriplets;
    numConnected = sum(counts, 'omitnan');
end
end

function metrics = motif_group_metrics(counts, numTriplets, numConnected, numHyperedges, nnzWPositive, nnzADirected)
allCount = sum(counts, 'omitnan');
simpleIdx = [1 2 3];
fflIdx = 6;
cycleIdx = 7;
reciprocalIdx = [4 5 8 9 10 11 12 13];
denseIdx = [9 10 11 12 13];
complexIdx = setdiff(1:13, simpleIdx);

simpleCount = sum(counts(simpleIdx), 'omitnan');
fflCount = sum(counts(fflIdx), 'omitnan');
cycleCount = sum(counts(cycleIdx), 'omitnan');
reciprocalCount = sum(counts(reciprocalIdx), 'omitnan');
denseCount = sum(counts(denseIdx), 'omitnan');
complexCount = sum(counts(complexIdx), 'omitnan');

metrics = [numHyperedges, nnzWPositive, nnzADirected, numTriplets, numConnected, ...
    allCount, simpleCount, fflCount, cycleCount, reciprocalCount, denseCount, complexCount, ...
    safe_div(simpleCount, allCount), safe_div(fflCount, allCount), safe_div(cycleCount, allCount), ...
    safe_div(reciprocalCount, allCount), safe_div(denseCount, allCount), safe_div(complexCount, allCount), ...
    safe_div(reciprocalCount, simpleCount), safe_div(denseCount, simpleCount), ...
    safe_div(fflCount, simpleCount), safe_div(complexCount, simpleCount)];
end

function x = safe_div(a, b)
if b == 0 || isnan(b)
    x = NaN;
else
    x = a ./ b;
end
end

function T = rows_to_table(rows)
varNames = {'case_id','patient_id','window_idx','phase5', ...
    'motif_census_mode','motif_triads_evaluated','motif_sampling_fraction', ...
    'num_hyperedges','nnz_W_positive','nnz_A_directed','num_triplets','num_connected_triplets', ...
    'all_motif_count','simple_count','ffl_count','cycle_count','reciprocal_count','dense_count','complex_count', ...
    'simple_fraction','ffl_fraction','cycle_fraction','reciprocal_fraction','dense_fraction','complex_fraction', ...
    'reciprocal_to_simple','dense_to_simple','ffl_to_simple','complex_to_simple'};
if isempty(rows)
    T = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
else
    T = cell2table(rows, 'VariableNames', varNames);
end
numericVars = varNames([3, 6:numel(varNames)]);
for v = 1:numel(numericVars)
    name = numericVars{v};
    if iscell(T.(name))
        T.(name) = cell2mat(T.(name));
    end
end
end

function stageT = aggregate_stage_metrics(winT)
metricNames = {'motif_triads_evaluated','motif_sampling_fraction', ...
    'num_hyperedges','nnz_W_positive','nnz_A_directed','num_triplets','num_connected_triplets', ...
    'all_motif_count','simple_count','ffl_count','cycle_count', ...
    'reciprocal_count','dense_count','complex_count', ...
    'simple_fraction','ffl_fraction','cycle_fraction','reciprocal_fraction', ...
    'dense_fraction','complex_fraction','reciprocal_to_simple', ...
    'dense_to_simple','ffl_to_simple','complex_to_simple'};
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
        'VariableNames', {'case_id','patient_id','phase5','n_windows'});
    row.exact_window_count = sum(strcmp(winT.motif_census_mode(idx), 'exact'));
    row.sampled_window_count = sum(strcmp(winT.motif_census_mode(idx), 'sampled'));
    for m = 1:numel(metricNames)
        vals = winT.(metricNames{m})(idx);
        row.([metricNames{m} '_mean']) = mean(vals, 'omitnan');
        row.([metricNames{m} '_std']) = std(vals, 'omitnan');
        row.([metricNames{m} '_median']) = median(vals, 'omitnan');
    end
    stageT = [stageT; row];
end
end

