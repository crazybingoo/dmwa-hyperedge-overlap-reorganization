%#ok<*AGROW,*SAGROW>
%% Phase 16: elbow-threshold macro-state reconstruction for Paper 1 Fig. 2
% Recomputes the Fig. 2 SOZ/PZ/NIZ macro-state occupancy from SEEG
% windows using the seizure-specific threshold-density elbow rule. The
% macro-state definition and possible-combination normalization follow the
% legacy W6.m logic, but hyperedges are re-extracted under the elbow threshold.

clear; clc;

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot)
    projectRoot = pwd;
end
dataDir = fullfile(projectRoot, 'data_tables');
figDir = fullfile(projectRoot, 'figures', 'Fig_2');
noteDir = fullfile(projectRoot, 'notes');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(noteDir, 'dir'), mkdir(noteDir); end

caseMetaPath = fullfile(dataDir, 'paper1_case_metadata.csv');
thresholdPath = fullfile(dataDir, 'paper2_threshold_elbow_case_summary.csv');
caseMeta = readtable(caseMetaPath, 'TextType', 'string');
thresholdSummary = readtable(thresholdPath, 'TextType', 'string');
requiredMeta = ["case_id","file_path","node_region_path","pre_start","pre_end", ...
    "ictal_start","ictal_end","post_start","post_end"];
missingMeta = setdiff(requiredMeta, string(caseMeta.Properties.VariableNames));
if ~isempty(missingMeta)
    error('Case metadata is missing required columns: %s', strjoin(missingMeta, ', '));
end

winLenSec = 3;
stageOrder = ["pre-ictal","early","mid","late","post-ictal"];
macroRaw = ["SOZ / SOZ+PZ", ...
    "SOZ+NIZ / SOZ+PZ+NIZ", ...
    "PZ+NIZ", "PZ", "NIZ"];
macroDisplay = ["NIZ", "PZ", "PZ+NIZ", ...
    "SOZ+NIZ / SOZ+PZ+NIZ", "SOZ / SOZ+PZ"];

windowRows = {};
sourceRows = {};
rowIdx = 0;
sourceIdx = 0;

fprintf('Phase16 elbow macro-state reconstruction for %d cases.\n', height(caseMeta));

for c = 1:height(caseMeta)
    caseId = string(caseMeta.case_id(c));
    patientId = get_patient_id(caseMeta, c);
    filePath = resolve_local_path(char(caseMeta.file_path(c)), projectRoot);
    nodeRegionPath = resolve_local_path(char(caseMeta.node_region_path(c)), projectRoot);

    match = thresholdSummary(thresholdSummary.case_id == caseId, :);
    if height(match) ~= 1
        error('Cannot find exactly one elbow threshold for %s', caseId);
    end
    threshold = match.selected_threshold(1);

    if ~exist(nodeRegionPath, 'file')
        error('Missing node-region file for %s: %s', caseId, nodeRegionPath);
    end
    L = load(nodeRegionPath, 'node_regions');
    if ~isfield(L, 'node_regions')
        error('Missing node_regions in %s', nodeRegionPath);
    end
    nodeRegions = double(L.node_regions(:));

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

    nChannel = size(X1, 1);
    nSample = size(X1, 2);
    nWindow = floor(nSample / fs) - winLenSec + 1;
    if numel(nodeRegions) < nChannel
        error('%s has %d channels but only %d node-region labels.', caseId, nChannel, numel(nodeRegions));
    end
    nodeRegions = nodeRegions(1:nChannel);

    stageLabels = assign_stage_labels(nWindow, ...
        caseMeta.pre_start(c):caseMeta.pre_end(c), ...
        caseMeta.ictal_start(c):caseMeta.ictal_end(c), ...
        caseMeta.post_start(c):caseMeta.post_end(c));

    nS = sum(nodeRegions == 1);
    nP = sum(nodeRegions == 2);
    nN = sum(nodeRegions == 3);
    rawCounts = zeros(nWindow, 5);
    normCounts = zeros(nWindow, 5);
    totalHyperedges = zeros(nWindow, 1);
    unclassifiedCounts = zeros(nWindow, 1);
    pairDensity = zeros(nWindow, 1);
    meanPlv = zeros(nWindow, 1);

    fprintf('[%02d/%02d] %s, threshold %.4f\n', c, height(caseMeta), caseId, threshold);

    for w = 1:nWindow
        startPoint = fs * (w - 1) + 1;
        endPoint = fs * (w - 1 + winLenSec);
        datanew = X1(:, startPoint:endPoint);

        [hyperEdges, pairDensity(w), meanPlv(w)] = extract_hyperedges_at_threshold(datanew, threshold);
        totalHyperedges(w) = numel(hyperEdges);

        for h = 1:numel(hyperEdges)
            nodes = hyperEdges{h};
            nodes = unique(nodes(nodes > 0 & nodes <= numel(nodeRegions)));
            if isempty(nodes)
                unclassifiedCounts(w) = unclassifiedCounts(w) + 1;
                continue;
            end
            edgeSize = numel(nodes);
            regions = sort(unique(nodeRegions(nodes)))';
            origLabel = original_region_label(regions);
            if origLabel == 0
                unclassifiedCounts(w) = unclassifiedCounts(w) + 1;
                continue;
            end
            macroLabel = map_original_to_macro(origLabel);
            rawCounts(w, macroLabel) = rawCounts(w, macroLabel) + 1;

            possibleCount = get_possible_combinations(edgeSize, origLabel, nS, nP, nN);
            if possibleCount > 0
                normCounts(w, macroLabel) = normCounts(w, macroLabel) + 1 / possibleCount;
            end
        end
    end

    denomRaw = sum(rawCounts, 2);
    denomRaw(denomRaw == 0) = 1;
    pctRaw = rawCounts ./ denomRaw * 100;

    denomNorm = sum(normCounts, 2);
    denomNorm(denomNorm == 0) = 1;
    pctNorm = normCounts ./ denomNorm * 100;
    pctNormReversed = pctNorm(:, end:-1:1);
    pctNormReversed = movmean(pctNormReversed, 5, 1);

    for w = 1:nWindow
        phase5 = string(stageLabels{w});
        if phase5 == "unassigned"
            continue;
        end
        for mi = 1:5
            rowIdx = rowIdx + 1;
            rawIdx = 6 - mi;
            windowRows(rowIdx, :) = {char(caseId), char(patientId), w, char(phase5), ...
                threshold, nChannel, pairDensity(w), meanPlv(w), totalHyperedges(w), ...
                unclassifiedCounts(w), mi, char(macroDisplay(mi)), ...
                rawCounts(w, rawIdx), normCounts(w, rawIdx), pctRaw(w, rawIdx), ...
                pctNormReversed(w, mi)};
        end
    end

    sourceIdx = sourceIdx + 1;
    sourceRows(sourceIdx, :) = {char(caseId), char(patientId), '<local SEEG input file>', '<local node-region input file>', ...
        nWindow, threshold, nChannel, nS, nP, nN, ...
        sum(stageLabels == "pre-ictal"), sum(stageLabels == "early"), ...
        sum(stageLabels == "mid"), sum(stageLabels == "late"), sum(stageLabels == "post-ictal")};
end

winVars = {'case_id','patient_id','window_idx','phase5','selected_threshold', ...
    'n_channels','pair_density_elbow','mean_plv','n_hyperedges','unclassified_count', ...
    'macro_state_index','macro_state','raw_count','normalized_count', ...
    'raw_occupancy_percent','occupancy_percent'};
windowTable = cell2table(windowRows, 'VariableNames', winVars);
windowTable = convert_numeric_columns(windowTable, winVars(3:end), {'phase5','macro_state'});

sourceVars = {'case_id','patient_id','source_file','node_region_file', ...
    'n_windows','selected_threshold','n_channels','n_soz','n_pz','n_niz', ...
    'n_preictal','n_early','n_mid','n_late','n_postictal'};
sourceTable = cell2table(sourceRows, 'VariableNames', sourceVars);
sourceTable = convert_numeric_columns(sourceTable, sourceVars(5:end), {});

patientRows = {};
prIdx = 0;
caseIds = unique(windowTable.case_id, 'stable');
for c = 1:numel(caseIds)
    caseId = string(caseIds{c});
    caseWin = windowTable(strcmp(windowTable.case_id, caseId), :);
    patientId = string(caseWin.patient_id{1});
    for s = 1:numel(stageOrder)
        phase5 = stageOrder(s);
        for mi = 1:5
            idx = strcmp(caseWin.phase5, phase5) & caseWin.macro_state_index == mi;
            if ~any(idx), continue; end
            vals = caseWin.occupancy_percent(idx);
            prIdx = prIdx + 1;
            patientRows(prIdx, :) = {char(caseId), char(patientId), char(phase5), ...
                mi, char(macroDisplay(mi)), sum(idx), mean(vals, 'omitnan')};
        end
    end
end
patientVars = {'case_id','patient_id','stage','macro_state_index','macro_state', ...
    'n_windows','occupancy_percent'};
patientTable = cell2table(patientRows, 'VariableNames', patientVars);
patientTable = convert_numeric_columns(patientTable, {'macro_state_index','n_windows','occupancy_percent'}, {'stage','macro_state'});

summaryTable = summarize_stage(patientTable, stageOrder, macroDisplay);
contrastTable = paired_contrasts(patientTable, stageOrder(2:end), macroDisplay);

sourceOut = fullfile(figDir, 'fig2_macro_state_source_files_elbow.csv');
windowOut = fullfile(figDir, 'fig2_macro_state_window_occupancy_elbow.csv');
patientOut = fullfile(figDir, 'fig2_macro_state_patient_occupancy_elbow.csv');
summaryOut = fullfile(figDir, 'fig2_macro_state_summary_elbow.csv');
contrastOut = fullfile(figDir, 'fig2_macro_state_paired_contrasts_elbow.csv');
noteOut = fullfile(noteDir, 'phase16_elbow_macro_state_reconstruction.md');

writetable(sourceTable, sourceOut);
writetable(windowTable, windowOut);
writetable(patientTable, patientOut);
writetable(summaryTable, summaryOut);
writetable(contrastTable, contrastOut);

fid = fopen(noteOut, 'w');
fprintf(fid, '# Phase16 elbow-threshold macro-state reconstruction\n\n');
fprintf(fid, '- Date: 2026-06-17\n');
fprintf(fid, '- Cases: %d\n', height(caseMeta));
fprintf(fid, '- Window length: %d s\n', winLenSec);
fprintf(fid, '- Sampling rate: read `fs` from `.mat` when present; otherwise fallback to 1024 Hz.\n');
fprintf(fid, '- Threshold rule: one seizure-specific threshold-density elbow applied to all windows from that seizure.\n');
fprintf(fid, '- Macro-state rule: legacy W6.m SOZ/PZ/NIZ region-combination labels and possible-combination normalization.\n');
fprintf(fid, '- Display order: NIZ, PZ, PZ+NIZ, SOZ+NIZ / SOZ+PZ+NIZ, SOZ / SOZ+PZ.\n');
fprintf(fid, '- Moving average: 5-window moving average after reversing the five normalized macro-state percentages, matching legacy Fig. 2.\n');
fprintf(fid, '- Window rows: %d\n', height(windowTable));
fprintf(fid, '- Patient-stage rows: %d\n', height(patientTable));
fprintf(fid, '\n## Outputs\n\n');
fprintf(fid, '- `%s`\n', sourceOut);
fprintf(fid, '- `%s`\n', windowOut);
fprintf(fid, '- `%s`\n', patientOut);
fprintf(fid, '- `%s`\n', summaryOut);
fprintf(fid, '- `%s`\n', contrastOut);
fclose(fid);

fprintf('Wrote Phase16 Fig. 2 elbow macro-state outputs to %s\n', figDir);
disp(contrastTable);

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

function patientId = get_patient_id(caseMeta, rowIndex)
hasPatient = ismember("patient_id", string(caseMeta.Properties.VariableNames));
if hasPatient && strlength(caseMeta.patient_id(rowIndex)) > 0
    patientId = char(caseMeta.patient_id(rowIndex));
    return;
end
caseId = char(caseMeta.case_id(rowIndex));
idx = strfind(caseId, '_cut');
if ~isempty(idx)
    patientId = caseId(1:idx(1)-1);
elseif endsWith(caseId, '_process')
    patientId = erase(string(caseId), "_process");
    patientId = char(patientId);
else
    patientId = caseId;
end
end

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

function label = original_region_label(regions)
label = 0;
if isequal(regions, 1)
    label = 1;
elseif isequal(regions, 2)
    label = 2;
elseif isequal(regions, 3)
    label = 3;
elseif isequal(regions, [1 2])
    label = 4;
elseif isequal(regions, [1 3])
    label = 5;
elseif isequal(regions, [2 3])
    label = 6;
elseif isequal(regions, [1 2 3])
    label = 7;
end
end

function macroLabel = map_original_to_macro(origLabel)
map = [1, 4, 5, 1, 2, 3, 2];
macroLabel = map(origLabel);
end

function possibleCount = get_possible_combinations(edgeSize, labelId, nS, nP, nN)
possibleCount = 0;
switch labelId
    case 1
        possibleCount = safe_nchoosek(nS, edgeSize);
    case 2
        possibleCount = safe_nchoosek(nP, edgeSize);
    case 3
        possibleCount = safe_nchoosek(nN, edgeSize);
    case 4
        for a = 1:(edgeSize - 1)
            b = edgeSize - a;
            possibleCount = possibleCount + safe_nchoosek(nS, a) * safe_nchoosek(nP, b);
        end
    case 5
        for a = 1:(edgeSize - 1)
            b = edgeSize - a;
            possibleCount = possibleCount + safe_nchoosek(nS, a) * safe_nchoosek(nN, b);
        end
    case 6
        for a = 1:(edgeSize - 1)
            b = edgeSize - a;
            possibleCount = possibleCount + safe_nchoosek(nP, a) * safe_nchoosek(nN, b);
        end
    case 7
        for a = 1:(edgeSize - 2)
            for b = 1:(edgeSize - a - 1)
                c = edgeSize - a - b;
                possibleCount = possibleCount + safe_nchoosek(nS, a) * safe_nchoosek(nP, b) * safe_nchoosek(nN, c);
            end
        end
end
end

function val = safe_nchoosek(n, k)
if k < 0 || n < 0 || k > n
    val = 0;
else
    val = nchoosek(n, k);
end
end

function T = convert_numeric_columns(T, names, skipNames)
for v = 1:numel(names)
    name = char(names{v});
    if any(strcmp(skipNames, name))
        continue;
    end
    if iscell(T.(name))
        T.(name) = cell2mat(T.(name));
    end
end
end

function stage = assign_stage_labels(nWindow, preIdx, ictalIdx, postIdx)
stage = strings(nWindow, 1);
stage(:) = "unassigned";
preIdx = sanitize_idx(preIdx, nWindow);
ictalIdx = sanitize_idx(ictalIdx, nWindow);
postIdx = sanitize_idx(postIdx, nWindow);
ictalIdx = setdiff(ictalIdx, preIdx, 'stable');
postIdx = setdiff(postIdx, union(preIdx, ictalIdx), 'stable');
splitPts = round(linspace(0, numel(ictalIdx), 4));
earlyIdx = ictalIdx(1:splitPts(2));
midIdx = ictalIdx(splitPts(2)+1:splitPts(3));
lateIdx = ictalIdx(splitPts(3)+1:end);
stage(preIdx) = "pre-ictal";
stage(earlyIdx) = "early";
stage(midIdx) = "mid";
stage(lateIdx) = "late";
stage(postIdx) = "post-ictal";
end

function idx = sanitize_idx(idx, nWindow)
idx = idx(:)';
idx = idx(idx >= 1 & idx <= nWindow);
idx = unique(idx, 'stable');
end

function summaryTable = summarize_stage(patientTable, stageOrder, macroDisplay)
rows = {};
ri = 0;
for s = 1:numel(stageOrder)
    for mi = 1:5
        idx = strcmp(patientTable.stage, stageOrder(s)) & patientTable.macro_state_index == mi;
        vals = patientTable.occupancy_percent(idx);
        vals = vals(isfinite(vals));
        n = numel(vals);
        if n == 0
            stats = nan(1, 8);
        else
            stats = [mean(vals,'omitnan'), median(vals,'omitnan'), std(vals,'omitnan'), ...
                std(vals,'omitnan') / sqrt(n), prctile(vals,25), prctile(vals,75), ...
                prctile(vals,2.5), prctile(vals,97.5)];
        end
        ri = ri + 1;
        rows(ri, :) = {char(stageOrder(s)), mi, char(macroDisplay(mi)), n, ...
            stats(1), stats(2), stats(3), stats(4), stats(5), stats(6), stats(7), stats(8)};
    end
end
vars = {'stage','macro_state_index','macro_state','n_cases','mean_occupancy_percent', ...
    'median_occupancy_percent','sd','sem','q25','q75','ci2_5','ci97_5'};
summaryTable = cell2table(rows, 'VariableNames', vars);
summaryTable = convert_numeric_columns(summaryTable, vars(2:end), {'macro_state'});
end

function contrastTable = paired_contrasts(patientTable, targetStages, macroDisplay)
rows = {};
ri = 0;
rawP = nan(numel(targetStages) * 5, 1);
for s = 1:numel(targetStages)
    target = targetStages(s);
    for mi = 1:5
        pre = patientTable(strcmp(patientTable.stage, "pre-ictal") & patientTable.macro_state_index == mi, :);
        tar = patientTable(strcmp(patientTable.stage, target) & patientTable.macro_state_index == mi, :);
        joined = innerjoin(pre, tar, 'Keys', 'case_id', ...
            'LeftVariables', {'case_id','patient_id','occupancy_percent'}, ...
            'RightVariables', {'occupancy_percent'});
        diffVals = joined.occupancy_percent_tar - joined.occupancy_percent_pre;
        try
            pVal = signrank(diffVals);
        catch
            pVal = NaN;
        end
        ri = ri + 1;
        rawP(ri) = pVal;
        rows(ri, :) = {char(target), mi, char(macroDisplay(mi)), height(joined), ...
            mean(joined.occupancy_percent_pre, 'omitnan'), ...
            mean(joined.occupancy_percent_tar, 'omitnan'), ...
            mean(diffVals, 'omitnan'), median(diffVals, 'omitnan'), ...
            sum(diffVals > 0) / max(numel(diffVals), 1), pVal};
    end
end
vars = {'target_stage','macro_state_index','macro_state','n_cases','preictal_mean', ...
    'target_mean','mean_diff_target_minus_preictal','median_diff_target_minus_preictal', ...
    'positive_fraction','signrank_p'};
contrastTable = cell2table(rows, 'VariableNames', vars);
contrastTable = convert_numeric_columns(contrastTable, vars(2:end), {'macro_state'});
contrastTable.signrank_p_fdr_bh = benjamini_hochberg(rawP);
end

function adjusted = benjamini_hochberg(pValues)
p = pValues(:);
adjusted = nan(size(p));
valid = ~isnan(p);
pv = p(valid);
m = numel(pv);
if m == 0
    return;
end
[sortedP, sortIdx] = sort(pv, 'ascend');
adjSorted = sortedP .* m ./ (1:m)';
for ii = m-1:-1:1
    adjSorted(ii) = min(adjSorted(ii), adjSorted(ii+1));
end
adjSorted = min(adjSorted, 1);
tmp = nan(m, 1);
tmp(sortIdx) = adjSorted;
adjusted(valid) = tmp;
end

