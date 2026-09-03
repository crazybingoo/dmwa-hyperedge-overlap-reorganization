%#ok<*AGROW,*SAGROW>
function models = compare_ablation_models(hyperEdges)
%COMPARE_ABLAhION_MODELS Build DMWA ablation variants for one window.
%
% The paper compares the full DMWA weighting rule against two controls:
%   1. Coverage-only: retains overlap coverage O_ij / |e_i|.
%   2. degree-bias-only: retains pairwise degree bias among overlapping
%      hyperedges.
%
% No case metadata or patient data are embedded in this function.

numEdges = numel(hyperEdges);
models = struct();
models.W_full = zeros(numEdges);
models.W_coverage_only = zeros(numEdges);
models.W_degree_bias_only = zeros(numEdges);
models.S_full = 0;
models.S_coverage_only = 0;
models.S_degree_bias_only = 0;

if numEdges < 2
    return;
end

sizes = cellfun(@numel, hyperEdges);
sizes = sizes(:);
overlap = zeros(numEdges, numEdges);

for i = 1:numEdges
    for j = (i + 1):numEdges
        overlapLen = numel(intersect(hyperEdges{i}, hyperEdges{j}));
        if overlapLen > 0
            overlap(i, j) = overlapLen;
            overlap(j, i) = overlapLen;
        end
    end
end

degree = sum(overlap > 0, 2);
sizesI = repmat(sizes, 1, numEdges);
sizesJ = repmat(sizes.', numEdges, 1);
degreeI = repmat(degree, 1, numEdges);
degreeJ = repmat(degree.', numEdges, 1);

coverage = overlap ./ sizesI;
degreeSum = degreeI + degreeJ;
degreeBias = degreeJ ./ degreeSum;
degreeBias(degreeSum == 0) = 0.5;

WFull = coverage;
sameOrder = (sizesI == sizesJ) & (overlap > 0);
WFull(sameOrder) = WFull(sameOrder) .* degreeBias(sameOrder);
WFull(overlap == 0) = 0;
WFull(1:numEdges+1:end) = 0;

WCoverage = coverage;
WCoverage(overlap == 0) = 0;
WCoverage(1:numEdges+1:end) = 0;

WDegree = degreeBias;
WDegree(overlap == 0) = 0;
WDegree(1:numEdges+1:end) = 0;

models.W_full = WFull;
models.W_coverage_only = WCoverage;
models.W_degree_bias_only = WDegree;
models.S_full = calculate_S_from_W(WFull);
models.S_coverage_only = calculate_S_from_W(WCoverage);
models.S_degree_bias_only = calculate_S_from_W(WDegree);
end

function S = calculate_S_from_W(W)
n = size(W, 1);
if n < 2
    S = 0;
    return;
end
W = double(W);
W(~isfinite(W)) = 0;
W(W < 0) = 0;
W(1:n+1:end) = 0;
mask = triu(true(n), 1);
Wt = W.';
wij = W(mask);
wji = Wt(mask);
denominator = sum(max(wij, wji));
if denominator > 0
    S = sum(min(wij, wji)) / denominator;
else
    S = 0;
end
end
