%#ok<*AGROW,*SAGROW>
function W = build_W_from_hyperedges(HE)

num_HE = numel(HE);

if num_HE < 2
    W = zeros(num_HE);
    return;
end

sizes = cellfun(@numel, HE);
sizes = sizes(:);

O = zeros(num_HE, num_HE);
for i = 1:num_HE
    for j = (i + 1):num_HE
        overlap_len = numel(intersect(HE{i}, HE{j}));
        if overlap_len > 0
            O(i, j) = overlap_len;
            O(j, i) = overlap_len;
        end
    end
end

D = sum(O > 0, 2);

Sizes_i = repmat(sizes, 1, num_HE);
Sizes_j = repmat(sizes', num_HE, 1);
D_i     = repmat(D, 1, num_HE);
D_j     = repmat(D', num_HE, 1);

Coverage = O ./ Sizes_i;

DegreeSum  = D_i + D_j;
DegreeBias = D_j ./ DegreeSum;
DegreeBias(DegreeSum == 0) = 0.5;

W = Coverage;

EqualSizeMask = (Sizes_i == Sizes_j) & (O > 0);
W(EqualSizeMask) = W(EqualSizeMask) .* DegreeBias(EqualSizeMask);

W(1:num_HE+1:end) = 0;

end
