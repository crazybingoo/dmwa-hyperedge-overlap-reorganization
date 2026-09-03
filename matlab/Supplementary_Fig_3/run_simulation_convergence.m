%% Paper 1 simulation convergence and Monte Carlo uncertainty checks
clear; clc;
rng(20260714, 'twister');

projectRoot = getenv('PAPER1_PROJECT_ROOT');
if isempty(projectRoot), projectRoot = pwd; end
outDir = fullfile(projectRoot, 'stats', 'supplementary_sensitivity');
if ~exist(outDir, 'dir'), mkdir(outDir); end

models = ["density_only", "weight_strengthening", "reciprocity_only", ...
    "source_concentration", "overlap_synchronization", "role_reciprocity_hybrid"];
betas = (0:0.1:1)';
sizes = [30 45 60];
R = 80;
sampleTriples = 2500;

sizeOut = fullfile(outDir,'supp_fig3_size_sensitivity_raw.csv');
if ~exist(sizeOut, 'file')
    rows = cell(numel(models)*numel(betas)*numel(sizes)*R, 10);
    cursor = 0;
    for ni = 1:numel(sizes)
        n = sizes(ni);
        fprintf('Size sensitivity n=%d\n', n);
        for mi = 1:numel(models)
            for bi = 1:numel(betas)
                for r = 1:R
                    W = simulate_W(models(mi), betas(bi), n);
                    [Sval,gap,recFrac,denseFrac] = matrix_metrics(W,sampleTriples);
                    cursor = cursor + 1;
                    rows(cursor,:) = {char(models(mi)),betas(bi),n,r,sampleTriples,Sval,gap,recFrac,denseFrac,20260714};
                end
            end
        end
    end
    sizeT = cell2table(rows(1:cursor,:), 'VariableNames', {'model','beta','n_nodes','replicate', ...
        'sample_triples','S','source_sink_gini_gap','reciprocal_fraction','dense_fraction','seed'});
    for i=2:10
        name=sizeT.Properties.VariableNames{i};
        if iscell(sizeT.(name)), sizeT.(name)=cell2mat(sizeT.(name)); end
    end
    writetable(sizeT, sizeOut);
end

fit = readtable(fullfile(projectRoot,'stats','phase15_simulation_grid_observed_fit.csv'),'TextType','string');
Kgrid = [250 500 1000 2500 5000];
n = 45; Rmc = 120;
mcRows = cell(numel(models)*Rmc*numel(Kgrid), 11); cursor=0;
for mi=1:numel(models)
    fitModel = fit(fit.model==models(mi), :);
    [~, bestIdx] = min(fitModel.scaled_euclidean_distance);
    bestBeta = fitModel.beta(bestIdx);
    fprintf('Triad convergence %s beta=%.2f\n', models(mi), bestBeta);
    for r=1:Rmc
        W=simulate_W(models(mi),bestBeta,n); B=W>0;
        [exactRec,exactDense]=motif_exact_fractions(B);
        for ki=1:numel(Kgrid)
            rng(20260714 + mi*100000 + r*100 + ki,'twister');
            [sampleRec,sampleDense]=motif_sample_fractions(B,Kgrid(ki));
            cursor=cursor+1;
            mcRows(cursor,:)={char(models(mi)),bestBeta,n,r,Kgrid(ki),exactRec,sampleRec, ...
                abs(sampleRec-exactRec),exactDense,sampleDense,abs(sampleDense-exactDense)};
        end
    end
end
mcT=cell2table(mcRows(1:cursor,:), 'VariableNames', {'model','best_beta','n_nodes','replicate', ...
    'sample_triples','exact_reciprocal_fraction','sample_reciprocal_fraction','abs_error_reciprocal', ...
    'exact_dense_fraction','sample_dense_fraction','abs_error_dense'});
for i=2:11
    name=mcT.Properties.VariableNames{i};
    if iscell(mcT.(name)), mcT.(name)=cell2mat(mcT.(name)); end
end
writetable(mcT, fullfile(outDir,'supp_fig3_triad_sampling_convergence.csv'));
fprintf('Simulation convergence checks completed.\n');

function W = simulate_W(model,beta,n)
p0=0.045; baseMask=rand(n)<p0; baseMask(1:n+1:end)=false; W=zeros(n);
W(baseMask)=gamrnd(2,0.5,nnz(baseMask),1);
sourceCore=1:min(4,n); syncPool=1:round(0.60*n); sinkPool=1:n;
switch model
    case "density_only"
        pAdd=0.12*beta; M=rand(n)<pAdd; M(1:n+1:end)=false; W(M)=W(M)+gamrnd(2,0.5,nnz(M),1);
    case "weight_strengthening"
        W=W.*(1+3*beta);
    case "reciprocity_only"
        for i=1:n, for j=i+1:n, if rand<0.16*beta
                    w=gamrnd(2,0.5); W(i,j)=W(i,j)+w; W(j,i)=W(j,i)+w*(0.8+0.4*rand);
                end, end, end
    case "source_concentration"
        for s=sourceCore, for t=sinkPool, if s~=t && rand<0.28*beta
                    W(s,t)=W(s,t)+gamrnd(2.5,0.5);
                end, end, end
    case "overlap_synchronization"
        for ii=1:numel(syncPool), for jj=ii+1:numel(syncPool)
                i=syncPool(ii); j=syncPool(jj); if rand<0.075*beta
                    w=gamrnd(1.5,0.35); W(i,j)=W(i,j)+w; W(j,i)=W(j,i)+w*(0.85+0.3*rand);
                end
            end, end
        for s=sourceCore, for t=sinkPool, if s~=t && rand<0.26*beta
                    W(s,t)=W(s,t)+gamrnd(2.8,0.45);
                end, end, end
    case "role_reciprocity_hybrid"
        M=rand(n)<0.025*beta; M(1:n+1:end)=false; W(M)=W(M)+gamrnd(1.5,0.35,nnz(M),1);
        for ii=1:numel(syncPool), for jj=ii+1:numel(syncPool)
                i=syncPool(ii); j=syncPool(jj); if rand<0.026*beta
                    w=gamrnd(1.4,0.28); W(i,j)=W(i,j)+w; W(j,i)=W(j,i)+w*(0.8+0.35*rand);
                end
            end, end
        for s=sourceCore, for t=sinkPool, if s~=t && rand<0.34*beta
                    W(s,t)=W(s,t)+gamrnd(2.8,0.42);
                end, end, end
end
W(1:n+1:end)=0;
end

function [Sval,gap,recFrac,denseFrac]=matrix_metrics(W,k)
n=size(W,1); mask=triu(true(n),1); a=W(mask); WT=W'; b=WT(mask); den=sum(max(a,b));
if den==0, Sval=0; else, Sval=sum(min(a,b))/den; end
sout=sum(W,2); sinv=sum(W,1)'; q=(sout-sinv)./(sout+sinv+eps);
gap=gini_all(max(q,0))-gini_all(max(-q,0));
[recFrac,denseFrac]=motif_sample_fractions(W>0,k);
end

function g=gini_all(x)
x=sort(max(0,x(:))); n=numel(x); if isempty(x)||sum(x)==0, g=0; return; end
g=2*sum((1:n)'.*x)/(n*sum(x))-(n+1)/n;
end

function [recFrac,denseFrac]=motif_sample_fractions(B,k)
n=size(B,1); total=nchoosek(n,3); k=min(k,total);
if k==total, triads=nchoosek(1:n,3); else, triads=sample_unique_triads(n,k); end
[recFrac,denseFrac]=fractions_on_triads(B,triads);
end

function [recFrac,denseFrac]=motif_exact_fractions(B)
triads=nchoosek(1:size(B,1),3); [recFrac,denseFrac]=fractions_on_triads(B,triads);
end

function triads=sample_unique_triads(n,k)
triads=zeros(0,3);
while size(triads,1)<k
    cand=sort(randi(n,max(1000,2*(k-size(triads,1))),3),2);
    cand=cand(cand(:,1)~=cand(:,2)&cand(:,2)~=cand(:,3),:);
    triads=unique([triads;cand],'rows','stable'); %#ok<AGROW>
end
triads=triads(1:k,:);
end

function [recFrac,denseFrac]=fractions_on_triads(B,triads)
n=size(B,1); i=triads(:,1); j=triads(:,2); k=triads(:,3);
eij=B(sub2ind([n n],i,j)); eji=B(sub2ind([n n],j,i));
eik=B(sub2ind([n n],i,k)); eki=B(sub2ind([n n],k,i));
ejk=B(sub2ind([n n],j,k)); ekj=B(sub2ind([n n],k,j));
uij=eij|eji; uik=eik|eki; ujk=ejk|ekj;
connected=(uij+uik+ujk)>=2 & (uij|uik) & (uij|ujk) & (uik|ujk);
if ~any(connected), recFrac=0; denseFrac=0; return; end
rec=(eij&eji)|(eik&eki)|(ejk&ekj);
directed=double(eij)+double(eji)+double(eik)+double(eki)+double(ejk)+double(ekj);
recFrac=mean(rec(connected)); denseFrac=mean(directed(connected)>=4);
end
