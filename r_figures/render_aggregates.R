# Public aggregate figures: reads only reviewed group-level CSV tables.
# Run from the repository root. Base R is sufficient.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
root <- normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
out <- file.path(root, 'outputs', 'aggregate_figures')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
catalog <- read.csv(file.path(root, 'data', 'catalog.csv'), check.names=FALSE)
pick <- function(d, candidates) { z <- candidates[candidates %in% names(d)]; if(length(z)) z[1] else NA_character_ }
stage_order <- c('pre-ictal','early','mid','late','post-ictal','combined-ictal','ictal-combined')

# Each tuple is estimate, lower bound and upper bound in the original table.
summaries <- list(c('mean_z','low_z','high_z'), c('mean_z','ci_low_z','ci_high_z'),
 c('mean_delta','ci_low','ci_high'), c('mean_delta_S','ci_low','ci_high'),
 c('mean_difference','ci_low','ci_high'), c('residual_z','ci_low_z','ci_high_z'), c('residual','ci_low','ci_high'),
 c('norm_mean_residual','ci_low','ci_high'), c('centroid_distance_z','bootstrap_ci_low_z','bootstrap_ci_high_z'),
 c('rms_distance_z','rms_ci_low_z','rms_ci_high_z'), c('distance','ci_low','ci_high'),
 c('mean','ci_low','ci_high'), c('mean','low','high'), c('mean','lo','hi'),
 c('mean','sem','sem'), c('mean_occupancy_percent','sem','sem'), c('mean_fraction','sem_fraction','sem_fraction'),
 c('mean_diff_target_minus_preictal',NA,NA), c('mean_diff','sem_diff','sem_diff'),
 c('mean_diff_z',NA,NA))

render <- function(d, id, unit) {
 if(id=='Fig5_held_out') {
  v<-c(d$base,d$extended); barplot(v,names.arg=c('Baseline','DMWA extended'),col='white',border='#38566B',ylab='Mean held-out Brier score',main=id,ylim=c(0,max(v)*1.2))
  mtext(sprintf('Mean improvement %.4f; %d of %d participants',d$mean_delta,d$improved,d$n_patients),side=1,line=3,cex=.8);return()
 }
 if(id=='SFig6_balance') {
  d<-rbind(data.frame(label=paste(d$target_stage,d$covariate,'Before'),value=d$median_abs_before),data.frame(label=paste(d$target_stage,d$covariate,'After'),value=d$median_abs_after))
  par(mar=c(4,14,3,1));plot(d$value,seq_len(nrow(d)),yaxt='n',pch=21,bg='white',xlab='Median absolute standardized difference',ylab='',main=id);axis(2,at=seq_len(nrow(d)),labels=d$label,las=1,cex.axis=.65);return()
 }
 if(id=='SFig8_density_means') {
  d<-d[match(stage_order[stage_order %in% d$phase5],d$phase5),];matplot(seq_len(nrow(d)),d[c('density2','density3')],type='b',lty=1,pch=c(1,2),col=c('#2B7A78','#746080'),xaxt='n',xlab='',ylab='Mean density',main=id)
  axis(1,at=seq_len(nrow(d)),labels=d$phase5,cex.axis=.8);legend('topright',c('Pairs','Triplets'),col=c('#2B7A78','#746080'),pch=c(1,2),lty=1,bty='n');return()
 }
 spec<-NULL
 for(s in summaries) if(s[1] %in% names(d) && all(is.na(s[2:3]) | s[2:3] %in% names(d))) {spec<-s;break}
 if(is.null(spec)) stop(paste('No explicit summary mapping:',id))
 facet<-pick(d,c('metric','feature','motif','role','macro_state'))
 facets<-if(is.na(facet)) 'All' else unique(as.character(d[[facet]]))
 n<-length(facets);par(mfrow=c(ceiling(n/2),min(n,2)),mar=c(4,13,3,1),oma=c(1,0,2,0))
 for(f in facets) {
  z<-if(is.na(facet)) d else d[as.character(d[[facet]])==f,,drop=FALSE]
  fields<-intersect(c('target_stage','stage','phase5','label','clinical_class','cutoff','feature_set','method','contrast','reference','model','branch','profile','variant','ensemble','beta','strength','series'),names(z))
  fields<-fields[vapply(z[fields],function(v)length(unique(v[!is.na(v)]))>1,logical(1))]
  # Avoid redundant human-readable and internal branch labels.
  if('label' %in% fields) fields<-setdiff(fields,c('variant'))
  if('series' %in% fields) fields<-setdiff(fields,c('profile','model','branch'))
  label<-if(length(fields)) apply(z[fields],1,paste,collapse=' | ') else rep(f,nrow(z))
  value<-z[[spec[1]]];lo<-hi<-value;interval<-'point estimate'
  if(!is.na(spec[2])) {
   if(grepl('^sem',spec[2])) {lo<-value-z[[spec[2]]];hi<-value+z[[spec[3]]];interval<-'mean +/- SEM'}
   else {lo<-z[[spec[2]]];hi<-z[[spec[3]]];interval<-'reported 95% interval'}
  }
  good<-is.finite(value);value<-value[good];lo<-lo[good];hi<-hi[good];label<-label[good]
  yy<-rev(seq_along(value));lim<-range(c(value,lo,hi),finite=TRUE);pad<-diff(lim)*.08;if(pad==0)pad<-.05
  plot(value,yy,xlim=lim+c(-pad,pad),ylim=c(.4,max(yy)+.6),yaxt='n',pch=21,bg='white',col='#2B6470',xlab=paste(spec[1],paste0('(',interval,')')),ylab='',main=f,cex.main=.9,cex.lab=.75)
  segments(lo,yy,hi,yy,col='#66838B');points(value,yy,pch=21,bg='white',col='#2B6470',cex=.8)
  axis(2,at=yy,labels=label,las=1,tick=FALSE,cex.axis=max(.37,min(.75,15/length(yy))))
  if(lim[1]<0 && lim[2]>0) abline(v=0,lty=3,col='grey60')
 }
 mtext(paste(id,'- aggregate results'),outer=TRUE,line=.5,font=2)
}

for(i in seq_len(nrow(catalog))) {
 id<-catalog$dataset[i];d<-read.csv(file.path(root,catalog$file[i]),check.names=FALSE)
 height<-if(id=='SFig5_density') 16 else if(nrow(d)>70) 12 else 8
 pdf(file.path(out,paste0(id,'.pdf')),width=12,height=height,useDingbats=FALSE,title=paste(id,'aggregate results'),author='')
 render(d,id,catalog$analysis_unit[i]);dev.off()
 png(file.path(out,paste0(id,'.png')),width=2400,height=height*200,res=200)
 render(d,id,catalog$analysis_unit[i]);dev.off()
}
write.csv(catalog,file.path(out,'rendered_catalog.csv'),row.names=FALSE)
cat(sprintf('Rendered %d aggregate datasets as PDF and PNG.\n',nrow(catalog)))
