source('scripts/common.R')
suppressPackageStartupMessages(library(phylolm))
OUT <- file.path(RUN,'figures','figure_10'); ANALYSIS <- file.path(OUT,'analysis'); SNAP <- file.path(RUN,'inputs','nuclear_allometry')
dir.create(ANALYSIS,recursive=TRUE,showWarnings=FALSE); dir.create(SNAP,recursive=TRUE,showWarnings=FALSE)
repository <- normalizePath(file.path(RUN,'..','..'))
source_files <- c(paired_objects=file.path(DATA,'Supplementary_Data_S18_frozen_cell_nucleus_objects.csv'),species_estimates=file.path(DATA,'Supplementary_Data_S19_cell_nucleus_species_estimates.csv'),prior_pairwise=file.path(DATA,'Supplementary_Data_S26_pairwise_pgls_sensitivity.csv'),identity=file.path(repository,'data','identity','species_taxonomy_crosswalk.csv'),tree=file.path(DATA,'trees','desmo900dated_source46.tre'))
sha256 <- function(path) substr(system2('sha256sum',shQuote(path),stdout=TRUE),1,64)
manifest <- bind_rows(lapply(names(source_files),function(id){
  source <- source_files[[id]]; destination <- file.path(SNAP,basename(source))
  if(!file.exists(destination)) stopifnot(file.copy(source,destination))
  stopifnot(identical(sha256(source),sha256(destination)))
  tibble(input_id=id,source=source,snapshot=destination,sha256=sha256(source))
}))
write_csv(manifest,file.path(SNAP,'SOURCE_MANIFEST.csv'))
get_input <- function(id) manifest$snapshot[match(id,manifest$input_id)]
objects <- read_csv(get_input('paired_objects'),show_col_types=FALSE)
released <- read_csv(get_input('species_estimates'),show_col_types=FALSE)
identity <- read_csv(get_input('identity'),show_col_types=FALSE) %>% transmute(species=current_species,morphology_label=species_label(current_species))
stopifnot(!anyDuplicated(identity$species),!anyDuplicated(identity$morphology_label))
identify <- function(data){
  names(data)[names(data)=='species'] <- 'morphology_label'
  data <- left_join(data,identity,by='morphology_label',relationship='many-to-one')
  stopifnot(!anyNA(data$species)); data
}
objects <- identify(objects); released <- identify(released)
stopifnot(nrow(objects)==1152L,n_distinct(objects$species)==24L,n_distinct(objects$specimen_id)==42L,nrow(released)==24L,
  !anyDuplicated(objects$cell_object_id),!anyDuplicated(objects$nucleus_object_id),
  all(is.finite(objects$cell_area_um2)),all(is.finite(objects$nuc_area_um2)),
  all(objects$cell_area_um2>objects$nuc_area_um2),all(objects$nuc_area_um2>0),
  max(abs(objects$nc_area_ratio-objects$nuc_area_um2/objects$cell_area_um2))<5.1e-7,
  max(abs(objects$cytoplasm_area_um2-(objects$cell_area_um2-objects$nuc_area_um2)))<1e-5,
  setequal(objects$species,released$species))
# Each fraction comes from an actual nucleus-cell pair, never a ratio of aggregate medians.
objects <- objects %>% mutate(nuclear_fraction=nuc_area_um2/cell_area_um2,paired_cytoplasm_um2=cell_area_um2-nuc_area_um2,
  development=if_else(species %in% c('aeneus','wrighti'),'Direct development','Aquatic larvae'))
animals <- objects %>% group_by(species,morphology_label,specimen_id,development) %>%
  summarise(n_paired_cells=n(),cell_area_um2=median(cell_area_um2),nucleus_area_um2=median(nuc_area_um2),
    nuclear_fraction=median(nuclear_fraction),cytoplasm_area_um2=median(paired_cytoplasm_um2),.groups='drop')
balanced <- animals %>% group_by(species,morphology_label,development) %>%
  summarise(n_animals=n(),n_paired_cells=sum(n_paired_cells),across(c(cell_area_um2,nucleus_area_um2,nuclear_fraction,cytoplasm_area_um2),mean),.groups='drop')
pooled <- objects %>% group_by(species) %>% summarise(cell=median(cell_area_um2),nucleus=median(nuc_area_um2),.groups='drop')
reproduction <- inner_join(released,pooled,by='species',relationship='one-to-one')
stopifnot(max(abs(reproduction$cell_area_um2-reproduction$cell))<1e-7,
  max(abs(reproduction$nucleus_area_um2-reproduction$nucleus))<1e-7,
  sum(balanced$n_animals)==42L,sum(balanced$n_paired_cells)==1152L,sum(balanced$n_animals==1)==7L,
  all(balanced$nuclear_fraction>0 & balanced$nuclear_fraction<1))
primary <- released %>% mutate(development=if_else(species %in% c('aeneus','wrighti'),'Direct development','Aquatic larvae'))
fit_allometry <- function(data,id){
  d <- as.data.frame(data); d$log_cell <- log10(d$cell_area_um2); d$log_nucleus <- log10(d$nucleus_area_um2); rownames(d)<-d$species
  caught <- character()
  fit <- withCallingHandlers(phylolm(log_nucleus~log_cell,data=d,phy=prune_source_tree(d$species),model='lambda',lower.bound=1e-7,upper.bound=1),
    warning=function(w){caught<<-c(caught,conditionMessage(w));invokeRestart('muffleWarning')})
  co <- summary(fit)$coefficients; slope <- unname(co['log_cell','Estimate']); se <- unname(co['log_cell','StdErr']); df <- nrow(d)-2L; critical <- qt(.975,df)
  stat <- tibble(analysis_id=id,n_species=nrow(d),model='Pagel_lambda_PGLS',formula='log10(nucleus_area_um2) ~ log10(cell_area_um2)',
    coefficient_scale='unstandardized_log10_log10',intercept=unname(coef(fit)[1]),slope=slope,standard_error=se,
    ci_low=slope-critical*se,ci_high=slope+critical*se,df=df,t_vs_zero=slope/se,p_vs_zero=2*pt(-abs(slope/se),df),
    t_vs_isometry=(slope-1)/se,p_vs_isometry=2*pt(-abs((slope-1)/se),df),lambda=fit$optpar,
    lambda_at_boundary=fit$optpar<=1.1e-7|fit$optpar>=1-1e-7,
    interval_method='95% t coefficient interval conditional on fitted lambda',measurement_error_corrected=FALSE,warnings=paste(caught,collapse=' | '))
  stopifnot(abs(stat$t_vs_zero-co['log_cell','t.value'])<1e-8,abs(stat$p_vs_zero-co['log_cell','p.value'])<1e-8,
    stat$ci_low<=stat$slope,stat$ci_high>=stat$slope)
  list(fit=fit,summary=stat)
}
fit_inputs <- list(released24=primary,equal_animal24=balanced,released22=filter(primary,!species %in% c('aeneus','wrighti')),equal_animal22=filter(balanced,!species %in% c('aeneus','wrighti')))
fits <- lapply(names(fit_inputs),function(id)fit_allometry(fit_inputs[[id]],id)); names(fits)<-names(fit_inputs)
models <- bind_rows(lapply(fits,`[[`,'summary'))
prior <- read_csv(get_input('prior_pairwise'),show_col_types=FALSE) %>% filter(comparison=='nucleus_area_vs_cell_area')
# Figure 5 used the reverse conditional regression; reproduce standardized association, not raw slope.
primary_standardized_beta <- models$slope[1]*sd(log10(primary$cell_area_um2))/sd(log10(primary$nucleus_area_um2))
stopifnot(nrow(prior)==1L,abs(primary_standardized_beta-prior$pgls_standardized_beta)<1e-6,
  abs(models$p_vs_zero[1]-prior$pgls_p_value)<1e-7,identical(models$n_species,c(24L,24L,22L,22L)))
# Resample animals, then paired cells in each sampled animal, with equal animal weights.
# Single-animal species cannot identify between-animal uncertainty.
set.seed(20260924); B<-1000L; species_order<-phylogenetic_species_order(objects$species)
object_groups<-lapply(species_order,function(s){d<-filter(objects,species==s);split(d,d$specimen_id)});names(object_groups)<-species_order
boot_rows<-vector('list',B);boot_slopes<-vector('list',B)
for(b in seq_len(B)){
  draw<-bind_rows(lapply(species_order,function(s){
    group<-object_groups[[s]];selected<-sample.int(length(group),length(group),replace=TRUE)
    values<-vapply(selected,function(j){
      animal<-group[[j]];picked<-sample.int(nrow(animal),nrow(animal),replace=TRUE);sample<-animal[picked,]
      c(cell=median(sample$cell_area_um2),nucleus=median(sample$nuc_area_um2),fraction=median(sample$nuclear_fraction),cytoplasm=median(sample$paired_cytoplasm_um2))
    },numeric(4));values<-rowMeans(matrix(values,nrow=4))
    tibble(replicate=b,species=s,cell_area_um2=values[1],nucleus_area_um2=values[2],nuclear_fraction=values[3],cytoplasm_area_um2=values[4])
  }))
  stopifnot(nrow(draw)==24L,setequal(draw$species,species_order),all(draw$nucleus_area_um2<draw$cell_area_um2))
  boot_rows[[b]]<-draw
  boot_slopes[[b]]<-tryCatch(fit_allometry(draw,paste0('bootstrap_',b))$summary %>% mutate(replicate=b,fit_status='success',failure=NA_character_),
    error=function(e)tibble(analysis_id=paste0('bootstrap_',b),replicate=b,fit_status='failure',failure=conditionMessage(e)))
  if(b%%250L==0L)message('Figure10: ',b,'/',B,' paired hierarchical refits')
}
boot_traits<-bind_rows(boot_rows);boot_models<-bind_rows(boot_slopes)
stopifnot(nrow(boot_models)==B,all(boot_models$fit_status=='success'))
boot_summary<-boot_models %>% summarise(n_replicates=n(),median_slope=median(slope),ci_low=unname(quantile(slope,.025)),ci_high=unname(quantile(slope,.975)),fraction_slope_below_one=mean(slope<1),n_lambda_lower_boundary=sum(lambda<=1.1e-7),interval_method='95% percentile range across selected animal and paired-cell resampling refits',measurement_error_corrected=FALSE)
allocation<-balanced %>% left_join(boot_traits %>% group_by(species) %>% summarise(fraction_low=unname(quantile(nuclear_fraction,.025)),fraction_high=unname(quantile(nuclear_fraction,.975)),.groups='drop'),by='species',relationship='one-to-one') %>%
  mutate(y=match(species,rev(species_order)),support=paste0(n_animals,' / ',n_paired_cells))
animals<-animals %>% left_join(select(allocation,species,y),by='species',relationship='many-to-one') %>% group_by(species) %>% mutate(animal_y=y+(row_number()-(n()+1)/2)*.15) %>% ungroup()
x<-seq(min(primary$cell_area_um2),max(primary$cell_area_um2),length.out=160);X<-cbind(1,log10(x));mean_log<-as.numeric(X%*%coef(fits$released24$fit));se_log<-sqrt(pmax(rowSums((X%*%vcov(fits$released24$fit))*X),0));critical<-qt(.975,nrow(primary)-2L)
isometric_intercept<-mean(log10(primary$nucleus_area_um2))-mean(log10(primary$cell_area_um2))
predictions<-tibble(cell_area_um2=x,fitted=10^mean_log,low=10^(mean_log-critical*se_log),high=10^(mean_log+critical*se_log),isometry=10^(isometric_intercept+log10(x)))
write_csv(primary,file.path(ANALYSIS,'released_species_data.csv'));write_csv(animals,file.path(ANALYSIS,'paired_animal_summaries.csv'));write_csv(allocation,file.path(ANALYSIS,'equal_animal_species_data.csv'))
write_csv(select(objects,species,specimen_id,cell_object_id,nucleus_object_id,cell_area_um2,nuc_area_um2,nuclear_fraction,paired_cytoplasm_um2),file.path(ANALYSIS,'paired_allocation_data.csv'))
write_csv(models,file.path(ANALYSIS,'allometry_model_summary.csv'));write_csv(predictions,file.path(ANALYSIS,'fitted_relationship.csv'))
write_csv(boot_traits,file.path(ANALYSIS,'paired_hierarchical_traits.csv.gz'));write_csv(boot_models,file.path(ANALYSIS,'paired_hierarchical_refits.csv'));write_csv(boot_summary,file.path(ANALYSIS,'paired_hierarchical_summary.csv'))
write.tree(prune_source_tree(species_order),file.path(ANALYSIS,'allometry24.tre'));write.tree(prune_source_tree(setdiff(species_order,c('aeneus','wrighti'))),file.path(ANALYSIS,'allometry22.tre'))
colors<-c('Aquatic larvae'='#666666','Direct development'=GBE_COLORS[['purple']]);shapes<-c('Aquatic larvae'=21,'Direct development'=24)
p_format<-function(p)ifelse(p<.001,formatC(p,format='e',digits=1),sprintf('%.3f',p))
allometry_plot<-function(label_all=FALSE){
  labels<-primary %>% mutate(label=if(label_all)morphology_label else if_else(species %in% c('aeneus','wrighti','campi','ochrophaeus','welteri','brimleyorum'),morphology_label,''))
  m<-models[1,];annotation<-sprintf('Slope = %.2f [%.2f, %.2f]\nP(slope = 1) = %s; n = 24',m$slope,m$ci_low,m$ci_high,p_format(m$p_vs_isometry))
  ggplot(primary,aes(cell_area_um2,nucleus_area_um2))+
    geom_ribbon(data=predictions,aes(y=NULL,ymin=low,ymax=high),fill='#E6E6E6',alpha=.85)+
    geom_line(data=predictions,aes(y=fitted,linetype='PGLS'),linewidth=.6)+
    geom_line(data=predictions,aes(y=isometry,linetype='Proportional scaling (slope = 1)'),colour='#888888',linewidth=.5)+
    geom_errorbar(aes(xmin=cell_area_um2_ci_low,xmax=cell_area_um2_ci_high),orientation='y',width=0,colour='#BBBBBB',linewidth=.25)+
    geom_errorbar(aes(ymin=nucleus_area_um2_ci_low,ymax=nucleus_area_um2_ci_high),width=0,colour='#BBBBBB',linewidth=.25)+
    geom_point(aes(colour=development,shape=development),size=2.3,stroke=.65,fill='white')+
    geom_label_repel(data=labels,aes(label=label),fill="white",linetype=0,label.size=NA,label.padding=grid::unit(.08,"lines"),label.r=grid::unit(0,"lines"),family=GBE_FONT,fontface='italic',size=7/ggplot2::.pt,seed=20260924,box.padding=.45,point.padding=.25,max.overlaps=Inf,max.time=Inf,max.iter=150000,force=8,force_pull=.02,min.segment.length=0,segment.color='#AAAAAA',segment.size=.2)+
    annotate('text',x=min(primary$cell_area_um2_ci_low),y=Inf,hjust=0,vjust=1.1,label=annotation,family=GBE_FONT,size=7.5/ggplot2::.pt,lineheight=1.1)+
    scale_x_log10(breaks=c(200,250,300,350),expand=expansion(mult=c(.08,.16)))+
    scale_y_log10(breaks=c(30,40,50,60),expand=expansion(mult=c(.1,.27)))+
    scale_colour_manual(values=colors,name=NULL)+scale_shape_manual(values=shapes,name=NULL)+
    scale_linetype_manual(name=NULL,values=c('PGLS'='solid','Proportional scaling (slope = 1)'='dashed'))+
    labs(x='Cell area (µm²; log scale)',y='Nucleus area (µm²; log scale)')+theme_gbe()+
    guides(colour=guide_legend(order=1),shape=guide_legend(order=1),linetype=guide_legend(order=2))+
    theme(legend.box='vertical',plot.margin=margin(5,7,5,7))
}
allocation_plot<-ggplot(allocation,aes(nuclear_fraction,y))+
  geom_errorbar(aes(xmin=fraction_low,xmax=fraction_high,colour=development),orientation='y',width=.16,linewidth=.6)+
  geom_point(data=animals,aes(nuclear_fraction,animal_y,colour=development),shape=1,size=1.6,stroke=.55)+
  geom_point(aes(colour=development),shape=18,size=2.3)+
  geom_text(aes(x=.237,label=support),colour='#555555',family=GBE_FONT,size=7/ggplot2::.pt,hjust=.5)+
  annotate('text',x=.237,y=25.25,label='Animals / cells',family=GBE_FONT,size=7/ggplot2::.pt)+
  scale_x_continuous(labels=label_percent(accuracy=1),breaks=seq(.12,.22,.02),limits=c(.115,.256),expand=expansion(mult=0))+
  scale_y_continuous(breaks=allocation$y,labels=allocation$morphology_label,limits=c(.35,25.8),expand=expansion(mult=0))+
  scale_colour_manual(values=colors,name=NULL)+labs(x='Nuclear fraction of cell area',y=NULL,
    caption='Open circles: animal medians. Diamonds: equal-animal species means.\nBars: 95% paired animal/cell resampling intervals; single-animal intervals omit between-animal uncertainty.')+
  theme_gbe()+theme(legend.position='none',axis.text.y=element_text(face='italic'),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(5,7,5,7))
forest<-models %>% mutate(label=c('Released medians; 24 species','Equal-animal means; 24 species','Released medians; 22 species','Equal-animal means; 22 species'),interval_type='PGLS 95% t interval') %>% select(analysis_id,label,slope,ci_low,ci_high,interval_type)
forest<-bind_rows(forest,tibble(analysis_id='paired_bootstrap24',label='Paired resampling; 24 species',slope=boot_summary$median_slope,ci_low=boot_summary$ci_low,ci_high=boot_summary$ci_high,interval_type='95% resampling range')) %>% mutate(y=5:1)
write_csv(forest,file.path(ANALYSIS,'displayed_slope_intervals.csv'))
slope_plot<-ggplot(forest,aes(slope,y))+geom_vline(xintercept=1,linetype='dashed',colour='#888888',linewidth=.4)+geom_vline(xintercept=0,colour='#CCCCCC',linewidth=.3)+
  geom_errorbar(aes(xmin=ci_low,xmax=ci_high,colour=interval_type),orientation='y',width=.12,linewidth=.6)+
  geom_point(aes(shape=interval_type,colour=interval_type),size=2.2,fill='white',stroke=.6)+
  scale_shape_manual(values=c('PGLS 95% t interval'=21,'95% resampling range'=18),name=NULL)+
  scale_colour_manual(values=c('PGLS 95% t interval'='#222222','95% resampling range'=GBE_COLORS[['blue']]),name=NULL)+
  scale_y_continuous(breaks=forest$y,labels=forest$label,limits=c(.5,5.6))+scale_x_continuous(breaks=seq(0,1.25,.25))+
  labs(x='Unstandardized log–log slope (1 = proportional scaling)',y=NULL,caption='The 22-species sensitivity excludes D. aeneus and D. wrighti.\nResampling ranges quantify selected-measurement sensitivity, not total model uncertainty.')+
  theme_gbe()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(5,7,5,7))
legend<-paste('Nuclear–cell allometry and allocation among selected large erythrocytes, using all 1,152 reviewed paired objects from 42 animals and 24 species.',
 '(a) Primary Pagel-lambda PGLS regresses log10 nucleus area on log10 cell area using released species medians, exactly reproducing the Figure 5 input summaries. The slope is unstandardized; slope one means proportional area scaling. The dashed slope-one reference passes through the log-scale species centroid. The solid line and gray ribbon show the fitted conditional mean and 95% t coefficient band; inference is conditional on fitted lambda.',
 'Point bars are released paired-object measurement intervals on both axes. The conditional regression does not correct predictor measurement error or establish a structural allometric exponent.',
 '(b) Open circles are animal medians of actual paired-object nucleus/cell fractions. Diamonds are species means of those animal medians, with equal animal weights. Species follow shared phylogenetic order; right-hand counts report animals / paired cells. Bars use 1,000 paired hierarchical resamples: sample animals within species, then paired cells within each sampled animal. Purple denotes the two direct developers; other species have aquatic larvae.',
 '(c) Conditional PGLS slopes and 95% t intervals compare released pooled medians with equal-animal means of animal medians, with and without D. aeneus and D. wrighti. The blue diamond and interval are the median and central95% range of 1,000 equal-animal bootstrap PGLS refits, re-estimating lambda; this is selected-measurement sensitivity, not a combined model confidence interval.',
 'Every model has species, not cells, as observational units. Seven species have one animal; their resampling intervals cannot represent among-animal uncertainty. One D. orestes animal contributes only one accepted cell, so within-animal dispersion is not estimable.',
 'Cells were selected for large area rather than randomly sampled from circulating erythrocytes; results describe that selected panel. Nuclear fraction and cytoplasmic area are derived paired quantities, not independent variables for causal path analysis. No cell, animal or species is removed from the primary analysis. Source-tree calibration/citation provenance remains unresolved. These exploratory complementary displays use unadjusted tests and do not establish physiological function or causal direction.')
panels<-list(a=list(plot=allometry_plot(TRUE),title='Nuclear–cell allometry',width_mm=180,height_mm=125,legend=legend,alt_text='Twenty-four labelled species show nucleus versus cell area on logarithmic axes, a conditional PGLS line and proportional-scaling reference.'),
 b=list(plot=allocation_plot,title='Paired nuclear allocation by animal and species',width_mm=180,height_mm=155,legend=legend,alt_text='Forty-two animal medians and twenty-four equal-animal species summaries show nuclear-area fractions in phylogenetic order, with support counts.'),
 c=list(plot=slope_plot,title='Allometry sensitivity to animal balance and species panel',width_mm=180,height_mm=90,legend=legend,alt_text='Five intervals compare allometric slopes from four conditional models and paired resampling against proportional scaling at slope one.'))
collage<-(allometry_plot(FALSE)/allocation_plot/slope_plot)+plot_layout(heights=c(1,1.38,.7))+plot_annotation(tag_levels='a')
save_candidate(10,'Nuclear–cell scaling and paired allocation',collage,panels,180,315,'S18;S19;S26;species_taxonomy_crosswalk;trees/desmo900dated_source46.tre',legend,'Three panels connect species nuclear–cell allometry, animal nuclear allocation and sensitivity to unequal animal representation.')
write_lines(c('# Figure 10: nuclear–cell scaling and paired allocation','',
 'Exact frozen inputs are copied with SHA-256 hashes to inputs/nuclear_allometry. Identity joins use the species crosswalk. No prior figure or frozen release is altered.','',
 'Primary areas are released pooled species medians. Equal-animal alternative: arithmetic mean of within-animal medians, separately for each trait. Nuclear fraction is computed for each actual pair, then summarized by animal median and equal-animal species mean; it is not a ratio of aggregate medians.','',
 'phylolm fits log10(nucleus area) ~ log10(cell area) by ML with Pagel lambda in [1e-7,1]. Slopes are unstandardized. Isometry test: t=(slope-1)/SE with n-2 degrees of freedom. Model intervals use t critical values conditional on fitted lambda. Four fixed sensitivities use released/equal-animal summaries and all24/excluding aeneus+wrighti (22).','',
 'Seed20260924; 1,000 draws resample animals within species, then paired cells within each sampled animal. Equal-animal means of sampled animal medians are fitted across all24 species, re-estimating lambda. All failures and warnings are retained; every fit must succeed before export. This is selected-measurement sensitivity, not a full errors-in-variables model or combined model confidence interval. Seven single-animal species cannot identify between-animal uncertainty.','',
 'Conditional response-on-predictor regression is directional; predictor error and axis choice can affect the exponent. Do not interpret the slope as a demonstrated mechanistic scaling law. Nuclear fraction is descriptive and is not tested against its own component variables as an independent outcome.','',
 'Outputs: allometry_model_summary.csv has four models and slope=0/slope=1 tests. paired_animal_summaries.csv and equal_animal_species_data.csv contain plotted data. paired_hierarchical_refits.csv contains all fits. paired_hierarchical_traits.csv.gz stores all24,000 draws. displayed_slope_intervals.csv distinguishes t confidence intervals from selected-measurement ranges.','',
 sprintf('Primary conditional slope %.4f [%.4f,%.4f], P(slope=1) %.6g; lambda %.7f.',models$slope[1],models$ci_low[1],models$ci_high[1],models$p_vs_isometry[1],models$lambda[1]),
 sprintf('Equal-animal full24 slope %.4f; bootstrap median %.4f, central95%% [%.4f,%.4f].',models$slope[2],boot_summary$median_slope,boot_summary$ci_low,boot_summary$ci_high)),file.path(ANALYSIS,'README.md'))
write_lines(c('PASS: unique nucleus-cell pairs; positive nested areas; fraction/cytoplasm arithmetic.',
 'PASS: 1,152 objects / 42 animals / 24 species; seven one-animal species; frozen species medians reproduced.',
 'PASS: Figure5 standardized association and slope-zero P reproduced after reversing conditional axes.',
 'PASS: unstandardized slopes, n-2 df, t coefficient intervals and slope-one tests.',
 'PASS: 1,000 paired hierarchical draws preserve24 species; 1,000 PGLS fits succeeded.',
 'PASS: frozen snapshot/source SHA-256 agreement.'),file.path(ANALYSIS,'VERIFICATION.txt'))
message('Figure10 complete'); print(models %>% select(analysis_id,n_species,slope,ci_low,ci_high,p_vs_isometry,lambda));print(boot_summary)
