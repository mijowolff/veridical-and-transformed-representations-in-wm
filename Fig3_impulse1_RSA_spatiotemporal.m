%% Concurrent maintenance of veridical and transformend representations in WM

% Figure 3, top row, sptatio-temporal 
% Impulse 1 cross-validated RSA

% written by Michael J. Wolff,2022

clear all;
close all;
clc;
%%
main_dir=pwd; % path to folder
addpath(genpath(main_dir))

test_chans = { 'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O2','O1','Oz'}';

nreps=100;% number of repeats for random subsampling and cross-validation
nreps_null=1000; % number of repeats for label randomisation for null-distribution
nfolds=8; % number of folds for cross-validation

span=5; % number of time-points to average over (5 = 10 ms for 500 Hz)
toi=[.1 .4001]; % time-window of interest

% load in previous results
fname_in='imp1_CV_RDM';
fname_in_null='imp1_CV_RDM_NULL';

fname_out='NEW_imp1_CV_RDM';
fname_out_null='NEW_imp1_CV_RDM_NULL';

% run new analysis or load in previous results
do_decoding      = false; 
do_null_decoding = false;
%%
if do_decoding||do_null_decoding
    for sub=1:30
        fprintf(['Subject ' num2str(sub) '\n']) % indicator
                
        load(['FT_Gen_2_WM_Rotation_imp1_' num2str(sub) '.mat']);
        
        Results=ft_imp1.Results;
        Results_header=ft_imp1.Results_header';
        
        incl=setdiff(1:size(Results,1),[ft_imp1.bad_trials]);
        
        % extract good trials, channels and time window of interest of
        % impulse 1 epoch and reformat
        dat_temp=ft_imp1.trial(incl,ismember(ft_imp1.label,test_chans),ft_imp1.time>toi(1)&ft_imp1.time<=toi(2));
        dat_temp=bsxfun(@minus,dat_temp,mean(dat_temp,3)); % take relative baseline
        dat_temp=movmean(dat_temp,span,3,'Endpoints','discard'); % downsample
        dat_temp=dat_temp(:,:,1:span:end);
        dat_imp1=reshape(dat_temp,[size(dat_temp,1),size(dat_temp,2)*size(dat_temp,3)]); %combine channel and time dimensions
        
        cue_cond=Results(incl,18); % cued location
        rot_cond=Results(incl,22); % rotation condition
        cued_rad=circ_dist(Results(incl,16),0); % cued item in radians
        uncued_rad=circ_dist(Results(incl,17),0); % uncued item in radians
        rot_rad=circ_dist(Results(incl,28),0); % rotated item in radians
        %%
        conditions=[cue_cond,cued_rad,uncued_rad]; % these are all conditions of interest
        
        if do_decoding % do cross-validated RSA
            
            cfg.nfolds  = nfolds;
            cfg.nreps   = nreps;
            cfg.avg     = true;
            cfg.null    = false;
            cfg.verbose = true;
            
            [RDM,cond_combs] = mahal_CV_RSA(dat_imp1,conditions,cfg);
            
            imp1_RDMs(sub,:,:)=RDM;
        end
        if do_null_decoding % same as above, but shuffles trial labels
            
            cfg_null.nfolds  = nfolds;
            cfg_null.nreps   = nreps_null;
            cfg_null.avg     = false;
            cfg_null.null    = true;
            cfg_null.verbose = true;
            
            RDM_null = mahal_CV_RSA(dat_imp1,conditions,cfg_null);
            
            imp1_RDMs_null(sub,:,:,:)=RDM_null;
        end
    end
    if exist('imp1_RDMs','var')
        save(fullfile([main_dir '\results'],[fname_out '_output.mat']),'imp1_RDMs','cond_combs','cfg')
    end
    if exist('imp1_RDMs_null','var')
        save(fullfile([main_dir '\results'],[fname_out_null '_output.mat']),'imp1_RDMs_null','cfg_null')
    end
end
if ~do_decoding
    load(fullfile([main_dir '\results'],[fname_in '_output.mat']))
end
if ~do_null_decoding
    load(fullfile([main_dir '\results'],[fname_in_null '_output.mat']))
end
%% make models

% cued location
cued_loc_model=abs(cond_combs(:,1)-cond_combs(:,1)');

% cued item
cued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,2),cond_combs(:,2)))));
% same cue
cued_same_model=cued_model;
cued_same_model(cued_loc_model==1)=nanmean(cued_model(cued_loc_model==0),1);
% different cue
cued_diff_model=cued_model;
cued_diff_model(cued_loc_model==0)=nanmean(cued_model(cued_loc_model==1),1);

% uncued item
uncued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,3),cond_combs(:,3)))));
% same cue
uncued_same_model=uncued_model;
uncued_same_model(cued_loc_model==1)=nanmean(uncued_model(cued_loc_model==0),1);
% different cue
uncued_diff_model=uncued_model;
uncued_diff_model(cued_loc_model==0)=nanmean(uncued_model(cued_loc_model==1),1);

% cued/uncued generalization
cued_uncued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,3),cond_combs(:,2)))));
uncued_cued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,2),cond_combs(:,3)))));
% same cue
cued_uncued_same_model=cued_uncued_model;
cued_uncued_same_model(cued_loc_model==1)=nanmean(cued_uncued_same_model(cued_loc_model==0),1);
uncued_cued_same_model=uncued_cued_model;
uncued_cued_same_model(cued_loc_model==1)=nanmean(uncued_cued_same_model(cued_loc_model==0),1);
% different cue
cued_uncued_diff_model=cued_uncued_model;
cued_uncued_diff_model(cued_loc_model==0)=nanmean(cued_uncued_diff_model(cued_loc_model==1),1);
uncued_cued_diff_model=uncued_cued_model;
uncued_cued_diff_model(cued_loc_model==0)=nanmean(uncued_cued_diff_model(cued_loc_model==1),1);
%% run GLM on RDAs, none of the models of interest are statistically related, so can all be run at once.

models=[];
models(:,:,1)=cued_loc_model;
models(:,:,2)=cued_same_model;
models(:,:,3)=cued_diff_model;
models(:,:,4)=uncued_same_model;
models(:,:,5)=uncued_diff_model;
models(:,:,6)=cued_uncued_same_model+uncued_cued_same_model;
models(:,:,7)=cued_uncued_diff_model+uncued_cued_diff_model;
%
RDM_res1=[];
RDM=imp1_RDMs;
RDM_null=imp1_RDMs_null;
betas=nan(size(RDM,1),size(models,3)+1);
betas_NULL=nan(size(RDM_null,1),size(RDM_null,2),size(models,3)+1);
X=[];
for m=1:size(models,3)
    model_temp=models(:,:,m);
    X(:,m)=zscore(model_temp(1:end));
end
X(:,end+1)=1;
for sub=1:size(RDM,1)
    temp_dists=squeeze(RDM(sub,:,:));
    temp_dists_null=squeeze(RDM_null(sub,:,:,:));
    Y = zscore(temp_dists(1:end));
    betas(sub,:)=pinv(X)*(Y)';
    for irep=1:size(temp_dists_null,1)
        Y = zscore(temp_dists_null(irep,1:end));
        betas_NULL(sub,irep,:)=pinv(X)*(Y)';
    end
end
beta_side=betas(:,1);
beta_cued_s=betas(:,2);
beta_cued_d=betas(:,3);
beta_uncued_s=betas(:,4);
beta_uncued_d=betas(:,5);
beta_cued_uncued_gen_s=betas(:,6);
beta_cued_uncued_gen_d=betas(:,7);

beta_side_NULL=betas_NULL(:,:,1);
beta_cued_s_NULL=betas_NULL(:,:,2);
beta_cued_d_NULL=betas_NULL(:,:,3);
beta_uncued_s_NULL=betas_NULL(:,:,4);
beta_uncued_d_NULL=betas_NULL(:,:,5);
beta_cued_uncued_gen_s_NULL=betas_NULL(:,:,6);
beta_cued_uncued_gen_d_NULL=betas_NULL(:,:,7);
%% compute t values for statistical significance tests
beta_side_T=FastTtest(beta_side);
beta_cued_s_T=FastTtest(beta_cued_s);
beta_cued_d_T=FastTtest(beta_cued_d);
beta_uncued_s_T=FastTtest(beta_uncued_s);
beta_uncued_d_T=FastTtest(beta_uncued_d);
beta_cued_uncued_gen_s_T=FastTtest(beta_cued_uncued_gen_s);
beta_cued_uncued_gen_d_T=FastTtest(beta_cued_uncued_gen_d);

beta_side_NT=FastTtest(beta_side_NULL);
beta_cued_s_NT=FastTtest(beta_cued_s_NULL);
beta_cued_d_NT=FastTtest(beta_cued_d_NULL);
beta_uncued_s_NT=FastTtest(beta_uncued_s_NULL);
beta_uncued_d_NT=FastTtest(beta_uncued_d_NULL);
beta_cued_uncued_gen_s_NT=FastTtest(beta_cued_uncued_gen_s_NULL);
beta_cued_uncued_gen_d_NT=FastTtest(beta_cued_uncued_gen_d_NULL);
%%
p_side=FastPvalue(beta_side_T,beta_side_NT,2); % two-sided

p_c_s=FastPvalue(beta_cued_s_T,beta_cued_s_NT,2);
p_c_d=FastPvalue(beta_cued_d_T,beta_cued_d_NT,2);
p_c_svd=FastPvalue(beta_cued_s_T-beta_cued_d_T,beta_cued_s_NT-beta_cued_d_NT,2);

p_u_s=FastPvalue(beta_uncued_s_T,beta_uncued_s_NT,2);
p_u_d=FastPvalue(beta_uncued_d_T,beta_uncued_d_NT,2);
p_u_svd=FastPvalue(beta_uncued_s_T-beta_uncued_d_T,beta_uncued_s_NT-beta_uncued_d_NT,2);

p_g_s=FastPvalue(beta_cued_uncued_gen_s_T,beta_cued_uncued_gen_s_NT,2);
p_g_d=FastPvalue(beta_cued_uncued_gen_d_T,beta_cued_uncued_gen_d_NT,2);
p_g_svd=FastPvalue(beta_cued_uncued_gen_s_T-beta_cued_uncued_gen_d_T,beta_cued_uncued_gen_s_NT-beta_cued_uncued_gen_d_NT,2);
%% CIs for plotting
beta_side_ci=bootci(10000,@mean,beta_side);
beta_cued_s_ci=bootci(10000,@mean,beta_cued_s);
beta_cued_d_ci=bootci(10000,@mean,beta_cued_d);
beta_uncued_s_ci=bootci(10000,@mean,beta_uncued_s);
beta_uncued_d_ci=bootci(10000,@mean,beta_uncued_d);
beta_cued_uncued_gen_s_ci=bootci(10000,@mean,beta_cued_uncued_gen_s);
beta_cued_uncued_gen_d_ci=bootci(10000,@mean,beta_cued_uncued_gen_d);
%%
clc
close all
figure('Renderer', 'painters', 'Position', [10 10 740 740])
subplot(2,2,1)
imagesc(squeeze(nanmean(imp1_RDMs,1)),[-0 2]); axis xy
set(gca,'Xtick',1:72)
set(gca,'Ytick',1:72)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
pbaspect([1 1 1])
colormap(flipud(autumn))
c=colorbar;
set(gca,'TickDir','out')
c.Label.String = 'squared mahalanobis distance';
c.Label.FontSize = 13;
title('Impulse 1, RDM')

%
subplot(2,2,2)
pos=[1];
hold all
plot([0 2],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([beta_side],...
    'positions',pos,'Widths',0.2,'Symbol','','Labels',{'cued side'});
set(gca,'FontSize',10,'XTickLabelRotation',-40)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(b1(:,1),'color',[.5 .5 .5]);
set(b1,'LineWidth', 2);
plot(pos(1),mean(beta_side,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
if p_side<0.05
    plot([pos(1) pos(1)],max(beta_side)+max(beta_side)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
plot([pos(1) pos(1)],beta_side_ci','Color',[0 0 0],'LineWidth',4)
s=swarmchart(pos(1).*ones(1,length(beta_side)),beta_side,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5]);
s.XJitter = 'density';
s.XJitterWidth = .2;
set(gca,'TickDir','out')
ylim([-.02 .35])
xlim([.5 1.5])
pbaspect([1 2 1])
ylabel('beta')

%
subplot(2,2,3:4)
pos=[1 2 3.5 4.5 6 7];
hold all
plot([0 8],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([beta_cued_s,beta_cued_d,beta_uncued_s,beta_uncued_d,beta_cued_uncued_gen_s,beta_cued_uncued_gen_d],...
    'positions',pos,'Widths',0.3,'Symbol','','Labels',{'cued, within','cued, across',...
    'uncued, within','uncued, across','gen., within','gen., across'});

set(gca,'FontSize',10,'XTickLabelRotation',-40)
set(b1,'LineWidth', 1.5);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(b1(:,1),'color',[0 0 1]);set(b1(:,2),'color',[0 0 1]);set(b1(:,3),'color',[1 0 0]);
set(b1(:,4),'color',[1 0 0]);set(b1(:,5),'color',[.5 0 .5]);set(b1(:,6),'color',[.5 0 .5]);
plot(pos(1),mean(beta_cued_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',9)
plot(pos(2),mean(beta_cued_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(3),mean(beta_uncued_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(4),mean(beta_uncued_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(5),mean(beta_cued_uncued_gen_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(6),mean(beta_cued_uncued_gen_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)

plot([pos(1) pos(1)],beta_cued_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(2) pos(2)],beta_cued_d_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(3) pos(3)],beta_uncued_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(4) pos(4)],beta_uncued_d_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(5) pos(5)],beta_cued_uncued_gen_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(6) pos(6)],beta_cued_uncued_gen_d_ci','Color',[0 0 0],'LineWidth',4)

s=swarmchart(pos(1)*ones(1,length(beta_side)),beta_cued_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(2)*ones(1,length(beta_side)),beta_cued_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(3)*ones(1,length(beta_side)),beta_uncued_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(4)*ones(1,length(beta_side)),beta_uncued_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(5)*ones(1,length(beta_side)),beta_cued_uncued_gen_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[.5 0 .5],'MarkerEdgeColor',[.5 0 .5]);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(6)*ones(1,length(beta_side)),beta_cued_uncued_gen_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[.5 0 .5],'MarkerEdgeColor',[.5 0 .5]);
s.XJitter = 'density';
s.XJitterWidth = .2;

if p_c_s<0.05
    plot([pos(1) pos(1)],max(beta_cued_s)+max(beta_cued_s)/5,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_c_d<0.05
    plot([pos(1) pos(1)],max(beta_cued_d)+max(beta_cued_d)/5,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_c_svd<0.05
    plot([1.5 1.5],max(beta_cued_s)+max(beta_cued_s)/2,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
    plot([1 2],[max(beta_cued_s)+max(beta_cued_s)/3 max(beta_cued_s)+max(beta_cued_s)/3 ],'Color','k','LineWidth',1,'LineStyle','-')
end
ylim([-.04 .08])
xlim([.5 7.5])
set(gca,'TickDir','out')
ylabel('beta')

axes('pos',[.70 .93 .07 .07])
imagesc([.1 .1], [.5 .5],cued_loc_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]); 
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.15 .457 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(cued_same_model,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.26 .457 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(cued_diff_model,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.427 .457 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(uncued_same_model,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.537 .457 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(uncued_diff_model,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.705 .457 .07 .07])
imagesc([.1 .1], [.5 .5],uncued_cued_same_model+cued_uncued_same_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.815 .457 .07 .07])
imagesc([.1 .1], [.5 .5],uncued_cued_diff_model+cued_uncued_diff_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

annotation('textbox',[.05 1 .0 .0],'String','A','EdgeColor','none','FontSize',24)
annotation('textbox',[.55 1 .0 .0],'String','B','EdgeColor','none','FontSize',24)
annotation('textbox',[.05 .55 .0 .0],'String','C','EdgeColor','none','FontSize',24)
