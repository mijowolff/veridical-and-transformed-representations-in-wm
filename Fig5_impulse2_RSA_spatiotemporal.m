%% Concurrent maintenance of veridical and transformend representations in WM

% Figure 5, top row, sptatio-temporal 
% Impulse 2 cross-validated RSA

% written by Michael J. Wolff, 2022

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
fname_in='imp2_CV_RDM';
fname_in_null='imp2_CV_RDM_NULL';

fname_out='NEW_imp2_CV_RDM'; 
fname_out_null='NEW_imp2_CV_RDM_NULL';

% run new analysis or load in previous results
do_decoding      = false; 
do_null_decoding = false;
%%
if do_decoding||do_null_decoding
    for sub=1:30
        fprintf(['Subject ' num2str(sub) '\n']) % indicator
        
        load(['FT_Gen_2_WM_Rotation_imp2_' num2str(sub) '.mat']);
        
        Results=ft_imp2.Results;
        Results_header=ft_imp2.Results_header';
        
        incl=setdiff(1:size(Results,1),[ft_imp2.bad_trials]);
        
        % extract good trials, channels and time window of interest of
        % impulse 2 epoch and reformat
        dat_temp=ft_imp2.trial(incl,ismember(ft_imp2.label,test_chans),ft_imp2.time>toi(1)&ft_imp2.time<=toi(2));
        dat_temp=bsxfun(@minus,dat_temp,mean(dat_temp,3)); % take relative baseline
        dat_temp=movmean(dat_temp,span,3,'Endpoints','discard'); % downsample
        dat_temp=dat_temp(:,:,1:span:end);
        dat_imp2=reshape(dat_temp,[size(dat_temp,1),size(dat_temp,2)*size(dat_temp,3)]); %combine channel and time dimensions
        
        cue_cond=Results(incl,18); % cued location
        rot_cond=Results(incl,22); % rotation condition
        cued_rad=circ_dist(Results(incl,16),0); % cued item in radians
        uncued_rad=circ_dist(Results(incl,17),0); % uncued item in radians
        rot_rad=circ_dist(Results(incl,28),0); % rotated item in radians
        %%
        conditions=[cue_cond,cued_rad,rot_rad]; % these are all conditions of interest
        
        if do_decoding % do cross-validated RSA
            
            cfg.nfolds  = nfolds;
            cfg.nreps   = nreps;
            cfg.avg     = true;
            cfg.null    = false;
            cfg.verbose = true;
            
            [RDM,cond_combs] = mahal_CV_RSA(dat_imp2,conditions,cfg);
            
            imp2_RDMs(sub,:,:)=RDM;
        end
        if do_null_decoding % same as above, but shuffles trial labels
            
            cfg_null.nfolds  = nfolds;
            cfg_null.nreps   = nreps_null;
            cfg_null.avg     = false;
            cfg_null.null    = true;
            cfg_null.verbose = true;
            
            RDM_null = mahal_CV_RSA(dat_imp2,conditions,cfg_null);
            
            imp2_RDMs_null(sub,:,:,:)=RDM_null;
        end
    end
    if exist('imp2_RDMs','var')
        save(fullfile([main_dir '\results'],[fname_out '_output.mat']),'imp2_RDMs','cond_combs','cfg')
    end
    if exist('imp2_RDMs_null','var')
        save(fullfile([main_dir '\results'],[fname_out_null '_output.mat']),'imp2_RDMs_null','cfg_null')
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
% same cued location
cued_same_model=cued_model;
cued_same_model(cued_loc_model==1)=nanmean(cued_model(cued_loc_model==0),1);
% different cued location
cued_diff_model=cued_model;
cued_diff_model(cued_loc_model==0)=nanmean(cued_model(cued_loc_model==1),1);

% rotation product
rp_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,3),cond_combs(:,3)))));
% same cued location
rp_same_model=rp_model;
rp_same_model(cued_loc_model==1)=nanmean(rp_model(cued_loc_model==0),1);
% different cued location
rp_diff_model=rp_model;
rp_diff_model(cued_loc_model==0)=nanmean(rp_model(cued_loc_model==1),1);

% cued/rp generalization
cued_rp_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,3),cond_combs(:,2)))));
rp_cued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,2),cond_combs(:,3)))));
% same cued location
cued_rp_same_model=cued_rp_model;
cued_rp_same_model(cued_loc_model==1)=nanmean(cued_rp_same_model(cued_loc_model==0),1);
rp_cued_same_model=rp_cued_model;
rp_cued_same_model(cued_loc_model==1)=nanmean(rp_cued_same_model(cued_loc_model==0),1);
% different cued location
cued_rp_diff_model=cued_rp_model;
cued_rp_diff_model(cued_loc_model==0)=nanmean(cued_rp_diff_model(cued_loc_model==1),1);
rp_cued_diff_model=rp_cued_model;
rp_cued_diff_model(cued_loc_model==0)=nanmean(rp_cued_diff_model(cued_loc_model==1),1);

% rotation condition models
rot_conds=round(circ_rad2ang(angle(exp(1i*cond_combs(:,3))./exp(1i*cond_combs(:,2)))));

u_rot_conds=unique(rot_conds);
for c=1:length(u_rot_conds)
    rot_cond_models(:,:,c)=-(abs((rot_conds==u_rot_conds(c))*(rot_conds==u_rot_conds(c))'));
    
    % same cued location
    temp=rot_cond_models(:,:,c);
    temp(cued_loc_model==0)=nanmean(temp(cued_loc_model==1),1);
    rot_cond_models_diff(:,:,c)=temp;
    
    % different cued location
    temp=rot_cond_models(:,:,c);
    temp(cued_loc_model==1)=nanmean(temp(cued_loc_model==0),1);
    rot_cond_models_same(:,:,c)=temp;
end
%% run GLM on RDAs, some models are statistically related, so they need to be regressed out first

% cued location and rotation instructions)
models=[];
models(:,:,1)=cued_loc_model;
models(:,:,2:6)=rot_cond_models_same;
models(:,:,7:11)=rot_cond_models_diff;
%
RDM=imp2_RDMs;
RDM_null=imp2_RDMs_null;
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
beta_rot_cond_s=mean(betas(:,2:6),2);
beta_rot_cond_d=mean(betas(:,7:11),2);
beta_side_NULL=betas_NULL(:,:,1);
beta_rot_cond_s_NULL=mean(betas_NULL(:,:,2:6),3);
beta_rot_cond_d_NULL=mean(betas_NULL(:,:,7:11),3);
%% cued item (regress out rotated)
RDM=imp2_RDMs;
RDM_null=imp2_RDMs_null;

models=[];
models(:,:,1)=cued_same_model;
models(:,:,2)=cued_diff_model;

betas=nan(size(RDM,1),size(models,3)+1);
betas_NULL=nan(size(RDM_null,1),size(RDM_null,2),size(models,3)+1);

X=[];
for m=1:size(models,3)
    model_temp=models(:,:,m);
    X(:,m)=zscore(model_temp(1:end));
end
X(:,end+1)=1;

X_res=[];
X_res(:,1)=zscore(rp_model(1:end));
X_res(:,end+1)=1;

for sub=1:size(RDM,1)
    temp_dists=squeeze(RDM(sub,:,:));
    temp_dists_null=squeeze(RDM_null(sub,:,:,:));
    Y = zscore(temp_dists(1:end));
    b=pinv(X_res)*(Y)';
    Y_res=zscore(Y'-X_res*b);
    betas(sub,:)=pinv(X)*(Y_res);
    for irep=1:size(temp_dists_null,1)
        Y = zscore(temp_dists_null(irep,1:end));
        b=pinv(X_res)*(Y)';
        Y_res=zscore(Y'-X_res*b);
        betas_NULL(sub,irep,:)=pinv(X)*(Y_res);
    end
end
beta_cued_s=betas(:,1);
beta_cued_d=betas(:,2);

beta_cued_s_NULL=betas_NULL(:,:,1);
beta_cued_d_NULL=betas_NULL(:,:,2);
%% rotated (regress out cued)
RDM=imp2_RDMs;
RDM_null=imp2_RDMs_null;

models=[];
models(:,:,1)=rp_same_model;
models(:,:,2)=rp_diff_model;

betas=nan(size(RDM,1),size(models,3)+1);
betas_NULL=nan(size(RDM_null,1),size(RDM_null,2),size(models,3)+1);

X=[];
for m=1:size(models,3)
    model_temp=models(:,:,m);
    X(:,m)=zscore(model_temp(1:end));
end
X(:,end+1)=1;

X_res=[];
X_res(:,1)=zscore(cued_model(1:end));
X_res(:,end+1)=1;

for sub=1:size(RDM,1)
    temp_dists=squeeze(RDM(sub,:,:));
    temp_dists_null=squeeze(RDM_null(sub,:,:,:));
    Y = zscore(temp_dists(1:end));
    b=pinv(X_res)*(Y)';
    Y_res=zscore(Y'-X_res*b);
    betas(sub,:)=pinv(X)*(Y_res);
    for irep=1:size(temp_dists_null,1)
        Y = zscore(temp_dists_null(irep,1:end));
        b=pinv(X_res)*(Y)';
        Y_res=zscore(Y'-X_res*b);
        betas_NULL(sub,irep,:)=pinv(X)*(Y_res);
    end
end
beta_rot_s=betas(:,1);
beta_rot_d=betas(:,2);

beta_rot_s_NULL=betas_NULL(:,:,1);
beta_rot_d_NULL=betas_NULL(:,:,2);
%% generalization (regress out cued and rotated)
RDM=imp2_RDMs;
RDM_null=imp2_RDMs_null;

models=[];
models(:,:,1)=cued_rp_same_model+rp_cued_same_model;
models(:,:,2)=cued_rp_diff_model+rp_cued_diff_model;

betas=nan(size(RDM,1),size(models,3)+1);
betas_NULL=nan(size(RDM_null,1),size(RDM_null,2),size(models,3)+1);

X=[];
for m=1:size(models,3)
    model_temp=models(:,:,m);
    X(:,m)=zscore(model_temp(1:end));
end
X(:,end+1)=1;

X_res=[];
X_res(:,1)=zscore(cued_model(1:end));
X_res(:,2)=zscore(rp_model(1:end));
X_res(:,end+1)=1;

for sub=1:size(RDM,1)
    temp_dists=squeeze(RDM(sub,:,:));
    temp_dists_null=squeeze(RDM_null(sub,:,:,:));
    Y = zscore(temp_dists(1:end));
    b=pinv(X_res)*(Y)';
    Y_res=zscore(Y'-X_res*b);
    betas(sub,:)=pinv(X)*(Y_res);
    for irep=1:size(temp_dists_null,1)
        Y = zscore(temp_dists_null(irep,1:end));
        b=pinv(X_res)*(Y)';
        Y_res=zscore(Y'-X_res*b);
        betas_NULL(sub,irep,:)=pinv(X)*(Y_res);
    end
end
beta_cued_rot_gen_s=betas(:,1);
beta_cued_rot_gen_d=betas(:,2);

beta_cued_rot_gen_s_NULL=betas_NULL(:,:,1);
beta_cued_rot_gen_d_NULL=betas_NULL(:,:,2);
%% compute t values for statistical significance tests

beta_side_T=FastTtest(beta_side);
beta_rot_cond_s_T=FastTtest(beta_rot_cond_s);
beta_rot_cond_d_T=FastTtest(beta_rot_cond_d);
beta_cued_s_T=FastTtest(beta_cued_s);
beta_cued_d_T=FastTtest(beta_cued_d);
beta_rot_s_T=FastTtest(beta_rot_s);
beta_rot_d_T=FastTtest(beta_rot_d);
beta_cued_rot_gen_s_T=FastTtest(beta_cued_rot_gen_s);
beta_cued_rot_gen_d_T=FastTtest(beta_cued_rot_gen_d);

beta_side_NT=FastTtest(beta_side_NULL);
beta_rot_cond_s_NT=FastTtest(beta_rot_cond_s_NULL);
beta_rot_cond_d_NT=FastTtest(beta_rot_cond_d_NULL);
beta_cued_s_NT=FastTtest(beta_cued_s_NULL);
beta_cued_d_NT=FastTtest(beta_cued_d_NULL);
beta_rot_s_NT=FastTtest(beta_rot_s_NULL);
beta_rot_d_NT=FastTtest(beta_rot_d_NULL);
beta_cued_rot_gen_s_NT=FastTtest(beta_cued_rot_gen_s_NULL);
beta_cued_rot_gen_d_NT=FastTtest(beta_cued_rot_gen_d_NULL);
%%
p_side=FastPvalue(beta_side_T,beta_side_NT,2); % two-sided

p_rco_s=FastPvalue(beta_rot_cond_s_T,beta_rot_cond_s_NT,2);
p_rco_d=FastPvalue(beta_rot_cond_d_T,beta_rot_cond_d_NT,2);
p_rco_svd=FastPvalue(beta_rot_cond_s_T-beta_rot_cond_d_T,beta_rot_cond_s_NT-beta_rot_cond_d_NT,2);

p_c_s=FastPvalue(beta_cued_s_T,beta_cued_s_NT,2);
p_c_d=FastPvalue(beta_cued_d_T,beta_cued_d_NT,2);
p_c_svd=FastPvalue(beta_cued_s_T-beta_cued_d_T,beta_cued_s_NT-beta_cued_d_NT,2);

p_r_s=FastPvalue(beta_rot_s_T,beta_rot_s_NT,2);
p_r_d=FastPvalue(beta_rot_d_T,beta_rot_d_NT,2);
p_r_svd=FastPvalue(beta_rot_s_T-beta_rot_d_T,beta_rot_s_NT-beta_rot_d_NT,2);

p_g_s=FastPvalue(beta_cued_rot_gen_s_T,beta_cued_rot_gen_s_NT,2);
p_g_d=FastPvalue(beta_cued_rot_gen_d_T,beta_cued_rot_gen_d_NT,2);
p_g_svd=FastPvalue(beta_cued_rot_gen_d_T-beta_cued_rot_gen_s_T,beta_cued_rot_gen_d_NT-beta_cued_rot_gen_s_NT,2);
%% CIs for plotting
beta_side_ci=bootci(10000,@mean,beta_side);
beta_rot_cond_s_ci=bootci(10000,@mean,beta_rot_cond_s);
beta_rot_cond_d_ci=bootci(10000,@mean,beta_rot_cond_d);
beta_cued_s_ci=bootci(10000,@mean,beta_cued_s);
beta_cued_d_ci=bootci(10000,@mean,beta_cued_d);
beta_rot_s_ci=bootci(10000,@mean,beta_rot_s);
beta_rot_d_ci=bootci(10000,@mean,beta_rot_d);
beta_cued_rot_gen_s_ci=bootci(10000,@mean,beta_cued_rot_gen_s);
beta_cued_rot_gen_d_ci=bootci(10000,@mean,beta_cued_rot_gen_d);
%%
side_color=[160/255 80/255 0];
rot_cond_s_color=[48/255 48/255 48/255];
rot_cond_d_color=[128/255 128/255 128/255];

cued_s_color=[0 0 128/255];
cued_d_color=[96/255 96/255 255/255];

rot_s_color=[0 128/255 0];
rot_d_color=[0 196/255 0];

gen_s_color=[0 128/255 128/255];
gen_d_color=[0 196/255 196/255];

close all
figure('Renderer', 'painters', 'Position', [10 10 740 740])
subplot(2,2,1)
imagesc(squeeze(nanmean(imp2_RDMs,1)),[0 2]); axis xy
set(gca,'Xtick',1:60)
set(gca,'Ytick',1:60)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
pbaspect([1 1 1])
colormap(flipud(autumn))
set(gca,'TickDir','out')
c=colorbar;
c.Label.String = 'squared mahalanobis distance';
c.Label.FontSize = 13;
title('Impulse 2, RDM')
%
subplot(2,2,2)
pos=[1 2.5 3.5];
hold all
plot([0 5],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([beta_side,beta_rot_cond_s,beta_rot_cond_d],...
    'positions',pos,'Widths',0.35,'Symbol','','Labels',{'cued side','rot. c., within','rot. c., across'});
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')

set(gca,'FontSize',10,'XTickLabelRotation',-40)
set(b1,'LineWidth', 1.5);
set(b1(:,1),'color',side_color);set(b1(:,2),'color',rot_cond_s_color);set(b1(:,3),'color',rot_cond_d_color);
plot(pos(1),mean(beta_side,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(2),mean(beta_rot_cond_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(3),mean(beta_rot_cond_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot([pos(1) pos(1)],beta_side_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(2) pos(2)],beta_rot_cond_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(3) pos(3)],beta_rot_cond_d_ci','Color',[0 0 0],'LineWidth',4)

s=swarmchart(pos(1).*ones(1,length(beta_side)),beta_side,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[side_color],'MarkerEdgeColor',[side_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(2).*ones(1,length(beta_side)),beta_rot_cond_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[rot_cond_s_color],'MarkerEdgeColor',[rot_cond_s_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(3).*ones(1,length(beta_side)),beta_rot_cond_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[rot_cond_d_color],'MarkerEdgeColor',[rot_cond_d_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

if p_side<0.05
    plot([pos(1) pos(1)],max(beta_side)+max(beta_side)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_rco_s<0.05
    plot([pos(2) pos(2)],max(beta_rot_cond_s)+max(beta_rot_cond_s)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_rco_d<0.05
    plot([pos(3) pos(3)],max(beta_rot_cond_d)+max(beta_rot_cond_d)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
set(gca,'TickDir','out')
ylim([-.04 .12])
xlim([.5 4])
ylabel('beta')
%
subplot(2,2,3:4)
pos=[1 2 3.5 4.5 6 7];
hold all
plot([0 8],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([beta_cued_s,beta_cued_d,beta_rot_s,beta_rot_d,beta_cued_rot_gen_s,beta_cued_rot_gen_d],...
    'positions',pos,'Widths',0.3,'Symbol','k.','Labels',{'cued, within','cued, across',...
    'rotated, within','rotated, across','gen., within','gen., across'});
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(gca,'FontSize',10,'XTickLabelRotation',-40)

set(b1(:,1),'color',cued_s_color);set(b1(:,2),'color',cued_d_color);set(b1(:,3),'color',rot_s_color);
set(b1(:,4),'color',rot_d_color);set(b1(:,5),'color',gen_s_color);set(b1(:,6),'color',gen_d_color);
set(b1,'LineWidth', 1.5);
plot(pos(1),mean(beta_cued_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(2),mean(beta_cued_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(3),mean(beta_rot_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(4),mean(beta_rot_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(5),mean(beta_cued_rot_gen_s,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(6),mean(beta_cued_rot_gen_d,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)

plot([pos(1) pos(1)],beta_cued_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(2) pos(2)],beta_cued_d_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(3) pos(3)],beta_rot_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(4) pos(4)],beta_rot_d_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(5) pos(5)],beta_cued_rot_gen_s_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(6) pos(6)],beta_cued_rot_gen_d_ci','Color',[0 0 0],'LineWidth',4)

s=swarmchart(pos(1)*ones(1,length(beta_side)),beta_cued_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',cued_s_color,'MarkerEdgeColor',cued_s_color);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(2)*ones(1,length(beta_side)),beta_cued_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',cued_d_color,'MarkerEdgeColor',cued_d_color);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(3)*ones(1,length(beta_side)),beta_rot_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',rot_s_color,'MarkerEdgeColor',rot_s_color);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(4)*ones(1,length(beta_side)),beta_rot_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',rot_d_color,'MarkerEdgeColor',rot_d_color);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(5)*ones(1,length(beta_side)),beta_cued_rot_gen_s,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',gen_s_color,'MarkerEdgeColor',gen_s_color);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(6)*ones(1,length(beta_side)),beta_cued_rot_gen_d,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',gen_d_color,'MarkerEdgeColor',gen_d_color);
s.XJitter = 'density';
s.XJitterWidth = .2;

if p_c_s<0.05
    plot([pos(1) pos(1)],max(beta_cued_s)+max(beta_cued_s)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_c_d<0.05
    plot([pos(2) pos(2)],max(beta_cued_d)+max(beta_cued_d)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_r_s<0.05
    plot([pos(3) pos(3)],max(beta_rot_s)+max(beta_rot_s)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_r_d<0.05
    plot([pos(4) pos(4)],max(beta_rot_d)+max(beta_rot_d)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_g_s<0.05
    plot([pos(5) pos(5)],max(beta_cued_rot_gen_s)+max(beta_cued_rot_gen_s)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
if p_g_d<0.05
    plot([pos(6) pos(6)],max(beta_cued_rot_gen_d)+max(beta_cued_rot_gen_d)/8,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
ylim([-.04 .08])
xlim([.5 7.5])

set(gca,'TickDir','out')
ylabel('beta')

axes('pos',[.617 .93 .07 .07])
imagesc([.1 .1], [.5 .5],cued_loc_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]); 
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.743 .93 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(rot_cond_models_same,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.829 .93 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(rot_cond_models_diff,3))); axis xy
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
imagesc([.1 .1], [.5 .5],squeeze(mean(rp_same_model,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.537 .457 .07 .07])
imagesc([.1 .1], [.5 .5],squeeze(mean(rp_diff_model,3))); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.705 .457 .07 .07])
imagesc([.1 .1], [.5 .5],cued_rp_same_model+rp_cued_same_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.815 .457 .07 .07])
imagesc([.1 .1], [.5 .5],cued_rp_diff_model+rp_cued_diff_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

annotation('textbox',[.05 1 .0 .0],'String','A','EdgeColor','none','FontSize',24)
annotation('textbox',[.55 1 .0 .0],'String','B','EdgeColor','none','FontSize',24)
annotation('textbox',[.05 .55 .0 .0],'String','C','EdgeColor','none','FontSize',24)
%



