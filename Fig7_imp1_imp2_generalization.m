%% Concurrent maintenance of veridical and transformend representations in WM

% Figure 7
% Impulse 1 cued item generaliation to impulse 2 cued and rotated

% written by Michael J. Wolff, 2022

clear all;
close all;
clc;
%%
main_dir=pwd;
addpath(genpath(main_dir))

test_chans = { 'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O2','O1','Oz'}';

nreps=100;% number of repeats for random subsampling and cross-validation
nreps_null=1000;
nfolds=8;

span=5; % number of time-points to average over (5 = 10 ms for 500 Hz)
toi=[.1 .4001];

fname_in='imp_trn1_test2';
fname_in_null='imp_trn1_test2_NULL';

fname_out='NEW_imp_trn1_test2';
fname_out_null='NEW_imp_trn1_test2_NULL';

do_decoding      = false;
do_null_decoding = false;
%%
if do_null_decoding
    cos_rot_rot_trls_null=zeros(nreps_null,30);
    cos_cued_rot_trls_null=zeros(nreps_null,30);
    cos_cued_null=zeros(nreps_null,30);
end
if do_decoding||do_null_decoding
    for sub=1:30
        fprintf(['Subject ' num2str(sub) '\n']) 
        
        load(['Behav_Veri_Trans_WM_' num2str(sub) '.mat']);

        load(['Imp1_Veri_Trans_WM_' num2str(sub) '.mat']);
        incl1=setdiff(1:size(Results,1),bad_trials);

        load(['Imp2_Veri_Trans_WM_' num2str(sub) '.mat']);
        incl2=setdiff(1:size(Results,1),bad_trials);

        incl=intersect(incl1,incl2);
        
        % extract good trials, channels and time window of interest of
        % impulse 1 epoch and reformat
        dat_temp=ft_imp1.trial(incl,ismember(ft_imp1.label,test_chans),ft_imp1.time>toi(1)&ft_imp1.time<=toi(2));
        dat_temp=bsxfun(@minus,dat_temp,mean(dat_temp,3)); % take relative baseline
        dat_temp=movmean(dat_temp,span,3,'Endpoints','discard'); % downsample
        dat_temp=dat_temp(:,:,1:span:end);
        dat_imp1=reshape(dat_temp,[size(dat_temp,1),size(dat_temp,2)*size(dat_temp,3)]); %combine channel and time dimensions
        
        % extract good trials, channels and time window of interest of
        % impulse 2 epoch and reformat
        dat_temp=ft_imp2.trial(incl,ismember(ft_imp2.label,test_chans),ft_imp2.time>toi(1)&ft_imp2.time<=toi(2));
        dat_temp=bsxfun(@minus,dat_temp,mean(dat_temp,3)); % take relative baseline
        dat_temp=movmean(dat_temp,span,3,'Endpoints','discard'); % downsample
        dat_temp=dat_temp(:,:,1:span:end);
        dat_imp2=reshape(dat_temp,[size(dat_temp,1),size(dat_temp,2)*size(dat_temp,3)]); %combine channel and time dimensions        

        cue_cond=Results(incl,4); % cued location
        rot_cond=Results(incl,5); % rotation condition
        cued_rad=Results(incl,1); % cued item in radians
        uncued_rad=Results(incl,2); % uncued item in radians
        rot_rad=Results(incl,3); % rotated item in radians
        
        ntrials=size(cue_cond,1);
        
        if do_decoding
            
            cos_left=nan(nreps,sum(cue_cond==1));
            dists_left=nan(nreps,6,sum(cue_cond==1));
            dists_orig_left=nan(nreps,6,sum(cue_cond==1));
            
            cos_right=nan(nreps,sum(cue_cond==2));
            dists_right=nan(nreps,6,sum(cue_cond==2));
            dists_orig_right=nan(nreps,6,sum(cue_cond==2));
            
            cos_all=nan(ntrials,1);
            dists_all=nan(ntrials,6);
            dists_orig_all=nan(ntrials,6);
            
            for r=1:nreps
                % left
                [cos_temp,distances_temp,distances_orig_temp] = mahal_theta_kfold_basis_b(dat_imp2(cue_cond==1,:,:),...
                    cued_rad(cue_cond==1,1),nfolds,dat_imp1(cue_cond==1,:,:));
                cos_left(r,:)=cos_temp;
                dists_left(r,:,:)=distances_temp;
                dists_orig_left(r,:,:)=distances_orig_temp;
                
                % right
                [cos_temp,distances_temp,distances_orig_temp] = mahal_theta_kfold_basis_b(dat_imp2(cue_cond==2,:,:),...
                    cued_rad(cue_cond==2,1),nfolds,dat_imp1(cue_cond==2,:,:));
                cos_right(r,:)=cos_temp;
                dists_right(r,:,:)=distances_temp;
                dists_orig_right(r,:,:)=distances_orig_temp;
            end
            
            cos_all(cue_cond==1,1)=mean(cos_left,1);
            cos_all(cue_cond==2,1)=mean(cos_right,1);
            
            dists_all(cue_cond==1,:)=squeeze(mean(dists_left,1))';
            dists_all(cue_cond==2,:)=squeeze(mean(dists_right,1))';
            
            dists_orig_all(cue_cond==1,:)=squeeze(mean(dists_orig_left,1))';
            dists_orig_all(cue_cond==2,:)=squeeze(mean(dists_orig_right,1))';
            
            % reorder to test for generalization with rotated
            u_theta=unique(rot_rad);
            theta=rot_rad;
            distances=dists_orig_all';
            theta_dist=circ_dist2(u_theta',theta)';
            distance_cos=-mean(bsxfun(@times,cos(theta_dist)',distances),1); % take cosine-weigthed mean of distances
            for c=1:length(u_theta)
                temp=round(circ_dist(u_theta,u_theta(c)),4);
                temp(temp==round(pi,4))=round(-pi,4);
                [~,i]=sort(temp);
                distances(:,theta==u_theta(c))=distances(i,theta==u_theta(c));
            end
            distances=-bsxfun(@minus,distances,mean(distances,1)); % mean-centre distances
            cos_rot=distance_cos';
            dists_rot=distances';
            
            %
            cos_cued_all(sub,1)=mean(cos_all,1);
            
            cos_cued_rot_trls(sub,1)=mean(cos_all(rot_cond~=0,1),1);
            cos_rot_rot_trls(sub,1)=mean(cos_rot(rot_cond~=0,1),1);
            
            dists_cued_rot_trls(sub,:)=mean(dists_all(rot_cond~=0,:),1);
            dists_rot_rot_trls(sub,:)=mean(dists_rot(rot_cond~=0,:),1);
            
        end
        if do_null_decoding
            
            dat_imp1_left=dat_imp1(cue_cond==1,:,:);
            dat_imp1_right=dat_imp1(cue_cond==2,:,:);
            
            dat_imp2_left=dat_imp2(cue_cond==1,:,:);
            dat_imp2_right=dat_imp2(cue_cond==2,:,:);
            
            for r=1:nreps_null
                 cos_all=nan(ntrials,1);
                dists_orig_all=nan(ntrials,6);
                                
                trl_rand_left=randperm(sum(cue_cond==1,1));
                trl_rand_right=randperm(sum(cue_cond==2,1));
                
                % left
                [cos_temp_left,distances_temp_left,distances_orig_temp_left] = mahal_theta_kfold_basis_b(dat_imp2_left(trl_rand_left,:,:),...
                    cued_rad(cue_cond==1,1),nfolds,dat_imp1_left(trl_rand_left,:,:));
                
                % right
                [cos_temp_right,distances_temp_right,distances_orig_temp_right] = mahal_theta_kfold_basis_b(dat_imp2_right(trl_rand_right,:,:),...
                    cued_rad(cue_cond==2,1),nfolds,dat_imp1_right(trl_rand_right,:,:));
                
                cos_all(cue_cond==1,1)=cos_temp_left;
                cos_all(cue_cond==2,1)=cos_temp_right;
                
                dists_orig_all(cue_cond==1,:)=distances_orig_temp_left';
                dists_orig_all(cue_cond==2,:)=distances_orig_temp_right';
                
                u_theta=unique(rot_rad);
                theta=rot_rad;
                distances=dists_orig_all';
                theta_dist=circ_dist2(u_theta',theta)';
                distance_cos=-mean(bsxfun(@times,cos(theta_dist)',distances),1); % take cosine-weigthed mean of distances
                cos_rot=distance_cos';
                
                cos_rot_rot_trls_null(r,sub)=mean(cos_rot(rot_cond~=0,1),1);
                cos_cued_rot_trls_null(r,sub)=mean(cos_all(rot_cond~=0,1),1);
                cos_cued_null(r,sub)=mean(cos_all,1);
            end
        end
    end
    if exist('cos_cued_all','var')
        save(fullfile([main_dir '\results'],[fname_out '_output.mat']),'cos_cued_all','cos_cued_rot_trls','cos_rot_rot_trls','dists_cued_rot_trls','dists_rot_rot_trls')
    end
    if exist('cos_rot_rot_trls_null','var')
        save(fullfile([main_dir '\results'],[fname_out_null '_output.mat']),'cos_rot_rot_trls_null','cos_cued_rot_trls_null','cos_cued_null')
    end
end
if ~do_decoding
    load(fullfile([main_dir '\results'],[fname_in '_output.mat']))
end
if ~do_null_decoding
    load(fullfile([main_dir '\results'],[fname_in_null '_output.mat']))
end
%%
cos_cued_T=FastTtest(cos_cued_all);
cos_cued_rot_trls_T=FastTtest(cos_cued_rot_trls);
cos_rot_rot_trls_T=FastTtest(cos_rot_rot_trls);

cos_cued_NT=FastTtest(cos_cued_null');
cos_cued_rot_trls_NT=FastTtest(cos_cued_rot_trls_null');
cos_rot_rot_trls_NT=FastTtest(cos_rot_rot_trls_null');
%%
p_cued_all_trls=FastPvalue(cos_cued_T,cos_cued_NT,2);
p_cued_rot_trls=FastPvalue(cos_cued_rot_trls_T,cos_cued_rot_trls_NT,2);
p_rot_rot_trls=FastPvalue(cos_rot_rot_trls_T,cos_rot_rot_trls_NT,2);

p_diff=FastPvalue(cos_cued_rot_trls_T-cos_rot_rot_trls_T,cos_cued_rot_trls_NT-cos_rot_rot_trls_NT,2);
%%
cued_all_ci=bootci(10000,@mean,cos_cued_all);
cued_ci=bootci(10000,@mean,cos_cued_rot_trls);
rot_ci=bootci(10000,@mean,cos_rot_rot_trls);
%%
cued_color=[0/255,0/255,200/255];
gen_color=[0/255,150/255,150/255];
close all
figure('Renderer', 'painters', 'Position', [10 10 200 350])
pos=[1 2];
hold all
plot([0 3],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([cos_cued_rot_trls,cos_rot_rot_trls],'positions',pos,'Widths',0.5,'Symbol','','Labels',{'cued to cued','cued to rotated'});
set(gca,'FontSize',10,'XTickLabelRotation',-40)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(b1(:,1),'color',cued_color)
set(b1(:,2),'color',gen_color)
set(b1,'LineWidth', 1.5);
plot(pos(1),mean(cos_cued_rot_trls,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(2),mean(cos_rot_rot_trls,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot([pos(1) pos(1)],cued_ci','Color',[0 0 0],'LineWidth',3)
plot([pos(2) pos(2)],rot_ci','Color',[0 0 0],'LineWidth',3)
pbaspect([1 2 1])

s=swarmchart(pos(1)*ones(1,length(cos_cued_rot_trls)),cos_cued_rot_trls,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',cued_color,'MarkerEdgeColor',cued_color);
s.XJitter = 'density';
s.XJitterWidth = .2;
s=swarmchart(pos(2)*ones(1,length(cos_rot_rot_trls)),cos_rot_rot_trls,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',gen_color,'MarkerEdgeColor',gen_color);
s.XJitter = 'density';
s.XJitterWidth = .2;

set(gca,'TickDir','out')
ylabel('Decoding accraucy - generalization betw. before and after rotation')
ylim([-0.002 0.002])
if p_cued_rot_trls<0.05
    plot([pos(1) pos(1)],max(cos_cued_rot_trls)+max(cos_cued_rot_trls)/6,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
end
