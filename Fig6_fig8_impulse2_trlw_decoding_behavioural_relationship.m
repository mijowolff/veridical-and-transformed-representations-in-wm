%% Concurrent maintenance of veridical and transformend representations in WM

% Figure 6 & figure 8
% Impulse 2:
% trial-wise decoding strength correlation of cued/original and rotated
% relationship between trial-wise decoding strengths and performance

% written by Michael J. Wolff, 2022

clear all;
close all;

%%
main_dir=pwd;
addpath(genpath(main_dir))

test_chans = { 'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O2','O1','Oz'}';

nreps=100;% number of repeats for random subsampling and cross-validation
nfolds=8;

span=5; % number of time-points to average over (5 = 10 ms for 500 Hz)
toi=[.1 .4001];% time-window of interest

% load in previous results
fname_in='imp2_decoding';
fname_in_null='imp2_behavioural_relationship_NULL';
% fname_in2_null='imp2_trlw_corr_NULL';

fname_out_null='NEW_imp2_behavioural_relationship_NULL';
% fname_out2_null='NEW_imp2_trlw_corr_NULL';

fname_out='NEW_imp2_decoding';

do_decoding= false;
do_null=false;

imp2_rot_trls=[];
%%
if do_decoding
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

        cued_rad=circ_dist(Results(incl,16),0);
        rot_rad=circ_dist(Results(incl,28),0);
        acc=Results(incl,4);
        RT=Results(incl,5);
        rot_cond=Results(incl,22);

        acc_r=acc(rot_cond~=0,1);
        RT_log_r=log(RT(rot_cond~=0,1));
        %%
        cos_temp_cued_r=nan(nreps,sum(rot_cond~=0));
        cos_temp_rot_r=nan(nreps,sum(rot_cond~=0));

        for r=1:nreps

            % rot trials only
            cos_temp = mahal_theta_kfold_basis_b(dat_imp2(rot_cond~=0,:),cued_rad(rot_cond~=0),nfolds);
            cos_temp_cued_r(r,:)=cos_temp;

            cos_temp = mahal_theta_kfold_basis_b(dat_imp2(rot_cond~=0,:),rot_rad(rot_cond~=0),nfolds);
            cos_temp_rot_r(r,:)=cos_temp;
        end

        cos_temp_cued_r=mean(cos_temp_cued_r,1)';
        cos_temp_rot_r=mean(cos_temp_rot_r,1)';
        %%
        imp2_rot_trls=cat(1,imp2_rot_trls,cat(2,repmat(sub,[length(cos_temp_cued_r),1]),cos_temp_cued_r,cos_temp_rot_r,...
            acc_r,RT_log_r));

    end
    save(fullfile([main_dir '\results'],[fname_out '_output.mat']),'imp2_rot_trls')

end
if ~do_decoding
    load(fullfile([main_dir '\results'],[fname_in '_output.mat']))
end
%%
if ~do_null
    load(fullfile([main_dir '\results'],[fname_in_null '_output.mat']))
end


%%
if do_null
    rhoz_NULL=nan(30,10000);
    betas_acc_rot_NULL=nan(30,10000);
    betas_acc_cued_NULL=nan(30,10000);
    betas_RT_rot_NULL=nan(30,10000);
    betas_RT_cued_NULL=nan(30,10000);
end
for sub=1:30
    fprintf(['Subject ' num2str(sub) '\n'])
    sub_data=imp2_rot_trls(imp2_rot_trls(:,1)==sub,:);
    
    incl=1:size(sub_data,1);
    
    cued_dec=sub_data(incl,2);
    rot_dec=sub_data(incl,3);
    
    acc=sub_data(incl,4);
    RT_log=sub_data(incl,5);
    
    cued_rot_rhoz(sub,1)=fisherz(corr(cued_dec,rot_dec));
    
    ntrials=size(cued_dec,1);
    if do_null
        for r=1:10000
            rhoz_NULL(sub,r)=fisherz(corr(rot_dec(randperm(ntrials)),cued_dec));
        end
    end
    
    % cued_res
    X_res=[];
    X_res(:,1)=(rot_dec);
    X_res(:,end+1)=1;
    Y = (cued_dec);
    b=pinv(X_res)*(Y);
    cued_dec_res=(Y-X_res*b);

    % rot_res
    X_res=[];
    X_res(:,1)=(cued_dec);
    X_res(:,end+1)=1;
    Y = (rot_dec);
    b=pinv(X_res)*(Y);
    rot_dec_res=(Y-X_res*b);
            
    % acc
    acc_cat=categorical(acc);
    
    B = mnrfit(cued_dec_res,categorical(acc));
    betas_acc_cued(sub,1)=-B(2);
    
    B = mnrfit(rot_dec_res,categorical(acc));
    betas_acc_rot(sub,1)=-B(2);
    
    if do_null
        for irep=1:10000
            acc_rnd=acc_cat(randperm(ntrials));
            B = mnrfit(rot_dec_res,acc_rnd);
            betas_acc_rot_NULL(sub,irep)=-B(2);
            
            B = mnrfit(cued_dec_res,acc_rnd);
            betas_acc_cued_NULL(sub,irep)=-B(2);
        end
    end
    
    % RT cued
    X=[];
    X(:,1)=(cued_dec_res);
    X(:,end+1)=1;
    Y = (RT_log);
    beta_temp=pinv(X)*(Y);
    betas_RT_cued(sub,1)=beta_temp(1);
    if do_null
        for irep=1:10000
            beta_temp=pinv(X)*Y(randperm(ntrials));
            betas_RT_cued_NULL(sub,irep)=beta_temp(1);
        end
    end
    
    % RT rot
    X=[];
    X(:,1)=(rot_dec_res);
    X(:,end+1)=1;
    Y = (RT_log);
    beta_temp=pinv(X)*(Y);
    betas_RT_rot(sub,1)=beta_temp(1);
    if do_null
        for irep=1:10000
            beta_temp=pinv(X)*Y(randperm(ntrials));
            betas_RT_rot_NULL(sub,irep)=beta_temp(1);
        end
    end
end
if do_null
    save(fullfile([main_dir '\results'],[fname_out_null '_output.mat']),...
        'betas_acc_cued_NULL','betas_acc_rot_NULL','betas_RT_cued_NULL','betas_RT_rot_NULL','rhoz_NULL')
%     save(fullfile([main_dir '\results'],[fname_out2_null '_output.mat']),...
%         'rhoz_NULL')
end
    %
rhoz_t=FastTtest(cued_rot_rhoz);
rhoz_NT=FastTtest(rhoz_NULL);

betas_acc_cued_T=FastTtest(betas_acc_cued);
betas_acc_cued_NT=FastTtest(betas_acc_cued_NULL);

betas_acc_rot_T=FastTtest(betas_acc_rot);
betas_acc_rot_NT=FastTtest(betas_acc_rot_NULL);

betas_RT_cued_T=FastTtest(betas_RT_cued);
betas_RT_cued_NT=FastTtest(betas_RT_cued_NULL);

betas_RT_rot_T=FastTtest(betas_RT_rot);
betas_RT_rot_NT=FastTtest(betas_RT_rot_NULL);

p_acc_cued=FastPvalue(betas_acc_cued_T,betas_acc_cued_NT,2);
p_acc_rot=FastPvalue(betas_acc_rot_T,betas_acc_rot_NT,2);

p_RT_cued=FastPvalue(betas_RT_cued_T,betas_RT_cued_NT,2);
p_RT_rot=FastPvalue(betas_RT_rot_T,betas_RT_rot_NT,2);

p_r=FastPvalue(rhoz_t,rhoz_NT,2);
%%
rhoz_ci=bootci(10000,@mean,cued_rot_rhoz);
%%
figure('Renderer', 'painters', 'Position', [10 10 540 320])
hold all
plot([0 2],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot(cued_rot_rhoz,'Widths',0.3);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(b1(:,1),'color',[0 0 0])
plot(1,mean(cued_rot_rhoz,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',5)
plot([1 1],rhoz_ci','Color',[0 0 0],'LineWidth',2)
s=swarmchart(1.*ones(1,length(cued_rot_rhoz)),cued_rot_rhoz,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0]);
s.XJitter = 'density';
s.XJitterWidth = .2;
pbaspect([1 2 1])
set(gca,'TickDir','out')
ylabel({'mean correlation (z)';'trlw. decoding strength'})
ylim([-0.12 0.12])
set(gca,'xticklabel',{[]})
%%

betas_acc_cued_ci=bootci(10000,@mean,betas_acc_cued);
betas_acc_rot_ci=bootci(10000,@mean,betas_acc_rot);

betas_RT_cued_ci=bootci(10000,@mean,betas_RT_cued);
betas_RT_rot_ci=bootci(10000,@mean,betas_RT_rot);
%%
cued_color=[0 0 200/255];
rot_color=[0 150/255 0];
% close all
figure('Renderer', 'painters', 'Position', [10 10 540 320])

subplot(1,2,1)
title('Accuracy')
pos=[1 2];
hold all
plot([.5 2.5],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([betas_acc_cued,betas_acc_rot],...
    'positions',pos,'Widths',0.35,'Symbol','','Labels',{'cued/orginal','rotated'});
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(gca,'FontSize',10,'XTickLabelRotation',-40)
set(b1,'LineWidth', 1.5);
set(b1(:,1),'color',cued_color);set(b1(:,2),'color',rot_color);
plot(pos(1),mean(betas_acc_cued,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(2),mean(betas_acc_rot,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot([pos(1) pos(1)],betas_acc_cued_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(2) pos(2)],betas_acc_rot_ci','Color',[0 0 0],'LineWidth',4)

s=swarmchart(pos(1).*ones(1,length(betas_acc_cued)),betas_acc_cued,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[cued_color],'MarkerEdgeColor',[cued_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(2).*ones(1,length(betas_acc_rot)),betas_acc_rot,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[rot_color],'MarkerEdgeColor',[rot_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

ylabel('beta')
ylim([-16 16])
set(gca,'TickDir','out')

subplot(1,2,2)
title('Reaction time')
pos=[1 2];
hold all
plot([.5 2.5],[0 0 ],'Color','k','LineWidth',.5,'LineStyle',':')
b1=boxplot([betas_RT_cued,betas_RT_rot],...
    'positions',pos,'Widths',0.35,'Symbol','');
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(gca,'FontSize',10,'XTickLabelRotation',-40)
set(b1,'LineWidth', 1.5);
set(b1(:,1),'color',cued_color);set(b1(:,2),'color',rot_color);
plot(pos(1),mean(betas_RT_cued,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot(pos(2),mean(betas_RT_rot,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
plot([pos(1) pos(1)],betas_RT_cued_ci','Color',[0 0 0],'LineWidth',4)
plot([pos(2) pos(2)],betas_RT_rot_ci','Color',[0 0 0],'LineWidth',4)

s=swarmchart(pos(1).*ones(1,length(betas_RT_cued)),betas_RT_cued,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[cued_color],'MarkerEdgeColor',[cued_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

s=swarmchart(pos(2).*ones(1,length(betas_RT_rot)),betas_RT_rot,9,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor',[rot_color],'MarkerEdgeColor',[rot_color]);
s.XJitter = 'density';
s.XJitterWidth = .2;

ylabel('beta')
ylim([-2 2])
set(gca,'TickDir','out')