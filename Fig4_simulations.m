%% Concurrent maintenance of veridical and transformend representations in WM

% Figure 4
% Simulation of different scenarios during the mental rotation task, and
% the resulting pattern of results

% written by Michael J. Wolff, 2022

clear all
close all
clc

main_dir=pwd;
addpath(genpath(main_dir))

load('Results_all.mat');
%%
% keep these the same for each scenario
cfg.noise_mu        = 10; % sets average noise level for each "subject"
cfg.noise_sd_trl    = 3; % sd of noise level variability between trials
cfg.nchans          = 20; % number of "channels"
cfg.reps            = 100; % number of "subjects"

scenarios={'rotated_only','partial_rotation','same_coding_schemes_both','same_coding_schemes_switch','unique_coding_schemes_both','unique_coding_schemes_switch'};

fname_in='SIM_rotation_scenarios';
fname_out='NEW_SIM_rotation_scenarios';

run_sim=false;

cfg_dec.nfolds  = 8;
cfg_dec.nreps   = 100;
cfg_dec.avg     = true;
cfg_dec.null    = false;
cfg_dec.verbose = false;

if run_sim
    for isce=1:length(scenarios)
        
        scenario=scenarios{isce};
        
        fprintf(['\n Scenario ' scenario '\n']) % indicator
        
        switch scenario

            % either the cued or the rotated item is present in a trial, never both,
            % they use the same pattern
            case 'same_coding_schemes_both'
                cfg.sc_Pc   = 0; % unique cued item pattern (can also be scaled)
                cfg.sc_Pr   = 0; % unique rotation product pattern (can also be scaled)
                cfg.sc_Pb   = 1; % common pattern between items (can also be scaled)
                cfg.c_r_switch  = 0; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 0; % whether the rotation is done only partway (half way of a way in this case)
            
            % either the cued or the rotated item is present in a trial, never both,
            % they use the same pattern
            case 'same_coding_schemes_switch'
                cfg.sc_Pc   = 0; % unique cued item pattern (can also be scaled)
                cfg.sc_Pr   = 0; % unique rotation product pattern (can also be scaled)
                cfg.sc_Pb   = 1; % common pattern between items (can also be scaled)
                cfg.c_r_switch  = 1; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 0; % whether the rotation is done only partway (half way of a way in this case)
                
                % either the cued or the rotated item is present in a trial, never both,
                % they use unique patterns
            case 'unique_coding_schemes_switch'
                cfg.sc_Pc   = 1; % unique cued item pattern(can also be scaled)
                cfg.sc_Pr   = 1; % unique rotation product pattern(can also be scaled)
                cfg.sc_Pb   = 0; % common pattern between items(can also be scaled)
                cfg.c_r_switch  = 1; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 0; % whether the rotation is done only partway (half way of a way in this case)
                
                % only the rotated item is present in each trial, but is only partially rotated,
            case 'partial_rotation'
                cfg.sc_Pc   = 0; % unique cued item pattern(can also be scaled)
                cfg.sc_Pr   = 1; % unique rotation product pattern(can also be scaled)
                cfg.sc_Pb   = 0; % common pattern between items(can also be scaled)
                cfg.c_r_switch  = 0; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 1; % whether the rotation is done only partway (half way of a way in this case)
                
                % both items are present in each trial, they use unique patterns
            case 'unique_coding_schemes_both'
                cfg.sc_Pc   = 1; % unique cued item pattern(can also be scaled)
                cfg.sc_Pr   = 1; % unique rotation product pattern(can also be scaled)
                cfg.sc_Pb   = 0; % common pattern between items(can also be scaled)
                cfg.c_r_switch  = 0; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 0; % whether the rotation is done only partway (half way of the way in this case)
            case 'cued_only'
                cfg.sc_Pc   = 1; % unique cued item pattern(can also be scaled)
                cfg.sc_Pr   = 0; % unique rotation product pattern(can also be scaled)
                cfg.sc_Pb   = 0; % common pattern between items(can also be scaled)
                cfg.c_r_switch  = 0; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 0; % whether the rotation is done only partway (half way of the way in this case)
            case 'rotated_only'
                cfg.sc_Pc   = 0; % unique cued item pattern(can also be scaled)
                cfg.sc_Pr   = 1; % unique rotation product pattern(can also be scaled)
                cfg.sc_Pb   = 0; % common pattern between items(can also be scaled)
                cfg.c_r_switch  = 0; % whether either only the cued or only the rotated item is present in a trial, never both
                cfg.part_rot    = 0; % whether the rotation is done only partway (half way of the way in this case)
        end
        %%
        noise_mu = cfg.noise_mu;
        noise_sd_trl = cfg.noise_sd_trl;
        nchans = cfg.nchans;
        sc_Pc = cfg.sc_Pc;
        sc_Pr = cfg.sc_Pr;
        sc_Pb = cfg.sc_Pb;
        c_r_switch = cfg.c_r_switch;
        part_rot = cfg.part_rot;
        reps=cfg.reps;
        %%
        fprintf(' Sim:'); m1 = '';
        RDMs_temp=zeros(reps,30,30);
        for irep=1:reps
            
            Results=data_all{randi([1 30],1)}.Results;
            
            cued_rad = circ_dist(Results(:,16),0);
            rot_rad=circ_dist(Results(:,28),0);
            rot_cond=Results(:,22);
            conditions=[cued_rad,rot_rad];
            ntrials=size(Results,1);
            
            m2 = sprintf(num2str(irep));
            fprintf([m1 m2]); m1 = repmat('\b',1,length([m2]));
            %% Set noise levels
            
            noiseLevels = noise_mu + noise_sd_trl.*randn(ntrials,1);
            noiseLevels(noiseLevels<0)=0;
            noiseLevels(noiseLevels>20)=20;
            %% made up patterns
            Pca = randn(nchans,1).*sc_Pc; % two patterns, one for sine, one for cosine
            Pcb = randn(nchans,1).*sc_Pc;
            
            Pra = randn(nchans,1).*sc_Pr;
            Prb = randn(nchans,1).*sc_Pr;
            
            Pba = randn(nchans,1).*sc_Pb;
            Pbb = randn(nchans,1).*sc_Pb;
            
            if part_rot
                rot_rad_partial=circ_dist(cued_rad+circ_dist(rot_rad,cued_rad).*.5,0);
                rot_rad=rot_rad_partial;
            end
            data=zeros(ntrials,nchans);
            for trl=1:ntrials
                if c_r_switch
                    tr=round(rand(1));
                    tc=abs(tr-1);
                else
                    tr=1; tc=1;
                end
                neuralSignal = tc*Pca.*cos(cued_rad(trl)) + tc*Pcb.*sin(cued_rad(trl))...
                    + tr*Pra.*cos(rot_rad(trl)) + tr*Prb.*sin(rot_rad(trl))...
                    + tc*Pba.*cos(cued_rad(trl)) + tc*Pbb.*sin(cued_rad(trl))...
                    + tr*Pba.*cos(rot_rad(trl)) + tr*Pbb.*sin(rot_rad(trl));
                
                noise = randn(nchans,1).*noiseLevels(trl);
                data(trl,:) = neuralSignal + noise; % add some extra noise
            end
            
            cued_rad = circ_dist(Results(:,16),0);
            rot_rad=circ_dist(Results(:,28),0);
            %%
            [RDM,cond_combs] = mahal_CV_RSA(data,conditions,cfg_dec);
            
            RDMs_temp(irep,:,:)=RDM;
            
            cos_temp_cued=nan(cfg_dec.nreps,length(cued_rad(rot_cond~=0,1)));
            cos_temp_rot=nan(cfg_dec.nreps,length(rot_rad(rot_cond~=0,1)));
            
            dat_temp=data(rot_cond~=0,:);
            cued_rad_temp=cued_rad(rot_cond~=0,1);
            rot_rad_temp=rot_rad(rot_cond~=0,1);
            for r=1:cfg_dec.nreps
                cos_temp = mahal_theta_kfold_basis_b(dat_temp,cued_rad_temp,cfg_dec.nfolds);
                cos_temp_cued(r,:)=cos_temp;
                cos_temp = mahal_theta_kfold_basis_b(dat_temp,rot_rad_temp,cfg_dec.nfolds);
                cos_temp_rot(r,:)=cos_temp;
            end
            cos_temp_cued=mean(cos_temp_cued,1)';
            cos_temp_rot=mean(cos_temp_rot,1)';
            imp2_cued_rot_rhoz(isce,irep)=fisherz(corr(cos_temp_cued,cos_temp_rot));
            imp2_cued_rot_rho(isce,irep)=(corr(cos_temp_cued,cos_temp_rot));
            imp2_cued_trlw{isce,irep}=cos_temp_cued;
            imp2_rot_trlw{isce,irep}=cos_temp_rot;
            RDMs_nreps(isce,irep,:,:)=RDM;
        end
        RDMs(isce,:,:)=mean(RDMs_temp,1);
    end
    
    %% make models
    % cued item
    cued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,1),cond_combs(:,1)))));
    
    % rotation product
    rp_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,2),cond_combs(:,2)))));
    
    % cued/rp generalization
    cued_rp_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,2),cond_combs(:,1)))));
    rp_cued_model=round(circ_rad2ang(abs(circ_dist2(cond_combs(:,1),cond_combs(:,2)))));
    
    if exist('RDMs_nreps','var')
        save(fullfile([main_dir '\results'],[fname_out '_output.mat']),'RDMs_nreps','RDMs','scenarios','cfg','cfg_dec',...
            'cued_model','rp_model','cued_rp_model','rp_cued_model','imp2_cued_rot_rhoz','imp2_cued_rot_rho','imp2_cued_trlw','imp2_rot_trlw')
    end
else
    load(fullfile([main_dir '\results'],[fname_in '_output.mat']))
end
%%
for isce=1:length(scenarios)
    for irep=1:size(RDMs_nreps,2)
        RDM=squeeze(RDMs_nreps(isce,irep,:,:));
        %% cued item (regress out rotated)
        X=[];
        X(:,1)=zscore(cued_model(1:end));
        X(:,2)=1;
        X_res=[];
        X_res(:,1)=zscore(rp_model(1:end));
        X_res(:,end+1)=1;
        Y = zscore(RDM(1:end));
        b=pinv(X_res)*(Y)';
        Y_res=zscore(Y'-X_res*b);
        betas=pinv(X)*(Y_res);
        beta_cued(isce,irep)=betas(1);
        %% rotated item (regress out cued)
        X=[];
        X(:,1)=zscore(rp_model(1:end));
        X(:,2)=1;
        X_res=[];
        X_res(:,1)=zscore(cued_model(1:end));
        X_res(:,end+1)=1;
        Y = zscore(RDM(1:end));
        b=pinv(X_res)*(Y)';
        Y_res=zscore(Y'-X_res*b);
        betas=pinv(X)*(Y_res);
        beta_rot(isce,irep)=betas(1);
        %% generalization (regress out cued and rotated)
        models=[];
        models(:,:,1)=cued_rp_model+rp_cued_model;
        X=[];
        X(:,1)=zscore(cued_rp_model(1:end)+rp_cued_model(1:end));
        X(:,2)=1;
        X(:,end+1)=1;
        X_res=[];
        X_res(:,1)=zscore(cued_model(1:end));
        X_res(:,2)=zscore(rp_model(1:end));
        X_res(:,end+1)=1;
        Y = zscore(RDM(1:end));
        b=pinv(X_res)*(Y)';
        Y_res=zscore(Y'-X_res*b);
        betas=pinv(X)*(Y_res);
        beta_gen(isce,irep)=mean(betas(1),1);
    end
end
%%
for isce=1:length(scenarios)
    beta_cued_ci(isce,:)=bootci(10000,@mean,beta_cued(isce,:));
    beta_rot_ci(isce,:)=bootci(10000,@mean,beta_rot(isce,:));
    beta_gen_ci(isce,:)=bootci(10000,@mean,beta_gen(isce,:));
    imp2_cued_rot_rhoz_ci(isce,:)=bootci(10000,@mean,imp2_cued_rot_rhoz(isce,:));
end
err_beta_cued=abs(bsxfun(@minus,beta_cued_ci,mean(beta_cued,2)));
err_beta_rot=abs(bsxfun(@minus,beta_rot_ci,mean(beta_rot,2)));
err_beta_gen=abs(bsxfun(@minus,beta_gen_ci,mean(beta_gen,2)));
err_corr=abs(bsxfun(@minus,imp2_cued_rot_rhoz_ci,mean(imp2_cued_rot_rhoz,2)));
%%
for isce=1:length(scenarios)
    p_cued(isce)=GroupPermTest(beta_cued(isce,:)',100000,1,'t');
    p_rot(isce)=GroupPermTest(beta_rot(isce,:)',100000,1,'t');
    p_gen(isce)=GroupPermTest(beta_gen(isce,:)',100000,1,'t');
    
    p_corr(isce)=GroupPermTest(imp2_cued_rot_rhoz(isce,:)',100000,1,'t');
end
%%
clc
close all
iptsetpref('ImshowBorder','tight');
figure('Renderer', 'painters', 'Position', [0 0 1000 1000])
zlim=[0 80];
for isce=1:length(scenarios)
    subplot(length(scenarios),4,((isce-1)*4)+1)
    zlims=[0 max(RDMs(isce,:,:),[],'all')];
    imagesc(squeeze(RDMs(isce,:,:)),zlims); axis xy
    set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
    
    colormap(flipud(autumn))
    pbaspect([1 1 1])
    title_string=strrep(scenarios{isce},'_',' ');
    if length(title_string)>20
        space_inds=find(isspace(title_string));
        title_string1=title_string(1:space_inds(end)-1);
        title_string2=title_string(space_inds(end)+1:end);
        title({title_string1,title_string2})
    else
        title(strrep(scenarios{isce},'_',' '))
    end
    
    set(gca,'TickDir','out')
    c=colorbar;
    c.Label.String = ['squared mahal. distance'];
    c.Label.FontSize = 9;
    
    subplot(length(scenarios),4,[((isce-1)*4)+2 ((isce-1)*4)+3])
    
    x = categorical({'cued','rotated','generalization'});
    x = reordercats(x,{'cued','rotated','generalization'});
    
    if isce==length(scenarios)
        ba=bar(x,[mean(beta_cued(isce,:),2) mean(beta_rot(isce,:),2) mean(beta_gen(isce,:),2)],'FaceColor','flat', 'BarWidth', .4);
    else
        ba=bar(x,[mean(beta_cued(isce,:),2) mean(beta_rot(isce,:),2) mean(beta_gen(isce,:),2)],'FaceColor','flat', 'BarWidth', .4);
        set(gca,'XTickLabel',[]);
        
    end
    clr=[0 0 1;0 .8 0;0 .6 .6];
    ba.CData = clr;
    set(gca,'FontSize',10,'XTickLabelRotation',-40)
    hold on
    er = errorbar(x,[mean(beta_cued(isce,:),2) mean(beta_rot(isce,:),2) mean(beta_gen(isce,:),2)],...
        [err_beta_cued(isce,1),err_beta_rot(isce,1),err_beta_gen(isce,1)],[err_beta_cued(isce,2),err_beta_rot(isce,2),err_beta_gen(isce,2)]);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    ylim([-0.05 0.65])
    pbaspect([2.2 1 1])
    hold on
    if p_cued(isce)<0.025
        plot(x(1),mean(beta_cued(isce,:),2)+0.06,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
    end
    if p_rot(isce)<0.025
        plot(x(2),mean(beta_rot(isce,:),2)+0.06,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
    end
    if p_gen(isce)<0.025
        plot(x(3),mean(beta_gen(isce,:),2)+0.06,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
    end
    
    set(gca,'TickDir','out')
    ylabel('model fit (beta)')
    
    subplot(length(scenarios),4,((isce-1)*4)+4)
    bar(1,[mean(imp2_cued_rot_rhoz(isce,:),2)],'FaceColor',[.5 .5 .5], 'BarWidth', .4);
    hold on
    er = errorbar(1,[mean(imp2_cued_rot_rhoz(isce,:),2)],...
        [err_corr(isce,1)],[err_corr(isce,2)]);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    ylim([-0.08 0.08])
    xlim([.5 1.5])
    pbaspect([1 1.5 1])
    set(gca,'XTickLabel',[]);
    set(gca,'TickDir','out')
    ylabel({'mean correlation (z)';'trlw. decoding strength'})
    hold on
    if p_corr(isce)<0.025
        plot([1 1],mean(imp2_cued_rot_rhoz(isce,:),2)+0.015,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
    end
    if p_corr(isce)>0.975
        plot([1 1],mean(imp2_cued_rot_rhoz(isce,:),2)-0.015,'*','MarkerSize',9,'MarkerEdgeColor','black','LineWidth',1)
    end
end
axes('pos',[.385 .93 .06 .06])
imagesc([.1 .1], [.4 .4],cued_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.488 .93 .06 .06])
imagesc([.1 .1], [.5 .5],rp_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])

axes('pos',[.591 .93 .06 .06])
imagesc([.1 .1], [.5 .5],cued_rp_model+rp_cued_model); axis xy
pbaspect([1 1 1])
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'xtick',[]);set(gca,'ytick',[])


