%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Screening Chapter %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Most of this has to be run in the D:\hctsa_quiesc_cluster\hctsa-master
%environment

%%%%
%Manual experiments
%%%%

%These are mostly simple bar charts with data manually entered from my
%notebooks.

N2Hs=[71,73,72,73,73,66,85]; %
N2Os=[44,47,31,20,31,57,40]; %

npr1H=[50,33]; %14,18=32
npr1O=[18,8]; %11,12=23
%The following one was only tested once: daf-3
HawH=[63,62]; %11,13=24
HawO=[43,27]; %7,11=18
%And now dbl-1
unc129H=[40,33]; %15,15=30
unc129O=[13,0]; %15,13=28
tig2H=[62,60]; %13,15=28
tig2O=[14,20]; %14,15=29

x1=1:6;
y1=[N2Hs(1),N2Os(1),77,60,0,0];
x2=1:14;
y2=[mean(N2Hs(2:end)),mean(N2Os(2:end)),mean(npr1H),mean(npr1O),74,15,mean(HawH),mean(HawO),0,0,...
    mean(unc129H),mean(unc129O),mean(tig2H),mean(tig2O)];
figure
bar(x1,y1)
figure
bar(x2,y2)

%We also plot the emergence time between N2 on HB101, npr-1 and unc-129.
N2Htimes=[12,15,15,15,13,21,22,10,30,20,26,...
    20,15,17.5,30,30,21.5,21.5,10,12,30,30,22.5,15,...
    20,10,30,10,24.5,30,13,16,21.5,30,...
    30,30,12,19,17,21,27,30,27,12,22]; 

npr1Htimes=[15,15,20,10,30,30,14,...
    24,11,22,12,20,10]; 

unc129Htimes=[20,10,20.5,10,13,30,...
    12.5,12,20,14.5,18];

forboxplot=[N2Htimes';npr1Htimes';unc129Htimes'];
forboxplot(:,2)=[ones(length(N2Htimes),1);ones(length(npr1Htimes),1)*2;ones(length(unc129Htimes),1)*3];
figure
boxplot(forboxplot(:,1),forboxplot(:,2))

%%%%
%Power analysis
%%%%
%Power analysis will have two steps, one using simulated data through the
%GWAS pipeline to see the p value/variance explained returned corresponding
%to different effect sizes.
%The second step will be a more traditional method calculating the sample
%size necessary for a substantial number of ad hoc tests when looking for
%effects that have sizes that are findable using GWAS.

%The first step will make use of a spared down version of SNPs that the
%GWAS script uses.

load('snps_columns.mat')
load('mapping_snps_array.mat')

%I'll choose 6 random locations (one from each chromosome) and 4 effect 
%sizes (small, medium, large and very large).
chr_ids=[1,1053,2560,3956,5433,7082,8370]; %These correspond to the first SNP on
                        %each chromosome in mapping_snps_array.mat, plus
                        %the last value.
random_loc=zeros(6,1);
rng(3)
for i=1:6
    tmp_random_loc=randperm((chr_ids(i+1)-chr_ids(i)),1);
    random_loc(i,1)=tmp_random_loc+chr_ids(i)-1;
end
effect_sizes=[0.2,0.5,0.8,1.2];

load('corresponding_columns.mat') %This contains each strain that we have in
    %our dataset and the corresponding location in snps_columns. This will
    %help in identifying which strain needs to have a simulated real value
    %added to it.

%As the values will have to be normalised, an HCTSA file must be created.
timeSeriesData=cell(194*20,1);
labels=cell(194*20,1);
keywords=cell(194*20,1);
n=1;
for i=1:194
    for j=1:20
        timeSeriesData{n,1}=randn(100,1); %This is just so that it's not empty.
        labels{n,1}=['ts_',int2str(n),'_',corresponding_columns{i,1}];
        keywords{n,1}=corresponding_columns{i,1};
        n=n+1;
    end
end
save('HCTSA_simul.mat','timeSeriesData','labels','keywords')
%Starting here, the script must be run in an hctsa environment.
startup
TS_init('HCTSA_simul.mat',[],[],[false,false,false],'HCTSA_simul_init.mat');
TS_subset('HCTSA_simul_init.mat',[],1:4,1,'HCTSA_simul_init2.mat');

simulated_features=cell(length(random_loc),1);
simulated_summaries=cell(length(random_loc),1);
rng(3)
for i=1:length(random_loc)
    simulated_summaries{i,1}=cell(length(effect_sizes),1);
    load('HCTSA_simul_init2.mat');
    TS_DataMattmp=randn(size(TS_DataMat,1),1); %This way, only the effect size changes later.
    for l=1:size(TS_DataMat,1)
        l2=1;
        while l2<195 && ~strcmp(TimeSeries(l).Keywords,corresponding_columns{l2,1})
            l2=l2+1;
        end
        tmp_mod=(abs(str2double(mapping_snps_array{random_loc(i),corresponding_columns{l2,2}}))+...
            str2double(mapping_snps_array{random_loc(i),corresponding_columns{l2,2}}))/2;
        for j2=1:size(TS_DataMat,2)
            TS_DataMat(l,j2)=TS_DataMattmp(l)+effect_sizes(j2)*tmp_mod;
        end
    end
    save('HCTSA_simul_init2.mat','fromDatabase','gitInfo','MasterOperations','Operations','TimeSeries','TS_CalcTime','TS_DataMat','TS_Quality')
    TS_LabelGroups([],'HCTSA_simul_init2.mat');
    TS_normalize('scaledRobustSigmoid',[0.01,1.0],'HCTSA_simul_init2.mat');
    load('HCTSA_simul_init2_N.mat','TS_DataMat')
    load('HCTSA_simul_init2_N.mat','groupNames')
    simulated_features{i,1}=cell(194,1);
    n=1;
    for k=1:194
        simulated_features{i,1}{k,1}=groupNames{1,k};
        simulated_features{i,1}{k,2}=zeros(20,length(effect_sizes));
        for k2=1:20
            simulated_features{i,1}{k,2}(k2,:)=TS_DataMat(n,:);
            n=n+1;
        end
    end
    for j=1:length(effect_sizes)
        simulated_summaries{i,1}{j,1}=cell(194,1);
    for k=1:194
        simulated_summaries{i,1}{j,1}{k,1}=corresponding_columns{k,1};
        simulated_summaries{i,1}{j,1}{k,2}=nanmean(simulated_features{i,1}{k,2}(:,j));
        simulated_summaries{i,1}{j,1}{k,3}=nanmedian(simulated_features{i,1}{k,2}(:,j));
        simulated_summaries{i,1}{j,1}{k,4}=prctile(simulated_features{i,1}{k,2}(:,j),90);
    end
    end
end

%Each is then saved in a neatly useable version for the CeNDR GWAS
%platform.
for j=1:length(random_loc)
    fid = fopen(['simulated_effects_norm_chr_',int2str(j),'.csv'], 'w') ;
    for i=1:194
        for l=1:length(effect_sizes)-1
            fprintf(fid, '%s,', simulated_summaries{j,1}{l,1}{i,1:end-1}) ;
            fprintf(fid, '%s,', simulated_summaries{j,1}{l,1}{i,end}) ;
        end
        fprintf(fid, '%s,', simulated_summaries{j,1}{length(effect_sizes),1}{i,1:end-1}) ;
        fprintf(fid, '%s\n', simulated_summaries{j,1}{length(effect_sizes),1}{i,end}) ;
    end
    fclose(fid) ;
end

pvalues_all={[6.36,6.08,6.36,6.7];...
             [6.77,12.02,20.59,23.41,9.71,14.44,12.57];...
             [22.84,36.75,41.85,9.37,13.75,22.48,27.9,19.69];...
             [52.5,29.76,36.39,56.15,63.02,15.03,20.35,36.07,44.82]};
variances_all={[12.66,9.81,11.67,12.75];...
               [12.95,21.95,30.42,31.41,13.4,15.25,13.74];...
               [35.15,40.39,40.91,14.48,19.75,20.05,20.09,17.81];...
               [28.96,40.2,41.27,43.2,42.19,16.97,21.31,21.37,21.4]};

full_length=0;
for i=1:4
    full_length=full_length+length(pvalues_all{i,1});
end
for_hist_pv=zeros(full_length,2);
n=1;
for i=1:4
    for j=1:length(pvalues_all{i,1})
        for_hist_pv(n,1)=pvalues_all{i,1}(j);
        for_hist_pv(n,2)=i;
        n=n+1;
    end
end
figure
boxplot(for_hist_pv(:,1),for_hist_pv(:,2))

bars_means=zeros(4,1);
for i=1:4
    bars_means(i,1)=median(pvalues_all{i,1});
end
figure
bar(bars_means)

full_length=0;
for i=1:4
    full_length=full_length+length(variances_all{i,1});
end
for_hist_var=zeros(full_length,2);
n=1;
for i=1:4
    for j=1:length(variances_all{i,1})
        for_hist_var(n,1)=variances_all{i,1}(j);
        for_hist_var(n,2)=i;
        n=n+1;
    end
end
figure
boxplot(for_hist_var(:,1),for_hist_var(:,2))

combinedhist=zeros(length(full_length)*2,2);
n=0;
for i=1:4
    for j=1:length(pvalues_all{i,1})
        n=n+1;
        combinedhist(n,1)=pvalues_all{i,1}(j);
        combinedhist(n,2)=(i-1)*2+1;
        combinedhist(n+length(pvalues_all{i,1}),1)=variances_all{i,1}(j);
        combinedhist(n+length(pvalues_all{i,1}),2)=i*2;
    end
    n=n+length(pvalues_all{i,1});
end

figure
boxplot(combinedhist(:,1),combinedhist(:,2))

bars_medians=zeros(8,1);
for i=1:4
    bars_medians((i-1)*2+1,1)=median(pvalues_all{i,1});
    bars_medians(i*2,1)=median(variances_all{i,1});
end
figure
bar(bars_medians)

%For the sake of simplicity, I use the regular midbody speed for the power
%analysis.

%First, I confirm that the following conditions make no difference:
% a) day of the experiment
% b) camera
% c) time of day, by the hour.
%This is done by comparing the N2 worms on HB101 across conditions.

load('sortedfilename.mat')
load('sorted_mbspeed_reg.mat')

%Rearranging by the day.

by_day_N2H=cell(1,2);
by_day_N2H{1,1}='2017_06_10';
m=1;
by_day_N2H{1,2}=cell(1,1);
by_day_N2H{1,2}{1,1}=trajectorysorted_midspeedreg{141,2}{1,1};
n=1;
for i=2:451
    tmp_reg=regexpi(trajectorysorted_filename{141,2}{i,1},'\');
    if strcmp(by_day_N2H{m,1},trajectorysorted_filename{141,2}{i,1}(tmp_reg(7)+1:tmp_reg(8)-1))
        n=n+1;
        by_day_N2H{m,2}{n,1}=trajectorysorted_midspeedreg{141,2}{i,1};
    else
        m=m+1;
        by_day_N2H{m,1}=trajectorysorted_filename{141,2}{i,1}(tmp_reg(7)+1:tmp_reg(8)-1);
        n=1;
        by_day_N2H{m,2}{n,1}=trajectorysorted_midspeedreg{141,2}{i,1};
    end
end
for i=1:size(by_day_N2H,1)
    by_day_N2H{i,3}=zeros(length(by_day_N2H{i,2}),1);
    for j=1:length(by_day_N2H{i,2})
        by_day_N2H{i,3}(j,1)=nanmean(by_day_N2H{i,2}{j,1});
    end
end
for_anova_by_day=zeros(1,1);
tmp_titles_cells=cell(1,1);
n=1;
for i=1:size(by_day_N2H,1)
    for j=1:length(by_day_N2H{i,3})
        for_anova_by_day(n,1)=by_day_N2H{i,3}(j,1);
        tmp_titles_cells{n,1}=[by_day_N2H{i,1}(6:7),'-',by_day_N2H{i,1}(9:end)];
        n=n+1;
    end
end
p1=anova1(for_anova_by_day,tmp_titles_cells,'off'); %p=0.0567
figure
boxplot(for_anova_by_day,tmp_titles_cells)

%Rearranging by the camera.

by_cam_N2H={'1';
            '2';
            '3';
            '4';
            '5';
            '6'};
n=ones(6,1);
for i=1:451
    tmp_reg=regexpi(trajectorysorted_filename{141,2}{i,1},'\');
    out_of_3=trajectorysorted_filename{141,2}{i,1}(tmp_reg(10)-1);
    tmp_ch=regexpi(trajectorysorted_filename{141,2}{i,1},'Ch');
    chn=trajectorysorted_filename{141,2}{i,1}(tmp_ch(2)+2);
    switch out_of_3
        case '1'
            if strcmp(chn,'1')
                by_cam_N2H{1,2}{n(1),1}=trajectorysorted_midspeedreg{141,2}{i,1};
                n(1)=n(1)+1;
            else
                by_cam_N2H{2,2}{n(2),1}=trajectorysorted_midspeedreg{141,2}{i,1};
                n(2)=n(2)+1;
            end
        case '2'
            if strcmp(chn,'1')
                by_cam_N2H{3,2}{n(3),1}=trajectorysorted_midspeedreg{141,2}{i,1};
                n(3)=n(3)+1;
            else
                by_cam_N2H{4,2}{n(4),1}=trajectorysorted_midspeedreg{141,2}{i,1};
                n(4)=n(4)+1;
            end
        case '3'
            if strcmp(chn,'1')
                by_cam_N2H{5,2}{n(5),1}=trajectorysorted_midspeedreg{141,2}{i,1};
                n(5)=n(5)+1;
            else
                by_cam_N2H{6,2}{n(6),1}=trajectorysorted_midspeedreg{141,2}{i,1};
                n(6)=n(6)+1;
            end
    end
end
for i=1:size(by_cam_N2H,1)
    by_cam_N2H{i,3}=zeros(length(by_cam_N2H{i,2}),1);
    for j=1:length(by_cam_N2H{i,2})
        by_cam_N2H{i,3}(j,1)=nanmean(by_cam_N2H{i,2}{j,1});
    end
end
for_anova_by_cam=zeros(1,2);
n=1;
for i=1:size(by_cam_N2H,1)
    for j=1:length(by_cam_N2H{i,3})
        for_anova_by_cam(n,1)=by_cam_N2H{i,3}(j,1);
        for_anova_by_cam(n,2)=i;
        n=n+1;
    end
end
p1=anova1(for_anova_by_cam(:,1),for_anova_by_cam(:,2),'off'); %p=0.7248
figure
boxplot(for_anova_by_cam(:,1),for_anova_by_cam(:,2))

%Rearranging by time of day

by_time_N2H={'6'
            '7';
            '8';
            '9';
            '10';
            '11';
            '12';
            '13';
            '14'};
n=ones(9,1);
for i=1:451
    tmp_reg=regexpi(trajectorysorted_filename{141,2}{i,1},'_');
    hour_sig=trajectorysorted_filename{141,2}{i,1}(tmp_reg(length(tmp_reg)-1)+1:tmp_reg(length(tmp_reg)-1)+2);
    switch hour_sig
        case '06'
            by_time_N2H{1,2}{n(1),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(1)=n(1)+1;
        case '07'
            by_time_N2H{2,2}{n(2),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(2)=n(2)+1;
        case '08'
            by_time_N2H{3,2}{n(3),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(3)=n(3)+1;
        case '09'
            by_time_N2H{4,2}{n(4),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(4)=n(4)+1;
        case '10'
            by_time_N2H{5,2}{n(5),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(5)=n(5)+1;
        case '11'
            by_time_N2H{6,2}{n(6),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(6)=n(6)+1;
        case '12'
            by_time_N2H{7,2}{n(7),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(7)=n(7)+1;
        case '13'
            by_time_N2H{8,2}{n(8),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(8)=n(8)+1;
        case '14'
            by_time_N2H{9,2}{n(9),1}=trajectorysorted_midspeedreg{141,2}{i,1};
            n(9)=n(9)+1;
    end
end
for i=1:size(by_time_N2H,1)
    by_time_N2H{i,3}=zeros(length(by_time_N2H{i,2}),1);
    for j=1:length(by_time_N2H{i,2})
        by_time_N2H{i,3}(j,1)=nanmean(by_time_N2H{i,2}{j,1});
    end
end
for_anova_by_time=zeros(1,2);
n=1;
for i=1:size(by_time_N2H,1)
    for j=1:length(by_time_N2H{i,3})
        for_anova_by_time(n,1)=by_time_N2H{i,3}(j,1);
        for_anova_by_time(n,2)=str2double(by_time_N2H{i,1});
        n=n+1;
    end
end
p1=anova1(for_anova_by_time(:,1),for_anova_by_time(:,2),'off'); %p=0.0685
figure
boxplot(for_anova_by_time(:,1),for_anova_by_time(:,2))

cohens=0.5:0.025:1.2;
alphas=zeros(100,1);
for i=1:100
    alphas(i,1)=0.05/i;
end
nout=zeros(length(cohens),length(alphas));
for i=1:length(cohens)
    for j=1:length(alphas)
        nout(i,j)=sampsizepwr('t',[mean(for_anova_by_time(:,1)) std(for_anova_by_time(:,1))],...
            (cohens(i)*std(for_anova_by_time(:,1))+mean(for_anova_by_time(:,1))),0.8,[],'Alpha',alphas(j));
    end
    i
end

%This plot is number of experiments at alpha<0.05 vs cohen's d.
%the red ares is where the sample size=20-24
colormap_full=colormap(parula(length(unique(nout))));
for i=20:24
    tmploc=find(i==sort(unique(nout)));
    colormap_full(tmploc,:)=[1,0,0];
end
figure
heatmap(nout,1:100,cohens,[],'ColorMap',colormap_full);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GWAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A quick set of boxplots can be produced based on the divergent strains
%that shows that there is indeed difference. I can report the p value from
%the ANOVA as well, just to show off a bit more.

load('sorted_mbspeed_abs.mat')

%These are the 12 divergent strains (I removed the 4 non-CeNDR strains
%manually).
divergents={trajectorysorted_midspeedabs{47,:};
            trajectorysorted_midspeedabs{13,:};
            trajectorysorted_midspeedabs{17,:};
            trajectorysorted_midspeedabs{33,:};
            trajectorysorted_midspeedabs{44,:};
            trajectorysorted_midspeedabs{48,:};
            trajectorysorted_midspeedabs{86,:};
            trajectorysorted_midspeedabs{117,:};
            trajectorysorted_midspeedabs{124,:};
            trajectorysorted_midspeedabs{127,:};
            trajectorysorted_midspeedabs{130,:};
            trajectorysorted_midspeedabs{139,:}};

divmedians=zeros(12,2);
for i=1:12
    divergents{i,3}=zeros(length(divergents{i,2}),1);
    for j=1:length(divergents{i,2})
        divergents{i,3}(j,1)=nanmean(divergents{i,2}{j,1});
    end
    divmedians(i,1)=i;
    divmedians(i,2)=nanmedian(divergents{i,3});
end
divmedians=sortrows(divmedians,2);

for_anova_by_day=zeros(1,1);
tmp_titles_cells=cell(1,1);
n=1;
for i=1:size(divergents,1)
    for j=1:length(divergents{divmedians(i,1),2})
        for_anova_by_day(n,1)=nanmean(divergents{divmedians(i,1),2}{j,1});
        tmp_titles_cells{n,1}=divergents{divmedians(i,1),1};
        n=n+1;
    end
end
figure
boxplot(for_anova_by_day,tmp_titles_cells) %This plot is used.
[p1,anov,stats]=anova1(for_anova_by_day,tmp_titles_cells); %p=9.57*10^-50
c2=multcompare(stats);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Time Series Chapter %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This short script can be used to create violin plots of the differences
%between subsets of 50 N2 on HB101 and 50 N2 on OP50.
rng(1)
random_classifier_subsets=cell(200,1);
for i=1:200
    random_classifier_subsets{i,1}=[randperm(200,50),randperm(200,50)+200]';
end
namefiles={'HCTSA_n2s_combnew_speed_extraf_combined.mat';
           'HCTSA_n2s_combnew_pc1_extraf_combined.mat';
           'HCTSA_n2s_combnew_ic1_extraf_combined.mat';
           'HCTSA_n2s_combnew_ic2_extraf_combined.mat';
           'HCTSA_n2s_combnew_nmf2_extraf_combined.mat';
           'HCTSA_n2s_combnew_nmf5_extraf_combined.mat'};

for j=1:1
    tmp_name=namefiles{j};
    for i=6:7
        TS_subset(tmp_name,random_classifier_subsets{i,1},[],1,[tmp_name,'_s.mat'])
        TS_LabelGroups([],[tmp_name,'_s.mat']);
        TS_normalize('scaledRobustSigmoid',[0,1],[tmp_name,'_s.mat']);
        TS_TopFeatures([tmp_name,'_s_N.mat'],'fast_linear','whatPlots',{'distributions'});
        %pause
    end
end

%However, it might be more interesting to use the full dataset for images.

namefiles={'HCTSA_n2s_combnew_speed_extraf_combined.mat';
           'HCTSA_n2s_combnew_pc1_extraf_combined.mat';
           'HCTSA_n2s_combnew_ic1_extraf_combined.mat';
           'HCTSA_n2s_combnew_ic2_extraf_combined.mat';
           'HCTSA_n2s_combnew_nmf2_extraf_combined.mat';
           'HCTSA_n2s_combnew_nmf5_extraf_combined.mat'};

second_names={'HCTSA_gwasextra_N2s_combnew_mbsp.mat';
           'HCTSA_gwasextra_N2s_combnew_pc1.mat';
           'HCTSA_gwasextra_N2s_combnew_ic1.mat';
           'HCTSA_gwasextra_N2s_combnew_ic2.mat';
           'HCTSA_gwasextra_N2s_combnew_nmf2.mat';
           'HCTSA_gwasextra_N2s_combnew_nmf5.mat'};
for j=1:6
    tmp_name=namefiles{j};
    load(second_names{j},'Operations')
    tmp_idsp=zeros(5,1);
    for i=1:5
        tmp_idsp(i,1)=Operations(i).ID;
    end
    TS_subset(tmp_name,[],tmp_idsp,1,[tmp_name,'_s.mat'])
    TS_LabelGroups([],[tmp_name,'_s.mat']);
    TS_normalize('scaledRobustSigmoid',[0,1],[tmp_name,'_s.mat']);
    TS_TopFeatures([tmp_name,'_s_N.mat'],'fast_linear','whatPlots',{'distributions'},'numTopFeatures',5);
    pause(1)
end
