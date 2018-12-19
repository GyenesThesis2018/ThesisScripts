%This script can be used to get the figures in the paper.

%First, we call in some of the basics. This has to be done for each figure.

colours={'b','r','g','m','cyan','yellow','black'};

addpath(fullfile(pwd,'functions'));
addpath(fullfile(pwd,'data'));
addpath(fullfile(pwd,'fig_data_scripts'));

load('data_wildtype.mat') %This contains a testing set and a 
%training set of wild type shapes only.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 1
%%%%%%%

plot(DATA_testing(1,:)) %These are the shapes.
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 2
%%%%%%%

%Fig. 2A

%We take 100 random samples and run the ICA algorithm.
%Bear in mind that ICA changes slightly between different runs even with
%the same dataset.

% addpath(fullfile(pwd,'fastICA'));
% ica_resampled_eigenshapes=cell(100,1);
% 
% for i=1:100
%     rng(i) %We seed the random number generator to be reproducible.
%     randomset=randperm(12600,3000);
%     newDATA_training=DATA_training(randomset,:);
%     [~, ica_resampled_eigenshapes{i,1}, W]=fastica(newDATA_training','approach','defl', ...
%         'numOfIC', 4, 'g', 'pow3', 'finetune', 'pow3', ...
%         'stabilization', 'on', 'lasteig', 4);
%     ica_resampled_eigenshapes{i,1}=ica_resampled_eigenshapes{i,1}';
% end

%Due to the difference between runs, we also saved the sets that we used to
%allow perfect recreation of the figures. Each run produces a set of basis
%shapes with the same properties, but the ordering or the actual shape 
%might be different.

load('ica_eigenshapes.mat')
load('ica_resampled_eigenshapes.mat')

%These ones are the eigenshapes, but we are plotting the x-y
%representations.

resampling_plot(ica_eigenshapes,ica_resampled_eigenshapes);

%Fig. 2B

%Then we use the eigenshapes to extend the projections to the testing set.
%We are then going to reproduce the original shapes and see how good they
%were.

[ica_cumulative,~]=cumulative_correlation(ica_eigenshapes,DATA_testing);


%Fig. 2C

% [pca_eigenshapes, ~]=pca(DATA_testing,'NumComponents',4);

load('pca_eigenshapes'); %Although PCA does not change in amplitude between runs.

load('trajectory_data') %This contains 100 wild type worm trajectories with
                        %both the angle representation and the speed over
                        %time.

%For this plot, we pooled the trajectories of 10 worms for aesthetic
%reasons, but individual worm trajectories have similar profiles.
figure
subplot(2,2,1)
plot_2d_eigenmap(trajectory_data,ica_eigenshapes,[1,2],1:10);
subplot(2,2,2)
plot_2d_eigenmap(trajectory_data,pca_eigenshapes,[2,3],1:10);
subplot(2,2,3)
plot_2d_eigenmap(trajectory_data,ica_eigenshapes,[1,2],1:10,50);
subplot(2,2,4)
plot_2d_eigenmap(trajectory_data,pca_eigenshapes,[2,3],1:10,50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 3
%%%%%%%

%Fig. 3A

% nmf_resampled_eigenshapes=cell(100,1);
% opts = statset('MaxIter',1000,'TolFun',1e-6,'TolX',1e-4);
% for j=1:100
%     rng(j)
%     randomset=randperm(12600,3000); %Picking a random sample.
%     newDATA_training=DATA_training(randomset,:);
%     [nmf_resampled_eigenshapes{j,1}, ~]=nnmf(newDATA_training'+pi,5,...
%         'alg','mult','rep',50,'options',opts);
%     m=max(nmf_resampled_eigenshapes{j,1}); %Normalizing.
%     for i=1:5
%         nmf_resampled_eigenshapes{j,1}(:,i)=nmf_resampled_eigenshapes{j,1}(:,i)/m(i);
%     end
% end
% Just like ICA, NMF also changes slightly between runs. We have saved the
% set used in the figures.

load('nmf_eigenshapes.mat')
load('nmf_resampled_eigenshapes.mat')

%These ones are the eigenshapes, but we are plotting the x-y
%representations.

resampling_plot(nmf_eigenshapes,nmf_resampled_eigenshapes);

%Fig. 3B

figure
for i=1:5
    plot(nmf_eigenshapes(:,i),colours{i})
    hold on
end

%Fig. 3C

[nmf_cumulative,~]=cumulative_correlation(nmf_eigenshapes,DATA_testing);

%Fig. 3D

load('trajectory_data')

figure
for i=1:4
    for j=i+1:5
        subplot(5,5,((i-1)*5+j))
        plot_2d_eigenmap(trajectory_data,nmf_eigenshapes,[i,j],1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 4
%%%%%%%

%This requires access to a database of trajectories. Fig4_data.mat contains
%a processed dataset, with genolookup serving as a catalogue of the 335
%genotypes and nmfavsabs containing the average absolute amplitudes
%corresponding to the 5 eigenshapes for each of the 9947 trajectories. A
%detailed script on how these were collected can be found in
%Fig4_data_script.m.

load('Fig4_data.mat')

%We use a Wilcoxon rank sum to compare the genotypes to wild type (N2). 
%N2 corresponds to 112 in genolookup.

mannnmfhabs=zeros(335,5);
mannnmfpabs=zeros(335,5);
for i=1:335
    for j=1:5
        [mannnmfpabs(i,j), mannnmfhabs(i,j)]=ranksum(nmfavsabs(genolookup{112,3}:genolookup{112,4},j),...
            nmfavsabs(genolookup{i,3}:genolookup{i,4},j),'alpha',0.01/335/5);
        %We use a significance level of 0.01, but also take into account
        %the multiple hypotheses being tested using Bonferroni correction.
        %Lowering the significance level to 0.005 does not change any of
        %the results in Figure 4.
    end
end

eigenshape_order=[5,3,1,4,2]; %The eigenshapes are not in order. This is
%according to Figure 3A.

mut_numbers=[236,180,52]; %236 - snf-6, 180 - nlp-1, 52 - egg-5.

fortheimage=cell(5,1);
for i=1:5
    fortheimage{i,1}=[nmfavsabs(genolookup{112,3}:genolookup{112,4},eigenshape_order(i)),...
        ones(1-genolookup{112,3}+genolookup{112,4},1)*1];
    for j=1:3
    fortheimage{i,1}=[fortheimage{i,1};...
        nmfavsabs(genolookup{mut_numbers(j),3}:genolookup{mut_numbers(j),4},eigenshape_order(i)),...
        ones(1-genolookup{mut_numbers(j),3}+genolookup{mut_numbers(j),4},1)*(j+1)];
    end
end

figure
for i=1:5
    subplot(1,5,i)
    boxplot(fortheimage{i,1}(:,1),fortheimage{i,1}(:,2))
    axis([0 5 0.2 1.2])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 5
%%%%%%%

%Figure 5 uses the pre-defined eigenshapes that are all sinusoidal.

se_eigenshapes=zeros(48,4);
s=0:1:47;
s_total=48;
for i=1:4
    y=cos(i*pi*s/s_total);
    se_eigenshapes(:,i)=y;
end

%Fig. 5A

sineig_xy=eigenshapes_xy(se_eigenshapes);

%Fig. 5B

%Here we use shapes that were collated from all of the mutants that we
%have. There is one per valid video from the dataset.

load('mutantshapes.mat')

[se_cumulative,~]=cumulative_correlation(se_eigenshapes,mutant_shapes);

%Fig. 5C

%This is a histogram representing how good the sinusoidal eigenshapes are
%compared to PCA.

load('pca_eigenshapes.mat')

[~,PCA_correlation_full]=cumulative_correlation(pca_eigenshapes,mutant_shapes,0);
[~,SE_correlation_full]=cumulative_correlation(se_eigenshapes,mutant_shapes,0);


% figure
% [N,~]=histcounts(nanmean(PCA_correlation_full{4,1},2),0:0.01:1.1);
% plot(N/12600) %Normalising by the number of shapes.
% hold on
% [N,~]=histcounts(nanmean(SE_correlation_full{4,1},2),0:0.01:1.1);
% plot(N/12600) %Normalising by the number of shapes.

figure
[N,~]=histcounts(nanmean(PCA_correlation_full{4,1},2),0:0.01:1);
N=N/9964;
[N2,~]=histcounts(nanmean(SE_correlation_full{4,1},2),0:0.01:1);
N2=N2/9964;
bar([N(80:100)',N2(80:100)'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 6
%%%%%%%

%jPCA has to be performed on actual time series of data.
load('trajectory_data')

%This method is deterministic.

%Please see details of the jPCA package and the copyright information
%pertaining to it in readme.txt.
% addpath(fullfile(pwd,'jPCA'));
% addpath(fullfile(pwd,'jPCA/fromMarksLibraries'));
% addpath(fullfile(pwd,'jPCA/CircStat2010d'));
% 
% jPCA_params.softenNorm = 5;  %These are values suggested by the original script.
% jPCA_params.suppressBWrosettes = true;  %These are for convenience
% jPCA_params.suppressHistograms = true;
% jPCA_params.numPCs = 12;
% jPCA_resampled_eigenshapes=cell(100,1);
% for j=1:100
%     rng(j)
%     randomset=randperm(99,10);
%     Data_jPCA=struct;
%     for i=1:10
%         Data_jPCA(1,i).A=trajectory_data{randomset(i),1}.angleArray';
%     end
%     [Projection12, Summary12] = jPCA(Data_jPCA, [], jPCA_params);
%     jPCA_resampled_eigenshapes{j,1}=Summary12.jPCs_highD(:,1:6);
%     j
% end

load('jPCA_eigenshapes.mat')
load('jPCA_resampled_eigenshapes')

%Fig. 6A

resampling_plot(jPCA_eigenshapes,jPCA_resampled_eigenshapes)

%Fig. 6B

figure
for i=1:5
    for j=i+1:6
        subplot(6,6,((i-1)*6+j))
        plot_2d_eigenmap(trajectory_data,jPCA_eigenshapes,[i,j],100);
        %The 100th trajectory was withheld from analysis to be used here.
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%Figure 7
%%%%%%%


%This requires access to a database of trajectories. Fig7_data.mat contains
%a processed dataset, with genolookup serving as a catalogue of the 335
%genotypes and radiusesfull as a matrix containing a measure combining the
%first two eigenshapes when moving forward (first column), when moving
%in reverse (second column) and a ratio of the two (third column). A 
%detailed script on how these were collected can be found in 
%Fig7_data_script.m.

load('Fig7_data.mat')

mannjPCAh=zeros(335,3);
mannjPCAp=zeros(335,3);
for i=1:335
    for j=1:3
        if sum(isnan(radiusesfull(genolookup{i,3}:genolookup{i,4},j)))~=...
                length(radiusesfull(genolookup{i,3}:genolookup{i,4},j))
        [mannjPCAp(i,j), mannjPCAh(i,j)]=ranksum(radiusesfull(genolookup{112,3}:genolookup{112,4},j),...
            radiusesfull(genolookup{i,3}:genolookup{i,4},j),'alpha',0.01/335/5);
        %We use a significance level of 0.01, but also take into account
        %the multiple hypotheses being tested using Bonferroni correction.
        %Lowering the significance level to 0.005 does not change any of
        %the results in Figure 7B.
        end
    end
end

mut_numbers=[247,52]; %247 - tdc-1, 52 - egg-5.

fortheimage=cell(3,1);
for i=1:3
    fortheimage{i,1}=[radiusesfull(genolookup{112,3}:genolookup{112,4},i),...
        ones(1-genolookup{112,3}+genolookup{112,4},1)*1];
    for j=1:2
    fortheimage{i,1}=[fortheimage{i,1};...
        radiusesfull(genolookup{mut_numbers(j),3}:genolookup{mut_numbers(j),4},i),...
        ones(1-genolookup{mut_numbers(j),3}+genolookup{mut_numbers(j),4},1)*(j+1)];
    end
end

figure
for i=1:2
    subplot(1,3,i)
    boxplot(fortheimage{i,1}(:,1),fortheimage{i,1}(:,2))
    axis([0 4 0.5 2])
end
subplot(1,3,3)
boxplot(fortheimage{3,1}(:,1),fortheimage{3,1}(:,2))
axis([0 4 0.6 3])

%Then we can plot the radiuses_induced. The first column (as detailed in
%Fig7_data_script.m) is forward motion in N2, the second is N2 reverse
%motion, the third is tdc-1 forward and the fourth is tdc-1 reverse.

%This figure contains the values associated with induced reversals.
figure
boxplot(radiuses_induced([(16:30)';(46:60)'],1),radiuses_induced([(16:30)';(46:60)'],2))
axis([0 3 0.5 2])



