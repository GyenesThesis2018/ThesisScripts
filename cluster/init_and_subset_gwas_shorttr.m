corelabels={'mbsp','pc1','ic1','ic2','nmf2','nmf5'};
%corelabels={'mbsp','pc1','pc2','pc3','pc4','pc5'};
tmpname=['HCTSA_gwas_shorttr_5s_',corelabels{j},'.mat'];
%TS_init(tmpname,[],[],[false,false,false],[tmpname(1:end-4),'_init.mat']);
load('selectidids_5s.mat')
%for i=1:104
%TS_subset([tmpname(1:end-4),'_init.mat'],((i-1)*40+1):i*40,N2_feat_PCs_gwas_5s{j,1},1,['HCTSA_gwas_PCs_5s_',corelabels{j},'_sub_',int2str(i),'.mat']);
%end
TS_subset([tmpname(1:end-4),'_init.mat'],2601:2640,selectedOpIDs{j,1},1,['HCTSA_gwas_shorttr_5s_',corelabels{j},'_sub_',int2str(106),'.mat']);
exit