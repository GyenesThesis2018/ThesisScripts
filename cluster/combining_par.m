%corelabels={'mbsp','pc1','ic1','ic2','nmf2','nmf5'};
comb_names={'mbsp','pc1','pc2','pc3','pc4','pc5'};
TS_combine(['HCTSA_gwas_PCs_5s_',comb_names{j},'_sub2_',int2str(1),'.mat'],['HCTSA_gwas_PCs_5s_',comb_names{j},'_sub2_',int2str(2),'.mat'],false,false,['HCTSA_gwas_PCs_5s_',comb_names{j},'_comb2.mat']);
for i=3:105
TS_combine(['HCTSA_gwas_PCs_5s_',comb_names{j},'_comb2.mat'],['HCTSA_gwas_PCs_5s_',comb_names{j},'_sub2_',int2str(i),'.mat'],false,false,['HCTSA_gwas_PCs_5s_',comb_names{j},'_comb2.mat']);
end
exit