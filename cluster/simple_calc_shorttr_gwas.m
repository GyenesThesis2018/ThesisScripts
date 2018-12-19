corelabels={'mbsp','pc1','ic1','ic2','nmf2','nmf5'};
%corelabels={'mbsp','pc1','pc2','pc3','pc4','pc5'};
TS_compute(0,[],[],'missing',['HCTSA_gwas_shorttr_5s_',corelabels{j},'_sub_',int2str(ids),'.mat']);
exit
