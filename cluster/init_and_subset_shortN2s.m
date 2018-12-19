%corelabels={'mbsp','pc1','ic1','ic2','nmf2','nmf5'};
corelabels={'mbsp','pc1','pc2','pc3','pc4','pc5'};
tmpname=['HCTSA_PCs_N2_5s_',corelabels{j},'.mat'];
TS_init(tmpname,[],[],[false,false,false],[tmpname(1:end-4),'_init.mat']);
for i=1:80
TS_subset([tmpname(1:end-4),'_init.mat'],((i-1)*5+1):i*5,[],1,['HCTSA_PCs_N2_5s_',corelabels{j},'_sub_',int2str(i),'.mat']);
end
exit