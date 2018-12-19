function [ corrs_mean, corrs_full ] =cumulative_correlation( eigenshapes, testing, plotting )
%This function first gets the full projection using all of the eigenshapes.

%Then reproduces each shape in the testing set using an increasing number
%of eigenshapes (1, 1-2, 1-2-3, etc.). 

%It then compares each of these shapes with the original and plots the
%result.

%It plots the mean results by default.

if nargin<3
    plotting=1;
end

colours={'b','r','g','m','cyan','yellow','black'};

if size(eigenshapes,1)<size(eigenshapes,2) %This is to reorient the matrix so that it does not
                       %matter how it is put in.
    eigenshapes=eigenshapes';
end

if size(testing,2)~=size(eigenshapes,1) %This is to reorient the matrix so that it does not
                       %matter how it is put in.
    testing=testing';
end

corrs_full=cell(size(eigenshapes,2),1);
corrs_mean=zeros(size(eigenshapes,2),size(eigenshapes,1));

for i=1:size(eigenshapes,2)
    tmp=pinv(eigenshapes(:,1:size(eigenshapes,2)))*testing'; %Full projection
    newshapestmp=eigenshapes(:,1:i)*tmp(1:i,:); %Partial reconstruction
    corrs_full{i,1}=bsxfun(@minus, newshapestmp', testing).^2; %Comparison
    corrs_full{i,1}=1-corrs_full{i,1}; 
    corrs_mean(i,:)=nanmean(corrs_full{i,1});
end

if plotting==1
    figure
    for i=1:size(eigenshapes,2)
        plot(corrs_mean(i,:),colours{i})
        hold on
    end
    
    legendnumbers=cell(size(eigenshapes,2),1);
    for i=1:size(eigenshapes,2)
        legendnumbers{i,1}=int2str(i);
    end
    
    legend(legendnumbers,'Location','southeast')
end


end

