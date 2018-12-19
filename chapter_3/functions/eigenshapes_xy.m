function [ xy ] = eigenshapes_xy( A )
%This function takes the matrix of eigenshapes and transforms them into
%xy coordinates. It also plots them.
%OUTPUT: xy is a cell containing the xy coordinates of each eigenshape.

colours={'b','r','g','m','cyan','yellow','black'}; %colour palette

if size(A,1)<size(A,2) %This is to reorient the matrix so that it does not
                       %matter how it is put in.
    A=A';
end

xy=cell(size(A,2),1);
figure
for i=1:size(A,2);
    xy{i,1}=zeros(size(A,1)+1,2);
    [xy{i,1}(:,1),xy{i,1}(:,2)]=angle2skel(A(:,i),0,1);
    if size(A,2)<8
    plot(xy{i,1}(:,1),xy{i,1}(:,2),colours{i})
    else
        plot(xy{i,1}(:,1),xy{i,1}(:,2))
    end
    hold on
end    


end

