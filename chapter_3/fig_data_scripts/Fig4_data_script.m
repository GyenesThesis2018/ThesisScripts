%This is the script used to collate the data found in Fig4_data.mat. This
%makes extensive use of a large dataset that might be in a very different
%format.

load('filelist.mat') %This is a cell array with the locations of every
                     %file that is to be checked.
                     
load('nmf_eigenshapes.mat') %This contains the eigenshapes.

angleproj=cell(0);
temphelp=cell(0);
infoarray=cell(0);

n=0;
for i=1:length(fileList)
    %This initial part is specific to the file structure we used.
    a2=regexpi(fileList{i,1},'\');
    try
        temphelp=load([fileList{i,1}(1,1:a2(length(a2))),'angleArray.mat']);
    end
    %
    if ~isempty(temphelp)
        %The size constraints are once again specific to our dataset. The
        %videos were 15 mins long, average framerate=25Hz.
        if length(temphelp.angleArray)>10000
            if length(temphelp.angleArray)<30000
                %
                n=n+1;
                %We only store the projection.
                angleproj{n,1}=pinv(nmf_eigenshapes)*temphelp.angleArray';
                %
                %We want to keep the filename. This once again will be
                %different in different datasets.
                temphelp=cell(0);
                a=regexpi(fileList{i,1},'results-12-05-10');
                infoarray{n,1}=fileList{i,1}(1,(a+17):(end-13));
                infoarray{n,1}=infoarray{n,1}(~isspace(infoarray{n,1}));
                %
            end
        end
    end
end

%The construction of genolookup is extremely dependent on the file
%structure that is used. We include one based on ours, but please note that
%it is heuristic.

genonames=cell(length(infoarray),1); %First, the name of the genotype must
%be extracted from the filename.

for i=1:length(infoarray) %There were three types of filenames.
    tmp2=regexpi(infoarray{i,1},'on_food');
    tmp1=regexpi(infoarray{i,1},'Grundy');
    if ~isempty(tmp1)
        genonames{i,1}=infoarray{i,1}(tmp1+7:tmp2-2);
    else
        tmp1=regexpi(infoarray{i,1},'isolates');
        if ~isempty(tmp1)
            genonames{i,1}=infoarray{i,1}(tmp1+9:tmp2-2);
        else
            genonames{i,1}='AQ2947';
        end
    end
end

%Then we create a lookup table that also includes pointers on which
%trajectories belong to which genotype.
genonumbers=zeros(length(infoarray),1);
genolookup={};
genonumbers(1,1)=1;
genolookup{1,1}=1;
genolookup{1,2}=genonames{1,1};
genolookup{1,3}=1;
genolookup{1,4}=[];
n=1;
for i=1:length(infoarray)
    if strcmp(genolookup{n,2},genonames{i,1})
        genonumbers(i,1)=genolookup{n,1};
    else
        genolookup{n,4}=i-1;
        n=n+1;
        genolookup{n,1}=n;
        genolookup{n,2}=genonames{i,1};
        genolookup{n,3}=i;
        genonumbers(i,1)=genolookup{n,1};
    end
end
genolookup{n,4}=i;

%The next step is to create the matrix with all the average absolute
%projections.

nmfavsabs=zeros(9947,5);
for i=1:9947
    for j=1:5
        nmfavsabs(i,j)=nanmean(abs(angleproj{i,1}(j,:)));
    end
end
nmfavsabs(:,6)=genonumbers;