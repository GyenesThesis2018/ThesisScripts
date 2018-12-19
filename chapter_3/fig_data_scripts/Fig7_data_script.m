%This is the script used to collate the data found in Fig7_data.mat. This
%makes extensive use of a large dataset that might be in a very different
%format.

load('filelist.mat') %This is a cell array with the locations of every
                     %file that is to be checked.
                     
load('jPCA_eigenshapes.mat') %This contains the eigenshapes.

angleproj=cell(0);
temphelp=cell(0);
infoarray=cell(0);
wormarrayreverse=cell(0);
wormarrayforward=cell(0);

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
        %videos were 15 mins long, average framerate=30Hz.
        if length(temphelp.angleArray)>10000
            if length(temphelp.angleArray)<30000
                %
                n=n+1;
                %We only store the projection.
                tmp=pinv(jPCA_eigenshapes)*temphelp.angleArray';
                angleproj{n,1}=tmp(1:2,:);
                %
                %We want to keep the filename. This once again will be
                %different in different datasets.
                temphelp=cell(0);
                a=regexpi(fileList{i,1},'results-12-05-10');
                infoarray{n,1}=fileList{i,1}(1,(a+17):(end-13));
                infoarray{n,1}=infoarray{n,1}(~isspace(infoarray{n,1}));
                %
                %In this case, we also need information whether the worm is
                %moving forward or backwards.
                temphelp2=load(fileList{i,1},'worm');
                wormarrayforward{n,1}=temphelp2.worm.locomotion.motion.forward.frames;
                wormarrayreverse{n,1}=temphelp2.worm.locomotion.motion.backward.frames;
                %
            end
        end
    end
    i
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

%First, we need to isolate the forward values.

forward_proj=cell(length(wormarrayforward),1);
for i=1:length(wormarrayforward)
    try %Some of the worms do not have recorded forward bouts.
    if wormarrayforward{i,1}(1,1).start<1 %This is specific to our dataset,
        %some start values are smaller than one.
        wormarrayforward{i,1}(1,1).start=1;
    end
    %We group all of the projection values during forward motion in
    %individual worms together.
    forward_proj{i,1}=angleproj{i,1}(1:2,...
        wormarrayforward{i,1}(1,1).start:wormarrayforward{i,1}(1,1).end);
    for j=2:length(wormarrayforward{i,1})
        forward_proj{i,1}=[forward_proj{i,1},...
            angleproj{i,1}(1:2,wormarrayforward{i,1}(1,j).start:wormarrayforward{i,1}(1,j).end)];
    end
    forward_proj{i,1}(:,isnan(forward_proj{i,1}(1,:)))=[];
    end
end

radiusesf=zeros(length(forward_proj),2);
for i=1:length(forward_proj)
    tmp=zeros(1,length(forward_proj{i,1}));
    for j=1:size(forward_proj{i,1},2)
        tmp(1,j)=sqrt(forward_proj{i,1}(1,j)^2+forward_proj{i,1}(2,j)^2);
    end
    radiusesf(i,1)=nanmean(abs(tmp),2);
end
radiusesf(:,2)=genonumbers;


reverse_proj=cell(length(wormarrayreverse),1);
for i=1:length(wormarrayreverse)
    try %Some of the worms do not have recorded reverse bouts.
    if wormarrayreverse{i,1}(1,1).start<1 %This is specific to our dataset,
        %some start values are smaller than one.
        wormarrayreverse{i,1}(1,1).start=1;
    end
    %We group all of the projection values during forward motion in
    %individual worms together.
    reverse_proj{i,1}=angleproj{i,1}(1:2,...
        wormarrayreverse{i,1}(1,1).start:wormarrayreverse{i,1}(1,1).end);
    for j=2:length(wormarrayreverse{i,1})
        reverse_proj{i,1}=[reverse_proj{i,1},...
            angleproj{i,1}(1:2,wormarrayreverse{i,1}(1,j).start:wormarrayreverse{i,1}(1,j).end)];
    end
    reverse_proj{i,1}(:,isnan(reverse_proj{i,1}(1,:)))=[];
    end
end

radiusesr=zeros(length(reverse_proj),2);
for i=1:length(reverse_proj)
    tmp=zeros(1,length(reverse_proj{i,1}));
    for j=1:size(reverse_proj{i,1},2)
        tmp(1,j)=sqrt(reverse_proj{i,1}(1,j)^2+reverse_proj{i,1}(2,j)^2);
    end
    radiusesr(i,1)=nanmean(abs(tmp),2);
end
radiusesr(:,2)=genonumbers;

radiusesrat=zeros(length(reverse_proj),2);
for i=1:length(reverse_proj)
    try %Some of the worms have no reverse or forward bouts.
        radiusesrat(i,1)=radiusesf(i,1)/radiusesr(i,1);
        radiusesrat(i,2)=radiusesf(i,2);
    end
end

radiusesfull=[radiusesf(:,1),radiusesr(:,1),radiusesrat(:,1),genonumbers];

%%%
%For the induced reversal dataset, we created our own videos, therefore the
%analysis procedure is different.

%First we call in the skeleton and some extra data. The first in each
%corresponds to tdc-1 and the second to N2.

load('tdc_1_induced.mat')

%skeleton contains the x-y coordinates of the centreline of the worm.
%info_induced contains cell arrays with worm index and frame number data.
%The first column corresponds to the id of the worm in the video, the
%second corresponds to the frames of the video where the worms with the
%appropriate ids appear and the third column helps in indexing into the
%skeleton arrays.

%We closely followed the videos and recorded short, but exact time frames
%before and after the touch evoked response, corresponding to forward and
%reverse motion.

%The first column corresponds to the id of the worm when it was moving
%forward, the second column to the id when it reversed, the third to the
%beginning of the forward motion, the fourth when the forward motion ended
%(normally when it was touched), the fifth when the reverse motion started
%and the sixth when the reverse motion stopped. A zero in the data means
%that the value used should be the first/last (as appropriate) value
%available for the worm.

listofwormstdc=[2,12,1,330,0,455;...
    4,41,540,670,0,810;...
    43,81,1360,0,0,1610;...
    28,28,1990,2400,2420,2500;...
    28,28,2760,2990,3010,3070;...
    144,220,3340,0,0,3940;...
    239,239,0,4660,4700,4800;...
    273,273,0,6470,6490,6580;...
    274,274,6370,6570,6610,6670;...
    274,349,7000,0,0,7610;...
    349,366,7790,0,0,8220;...
    362,362,0,8450,8450,8520;...
    362,362,8700,9110,9120,9170;...
    387,387,10000,10320,10330,10410;...
    553,553,0,11900,12030,12090];

listofwormsn2=[1,1,0,340,410,490;...
    3,3,100,450,510,580;...
    3,58,610,0,0,800;...
    5,5,500,820,850,950;...
    58,58,1030,1230,1240,1300;...
    80,80,1550,1750,1780,1860;...
    217,302,1840,0,0,2210;...
    80,80,2550,2700,2770,2870;...
    427,427,0,3530,3550,3650;...
    428,520,4500,0,0,5150;...
    521,571,0,0,0,5520;...
    571,583,5750,0,0,6170;...
    520,607,6900,0,0,7050;...
    583,583,6900,7060,7100,7200;...
    607,607,7400,7580,7610,7660];

%We first extract the skeletons.
forwardskeltdc=cell(0);
reverseskeltdc=cell(0);
forwardskeln2=cell(0);
reverseskeln2=cell(0);


for i=1:size(listofwormstdc,1)
    %Forward movement in tdc-1
    tmpfor=info_induced{1,1}((info_induced{1,1}(:,1)==listofwormstdc(i,1)),:);
    if listofwormstdc(i,3)~=0;
        tmpskelb=tmpfor(tmpfor(:,2)==listofwormstdc(i,3),3);
    else
        tmpskelb=tmpfor(1,3);
    end
    if listofwormstdc(i,4)~=0;
        tmpskele=tmpfor(tmpfor(:,2)==listofwormstdc(i,4),3);
    else
        tmpskele=tmpfor(end,3);
    end
    forwardskeltdc{i,1}=skeleton{1,1}(:,:,tmpskelb:tmpskele);
    %And then reversals in tdc-1
    tmprev=info_induced{1,1}((info_induced{1,1}(:,1)==listofwormstdc(i,2)),:);
    if listofwormstdc(i,5)~=0;
        tmpskelb=tmprev(tmprev(:,2)==listofwormstdc(i,5),3);
    else
        tmpskelb=tmprev(1,3);
    end
    if listofwormstdc(i,6)~=0;
        tmpskele=tmprev(tmprev(:,2)==listofwormstdc(i,6),3);
    else
        tmpskele=tmprev(end,3);
    end
    reverseskeltdc{i,1}=skeleton{1,1}(:,:,tmpskelb:tmpskele);
end

for i=1:size(listofwormsn2,1)
    %Forward movement in N2
    tmpfor=info_induced{2,1}((info_induced{2,1}(:,1)==listofwormsn2(i,1)),:);
    if listofwormsn2(i,3)~=0;
        tmpskelb=tmpfor(tmpfor(:,2)==listofwormsn2(i,3),3);
    else
        tmpskelb=tmpfor(1,3);
    end
    if listofwormsn2(i,4)~=0;
        tmpskele=tmpfor(tmpfor(:,2)==listofwormsn2(i,4),3);
    else
        tmpskele=tmpfor(end,3);
    end
    forwardskeln2{i,1}=skeleton{2,1}(:,:,tmpskelb:tmpskele);
    %And then reversals in N2
    tmprev=info_induced{2,1}((info_induced{2,1}(:,1)==listofwormsn2(i,2)),:);
    if listofwormsn2(i,5)~=0;
        tmpskelb=tmprev(tmprev(:,2)==listofwormsn2(i,5),3);
    else
        tmpskelb=tmprev(1,3);
    end
    if listofwormsn2(i,6)~=0;
        tmpskele=tmprev(tmprev(:,2)==listofwormsn2(i,6),3);
    else
        tmpskele=tmprev(end,3);
    end
    reverseskeln2{i,1}=skeleton{2,1}(:,:,tmpskelb:tmpskele);
end

%The skeletons will now have to be transformed into angle arrays.
forwardtdcangle=cell(15,1);
for i=1:15
    for j=1:size(forwardskeltdc{i,1},3)
        [forwardtdcangle{i,1}(j,:),~]=makeAngleArray(forwardskeltdc{i,1}(1,:,j),forwardskeltdc{i,1}(2,:,j));
    end
end

reversetdcangle=cell(15,1);
for i=1:15
    for j=1:size(reverseskeltdc{i,1},3)
        [reversetdcangle{i,1}(j,:),~]=makeAngleArray(reverseskeltdc{i,1}(1,:,j),reverseskeltdc{i,1}(2,:,j));
    end
end

forwardn2angle=cell(15,1);
for i=1:15
    for j=1:size(forwardskeln2{i,1},3)
        [forwardn2angle{i,1}(j,:),~]=makeAngleArray(forwardskeln2{i,1}(1,:,j),forwardskeln2{i,1}(2,:,j));
    end
end

reversen2angle=cell(15,1);
for i=1:15
    for j=1:size(reverseskeln2{i,1},3)
        [reversen2angle{i,1}(j,:),~]=makeAngleArray(reverseskeln2{i,1}(1,:,j),reverseskeln2{i,1}(2,:,j));
    end
end

%Then we have to transform the angles into jPCA amplitudes and extract the
%anterior body oscillation measure.

forwardjPCAtdc=cell(15,1);
forwardjPCAradtdc=cell(15,1);
forwardjPCAradtdcmean=zeros(15,1);
for i=1:15
    forwardjPCAtdc{i,1}=pinv(jPCA_eigenshapes)*forwardtdcangle{i,1}';
    for j=1:size(forwardjPCAtdc{i,1},2)
        forwardjPCAradtdc{i,1}(j,1)=sqrt(forwardjPCAtdc{i,1}(1,j)^2+forwardjPCAtdc{i,1}(2,j)^2);
    end
    forwardjPCAradtdcmean(i,1)=nanmean(forwardjPCAradtdc{i,1});
end

reversejPCAtdc=cell(15,1);
reversejPCAradtdc=cell(15,1);
reversejPCAradtdcmean=zeros(15,1);
for i=1:15
    reversejPCAtdc{i,1}=pinv(jPCA_eigenshapes)*reversetdcangle{i,1}';
    for j=1:size(reversejPCAtdc{i,1},2)
        reversejPCAradtdc{i,1}(j,1)=sqrt(reversejPCAtdc{i,1}(1,j)^2+reversejPCAtdc{i,1}(2,j)^2);
    end
    reversejPCAradtdcmean(i,1)=nanmean(reversejPCAradtdc{i,1});
end


forwardjPCAn2=cell(15,1);
forwardjPCAradn2=cell(15,1);
forwardjPCAradmeann2=zeros(15,1);
for i=1:15
    forwardjPCAn2{i,1}=pinv(jPCA_eigenshapes)*forwardn2angle{i,1}';
    for j=1:size(forwardjPCAn2{i,1},2)
        forwardjPCAradn2{i,1}(j,1)=sqrt(forwardjPCAn2{i,1}(1,j)^2+forwardjPCAn2{i,1}(2,j)^2);
    end
    forwardjPCAradmeann2(i,1)=nanmean(forwardjPCAradn2{i,1});
end

reversejPCAn2=cell(15,1);
reversejPCAradn2=cell(15,1);
reversejPCAradmeann2=zeros(15,1);
for i=1:15
    reversejPCAn2{i,1}=pinv(jPCA_eigenshapes)*reversen2angle{i,1}';
    for j=1:size(reversejPCAn2{i,1},2)
        reversejPCAradn2{i,1}(j,1)=sqrt(reversejPCAn2{i,1}(1,j)^2+reversejPCAn2{i,1}(2,j)^2);
    end
    reversejPCAradmeann2(i,1)=nanmean(reversejPCAradn2{i,1});
end

%Then we save it.
radiuses_induced(:,1)=[forwardjPCAradmeann2;reversejPCAradmeann2;...
    forwardjPCAradtdcmean;reversejPCAradtdcmean];
radiuses_induced(:,2)=[ones(15,1);ones(15,1)*2;ones(15,1)*3;ones(15,1)*4];











