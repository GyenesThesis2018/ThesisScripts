function plot_2d_eigenmap( fullarray, eigenshapes, eigenshapeids,idsoftraj, speedthreshold )
%This function takes an array of structures each containing two fields:
%-angleArray (angle representations over time)
%-speed (the speed of the worm over time),
%each effectively being a single worm trajectory.

%It gets the amplitudes corresponding to the values in 'eigenshapes' of the
%trajectories specified in 'idsoftraj'.

%Then it plots a two dimensional histogram (using a heavily modified
%version of ndhist) for the eigenshapes in 'eigenshapeids'.

%It can be set to only plot values that are at the time above a certain
%speed threshold 'speedthreshold'.

%By default the function plots the first trajectory in the array and does
%not have a speed threshold.

if nargin<5
    speedthreshold=NaN;
end

if nargin<4
    idsoftraj=1;
end

%First, we need to confirm the orientation of the matrices.

if size(eigenshapes,1)<size(eigenshapes,2) 
    %This is to reorient the matrix so that it does not
    %matter how it is put in.
    eigenshapes=eigenshapes';
end

for i=1:length(idsoftraj)
    if size(fullarray{idsoftraj(i),1}.angleArray,2)~=size(eigenshapes,1) 
        %This is to reorient the matrix so that it does not
        %matter how it is put in.
        fullarray{idsoftraj(i),1}.angleArray=fullarray{idsoftraj(i),1}.angleArray';
    end
end

%Then we do the extensions.

extension=cell(length(idsoftraj),1);
for i=1:length(idsoftraj)
    extension{i,1}=pinv(eigenshapes(:,1:size(eigenshapes,2)))*fullarray{idsoftraj(i),1}.angleArray';
end

%If there is a speed threshold, it is applied now.

if ~isnan(speedthreshold)
    for i=1:length(idsoftraj)
        extension{i,1}=extension{i,1}(:,fullarray{idsoftraj(i),1}.speed>speedthreshold);
    end
end

%For the plotting to work, this is going to have to be put into a single
%matrix.

fullextension=[];
for i=1:length(extension)
    fullextension=[fullextension,extension{i,1}(eigenshapeids,:)];
end

%Finally, the plotting.
ndhist(fullextension(1,:),fullextension(2,:));

end

