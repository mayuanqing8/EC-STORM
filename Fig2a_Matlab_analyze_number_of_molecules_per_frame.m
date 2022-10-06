%% load all the exported PALM text files in the EC or UV folder
clear all
close all
clc

temp=dir('*.txt');
all_files=cell(length(temp),2);
for j=1:length(temp)
    filename=temp(j).name;
    all_files{j,1}=importdata(filename);
    all_files{j,2}=temp(j).name;
end
for j=1:length(temp)
all_files{j,1}.data(:,[3,4,11,13])=[];
all_files{j,1}.textdata=string(all_files{j,1}.textdata);
all_files{j,1}.textdata([3,4,11,13])=[];
all_files{j,1}.textdata(2)='frame detected';
end

%% n should be set to the first index number of file of 
%interest in the builded cell array.
for x=1:3:21
for s=x:x+2

STORM=all_files{s,1}.data;

count_per_frame=zeros(max(STORM(:,2)),1);
count=0;
n=1;
for i=1:length(STORM)-1
    if STORM(i,2)==STORM(i+1,2)
        count=count+1;
        count_per_frame(n)=count+1;
    else n=n+1;
        count=0;
    end
end


figure (x);

plot(count_per_frame);
hold on
ylabel('number of molecules per frame');
xlabel('frame number');


figure (x+21);

histogram(STORM(:,6),30);
hold on 


xlim([0,150]);
ylabel('frequency');
xlabel('photons emitted per molecule');
end
hold off
end
%% save the figures
FolderName = tempdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end
