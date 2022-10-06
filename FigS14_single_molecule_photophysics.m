clc
clear all
% cd 'E:\single molecule counting 2021 may\STORM buffer';

[storm_table,field,STORM_image]=make_STORM_image();
x_edge=linspace(1,49400,6);
y_edge=linspace(1,49400,6);
% close all
%% this is intinial attempt to fish out the nano rulers;
inter_cluster_cutoff=60;
each_molecule_frame_idx={};
for i=1:5
    for j=1:5
        ROI_X=storm_table(:,3)>x_edge(j)& storm_table(:,3)<x_edge(j+1);
        ROI_Y=storm_table(:,4)>y_edge(i)& storm_table(:,4)<y_edge(i+1);
        ROI_XY=storm_table((ROI_X & ROI_Y),1:4);
        if numel(ROI_XY(:,1))<50
            j=j+1;
        else
            inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
            tree=linkage(inter_molecule_dist);
%             subplot (1,2,1)
%             dendrogram(tree);
            cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
%            
%                subplot (1,2,2) 
%              gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx)
%               axis square;
            
            for a=1:max(cluster_idx)
                if numel(find(cluster_idx==a))<6
                    a=a+1;
                else
                    centroid=mean(ROI_XY(cluster_idx==a,3:4));
                    dist2cen=mean(pdist2(ROI_XY(cluster_idx==a,3:4),centroid));
                    if dist2cen<35 && numel(find(cluster_idx==a))>=10
%                         scatter(ROI_XY(cluster_idx==a,3),ROI_XY(cluster_idx==a,4))
%                         axis square;
%                         xlim([mean(ROI_XY(cluster_idx==a,3))-200,mean(ROI_XY(cluster_idx==a,3))+200]);
%                         ylim([mean(ROI_XY(cluster_idx==a,4))-200,mean(ROI_XY(cluster_idx==a,4))+200]);
%                         pause (0.2);
                        each_molecule_frame_idx=cat(2,each_molecule_frame_idx,ROI_XY(cluster_idx==a,2));
                    else
                        a=a+1;
                    end
                end
                
            end
        end
        
    end
end

%%
T_off=zeros(numel(each_molecule_frame_idx),1);
T_on=zeros(numel(each_molecule_frame_idx),1);
duty_cycle=zeros(numel(each_molecule_frame_idx),1);
blink_times=zeros(numel(each_molecule_frame_idx),1);
trace=[];
intensity_trace=cell(numel(each_molecule_frame_idx),1);
for m=1:numel(each_molecule_frame_idx)
    for i=each_molecule_frame_idx{1,m}(1): each_molecule_frame_idx{1,m}(end)
        if ismember(i,each_molecule_frame_idx{1,m})
            intensity_incident=1;
        else
            intensity_incident=0;
        end
        trace=cat(1,trace,intensity_incident);
    end
    intensity_trace{m}=[0;trace;0];
     aa=diff(intensity_trace{m});

%   subplot(2,1,1)
%   plot(intensity_trace{m},'*-');
%     subplot(2,1,2)
%     plot(aa,'*-');
%     pause(0.5);

    total_on_time=sum(intensity_trace{m});
    duty_cycle(m)=total_on_time/numel(intensity_trace{m});
    blink_times(m)=sum(aa==-1);
    T_on(m)=total_on_time/sum(aa==-1);
    T_off(m)=sum(intensity_trace{m}==0)/(sum(aa==1)-1);
    trace=[];
end
%%
subplot(2,2,1)
histogram(T_on*0.06,40);
set(gca,'yscale','linear');
xlabel('T-on (second)');
ylabel('Frequency');

subplot(2,2,2)
histogram(T_off(~isinf(T_off))*0.06,40);
set(gca,'yscale','linear');
xlabel('T-off (second)');
ylabel('Frequency');


subplot(2,2,3)
histogram(duty_cycle,40);
set(gca,'yscale','linear');
xlabel('duty cycle');
ylabel('Frequency');

subplot(2,2,4)
histogram(blink_times,40);
set(gca,'yscale','linear');
xlabel('blink times');
ylabel('Frequency');




