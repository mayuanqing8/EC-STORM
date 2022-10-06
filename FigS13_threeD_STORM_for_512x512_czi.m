[STORM_raw,field_name,~]=make_STORM_image_512();
[STORM_dift_corrected,~,sup_final]=make_STORM_image_512();
clear sup_final;
%%
reader = bfGetReader('35k frames.czi');
omeMeta = reader.getMetadataStore();
weight=fspecial('gaussian',7,2);
circle_red=zeros(7);
circle_red_W=ones(7);
imsize=512;
temp_green=zeros(imsize,imsize);
circle_green=zeros(7);
circle_green_W=ones(7);
r=3;
X_Cors=round(STORM_raw(:,3)/100);
Y_Cors=round(STORM_raw(:,4)/100);
border_molecules=find(X_Cors-r <=3 | Y_Cors-r<=3 | X_Cors + r >= imsize-3 | Y_Cors + r >= imsize-3);
STORM_raw(border_molecules,:)=[];
STORM_dift_corrected(border_molecules,:)=[];
X_Cors=round(STORM_raw(:,3)/100);
Y_Cors=round(STORM_raw(:,4)/100);
frame=STORM_raw(:,2);

stack_size=1000;
loaded_red=[];
loaded_green=[];
red_green=zeros(512,512,stack_size*2);
stack_number=0;
circle_ratio_total=[];

tic
while stack_number*stack_size<=35000
    read_frame_start=stack_number*stack_size*2+1;
    read_frame_end=stack_size*2+(stack_number*stack_size*2);
    m=1;
    for i=read_frame_start:read_frame_end
        red_green(:,:,m)=bfGetPlane(reader, i);
        m=m+1;
    end
    loaded_red = red_green(:,:,1:2:end-1);
    loaded_green = red_green(:,:,2:2:end);
    loaded_red=gpuArray(loaded_red);
    loaded_green=gpuArray(loaded_green);
    red_green=zeros(512,512,stack_size*2);
    temp_ind_start=find(frame==stack_number*stack_size+1);
    temp_ind_end=find(frame==stack_size+stack_number*stack_size);
    index=temp_ind_start(1):temp_ind_end(end);
    circle_red_sum=ones(1,numel(index));
    circle_green_sum=ones(1,numel(index));

mm=1;
for j=index
    temp_red=loaded_red(:,:,frame(j)-stack_number*stack_size); 
    temp_green=loaded_green(:,:,frame(j)-stack_number*stack_size);
    circle_red=temp_red((X_Cors(j)-r):(X_Cors(j)+r),(Y_Cors(j)-r):(Y_Cors(j)+r));
    circle_green=temp_green((X_Cors(j))-r:(X_Cors(j))+r,(Y_Cors(j))-r:(Y_Cors(j))+r);
    circle_red_W=circle_red.*weight; 
    circle_green_W=circle_green.*weight;
    circle_red_sum(mm)=sum(circle_red_W(:));
    circle_green_sum(mm)=sum(circle_green_W(:));
    mm=mm+1;
end

stack_number=stack_number+1;

circle_ratio_stack=(circle_red_sum-circle_green_sum)./(circle_red_sum + circle_green_sum); 
    circle_ratio_total=cat(2,circle_ratio_total,circle_ratio_stack);
end
toc
%%
STORM_dift_corrected=STORM_dift_corrected(1:numel(circle_ratio_total),:);
STORM_dift_corrected=[STORM_dift_corrected,circle_ratio_total'];
%%
 [threeD_STORM_image]=make_3D_STORM_image_new_512(STORM_dift_corrected,200,0.8 ,1.39);
 %% corse filtering
x_edge=linspace(1,imsize*100,5);
y_edge=linspace(1,imsize*100,5);
filtered_STORM=[];

for i=1:4
    for j=1:4
        ROI_X=STORM_dift_corrected(:,3)>x_edge(j)& STORM_dift_corrected(:,3)<x_edge(j+1);
        ROI_Y=STORM_dift_corrected(:,4)>y_edge(i)& STORM_dift_corrected(:,4)<y_edge(i+1);
        ROI_XY=STORM_dift_corrected((ROI_X & ROI_Y),1:10);
        cluster_idx = dbscan(ROI_XY(:,[3,4]),250,120);
%           gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx)
%                  axis square;
%         fiber=mode(cluster_idx);
        filtered_STORM=cat(1,filtered_STORM,ROI_XY(cluster_idx~=-1,1:10));

    end
end

%% 
 [Filtered_threeD_STORM_image]=make_3D_STORM_image_smooth(filtered_STORM,10,500,0.85 ,1.25);

 
%% fine filtering
x_edge=linspace(1,imsize*100,5);
y_edge=linspace(1,imsize*100,5);
filtered_STORM2=[];

for i=1:4
    for j=1:4
        ROI_X=filtered_STORM(:,3)>x_edge(j)& filtered_STORM(:,3)<x_edge(j+1);
        ROI_Y=filtered_STORM(:,4)>y_edge(i)& filtered_STORM(:,4)<y_edge(i+1);
        ROI_XY=filtered_STORM((ROI_X & ROI_Y),1:10);
        cluster_idx = dbscan(ROI_XY(:,[3,4]),50,15);
%           gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx)
%                  axis square;
%         fiber=mode(cluster_idx);
        filtered_STORM2=cat(1,filtered_STORM2,ROI_XY(cluster_idx~=-1,1:10));

    end
end
%%

 [Filtered_threeD_STORM_image]=make_3D_STORM_image_new_512(filtered_STORM2,500,0.89 ,1.2);

