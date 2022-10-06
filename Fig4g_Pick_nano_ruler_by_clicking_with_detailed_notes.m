% How to use it: EC ramping up and down first; then at the same location
% collect one normal STORM movie, export the corresponding ZEN PALM table.
% using Image J_open STORM movie_STD.tif (this is to find the nanoruler
% location). Load the PALM table and .tif figure, manually pick the
% nanorulers, then open the raw EC ramping up and down movie, then we get
% the blinking times, brightness, and nanoruler location information and
% the results are saved.

clc
clear all
close all
[storm_table,field,STORM_image]=make_STORM_image();
% STORM_image=int16(STORM_image);
file_name=uigetfile('*.tif');
sum_image=imread(file_name);
sum_image2=imresize(sum_image,[5120,5120],'bicubic');
STORM_image2=imresize(STORM_image(145:end-10,15:4970),[5120,5120],'bicubic');

fft_one=fft2(STORM_image2);
fft_two=fft2(double(sum_image2));

[shift, ~] = dftregistration(fft_one,fft_two,1);

% STORM_image=STORM_image(150:end-30,30:end-150);
STORM_image2= circshift(STORM_image2,[-9,10]);
im_hight=size(STORM_image, 1);
im_wdith=size(STORM_image,2);
%
% imwrite(mat2gray(STORM_image2),'STORM_CA.tif');
%%
X_cors=[];
Y_cors=[];
% split the whole storm image into 3x3 9 regions;
x_edge=linspace(1,im_wdith,5);
y_edge=linspace(1,im_hight,5);
XY_cors=[];

for i=1:4
    for j=1:4
fh=figure (i);
fh.WindowState = 'maximized';
ROI_XY=STORM_image2(y_edge(i):y_edge(i+1), x_edge(j):x_edge(j+1));
imagesc(ROI_XY);
axis square;
colormap turbo;
colorbar;
caxis([0 0.5]);
str=sprintf('this is row %d Column %d',i,j);
title(str);
X_cors_temp=1;
Y_cors_temp=1;
finish=1;
uiwait(msgbox('Click between the two dots of nanoruler'));

while finish>0
    [X_cors_temp,Y_cors_temp] = ginput(1);
    finish= (X_cors_temp>0 && Y_cors_temp>0);
    X_cors_temp=round((X_cors_temp + x_edge(j))./10);
    Y_cors_temp=round((Y_cors_temp + y_edge(i))./10);
    X_cors=cat(1,X_cors,X_cors_temp);
    Y_cors=cat(1,Y_cors,Y_cors_temp);
end
finish=1;
close(figure(i));
pause(1);
    end
end
XY_cors=cat(2,Y_cors,X_cors);
XY_cors(Y_cors<0 | X_cors<0,:)=[];


%% test if the found nanoruler coordinates are correct;
XY_swap=cat(2,XY_cors(:,2),XY_cors(:,1));

sum_image=imread(file_name);

figure (2)
imagesc(sum_image);
axis square;
colormap turbo;
colorbar;
hold on;
for i=1:size(XY_cors,1)
% Then, from the help:
rectangle('Position',[XY_swap(i,:)-2,4,4])
end
hold off

%%
% load the raw CZI file 
hearder_movie=bfopen('*.czi');
movie=hearder_movie{1}(:,1);
im_high=size(movie{1,1},1);
im_width=size(movie{1,1},2);

%%
% edge defines the molecule crop size;
edge=2;
ROIs_int=zeros(length(movie),size(XY_cors,1));
i=1;
%  weight=[0 0 0.2 0 0; 0 0.4 0.6 0.4 0; 0.2 0.6 1 0.6 0.2; 0 0.4 0.6 0.4 0; 0 0 0.2 0 0];

while i<=size(XY_cors,1)
    if (XY_cors(i,1)-edge)>=1 &&(XY_cors(i,1)+ edge) <= im_high &&(XY_cors(i,2)-edge)>=1 && (XY_cors(i,2)+ edge) <= im_high
        for j=1:length(movie)
            frame=movie{j};
            ROI=frame((XY_cors(i,1)-edge):(XY_cors(i,1)+ edge),(XY_cors(i,2)-edge):(XY_cors(i,2)+ edge));
            ROIs_int(j,i)=sum(ROI,[1, 2]);
        end
        i=i+1;
    else
        i=i+1;
        
    end
end
%%
ruler_number=size(XY_cors,1);
for i=round(ruler_number/10)*3:round(ruler_number/10)*4
plot(ROIs_int(:,i),'*-');
hold on
end
hold off
xlim([50 500]);
ylim([3.9e4,6e4]);


%%

total_ontime=zeros(ruler_number,1);
average_ontime=zeros(ruler_number,1);
blink_times=zeros(ruler_number,1);
mean_brightness=zeros(ruler_number,1);
max_brightness=zeros(ruler_number,1);

for i=1:ruler_number
    bc=smoothdata(ROIs_int(:,i),'loess',10);
    on=bc>mean(bc)*1.02;
    if sum(on)>0
        mean_brightness(i)=mean(ROIs_int(on,i));
        max_brightness(i)=max(ROIs_int(on,i));
        total_ontime(i)=sum(on);
        blink_times(i)=sum(diff(on)==1);
    else
        mean_brightness(i)=NaN;
        max_brightness(i)=NaN;
        total_ontime(i)=NaN;
        blink_times(i)=NaN;
    end
% 
% 
%     plot(ROIs_int(:,i));
%     hold on
%     plot(bc);
%     plot(mean(bc)*0.995+on*2000);
%     hold off
%     ylim([3.8e4,5.2e4])
%      pause(1);

end
%%
subplot(2,2,1)
histogram(blink_times);
title('blink-times');
subplot(2,2,2)
histogram(max_brightness);
title('max-brightness');

subplot(2,2,3)
histogram(mean_brightness);
title('mean-brightness');

subplot(2,2,4)
histogram(total_ontime);
title('total-ontime');

%%
result.blink_time=blink_times;
result.mean_brightness=mean_brightness;
result.max_brightness=max_brightness;
result.totoal_ontime=total_ontime;
result.ruler_coordinates=XY_cors;
%%
save upanddown4_analzyed.mat result;

%% reploting the reloaded data;
subplot(2,2,1)
histfit(result.blink_time,22);
title('blink-times');
subplot(2,2,2)
histogram(result.max_brightness);
title('max-brightness');

subplot(2,2,3)
histogram(result.mean_brightness);
title('mean-brightness');

subplot(2,2,4)
histogram(result.totoal_ontime);
title('total-ontime');
%%
count_per_frame=zeros(max(storm_table(:,2)),1);
count=0;
n=1;
for i=1:length(storm_table)-1
    if storm_table(i,2)==storm_table(i+1,2)
        count=count+1;
        count_per_frame(n)=count+1;
    else n=n+1;
        count=0;
    end
end
filter=smoothdata(count_per_frame,'movmean',75);
flatted=count_per_frame-filter;
plot(count_per_frame);
hold on
plot(filter)
plot(flatted)
ylabel('number of molecules per frame');
xlabel('frame number');
hold off

%%
figure (1)
subplot(2,2,1)
histogram(result1.max_brightness);
title('max-brightness');
subplot(2,2,2)
histogram(result2.max_brightness);
title('max-brightness');

subplot(2,2,3)
histogram(result3.max_brightness);
title('max-brightness');

subplot(2,2,4)
histogram(result4.max_brightness);
title('max-brightness');

figure (2)
histogram([result1.max_brightness;result2.max_brightness;result3.max_brightness;result4.max_brightness]);


