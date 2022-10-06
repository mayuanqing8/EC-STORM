clc
clear all
close all
temp=dir('*.txt');
for file_number=1:length(temp)
    filename=temp(file_number).name;
    raw=importdata(filename);
    raw.data(:,[3,4,11,13])=[];
    names=string(raw.textdata);
    names([3,4,11,13])=[];
    names(2)='frame detected';
    storm_table=raw.data;
    
        
    x_edge=linspace(1,49400,21);
    y_edge=linspace(1,49400,21);
    new_storm_table=[];
    total_clusters=0;
    total_nano_ruler=0;
    inter_cluster_cutoff=200;
    
    for i=1:20
        for j=1:20
            ROI_X=storm_table(:,3)>x_edge(j)& storm_table(:,3)<x_edge(j+1);
            ROI_Y=storm_table(:,4)>y_edge(i)& storm_table(:,4)<y_edge(i+1);
            ROI_XY=storm_table((ROI_X & ROI_Y),1:4);
            if numel(ROI_XY(:,1))<80
                j=j+1;
            else
                inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
                tree=linkage(inter_molecule_dist);
                %             subplot (1,2,1)
                %             dendrogram(tree);
                cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
                total_clusters=total_clusters+max(cluster_idx);
                %                subplot (1,2,2)
                %              gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx)
                %               axis square;
                %
                for a=1:max(cluster_idx)
                    if numel(find(cluster_idx==a))<10
                        a=a+1;
                    else
                        centroid=mean(ROI_XY(cluster_idx==a,3:4));
                        dist2cen=mean(pdist2(ROI_XY(cluster_idx==a,3:4),centroid));
                        if dist2cen>25 && numel(find(cluster_idx==a))>15
                            new_storm_table=cat(1,new_storm_table,ROI_XY(cluster_idx==a,1:4));
                            total_nano_ruler=total_nano_ruler+1;
                        else
                            a= a+1;
                        end
                    end
                    
                end
            end
            
        end
    end
    
%     
%     ruler_indx=ismember(storm_table(:,1),new_storm_table(:,1));
%     non_ruler=find(ruler_indx==0);
    
    x_edge=linspace(1,49400,6);
    y_edge=linspace(1,49400,6);
    all_spot_2_spot=[];
    each_spot_1_frame_ind={};
    each_spot_2_frame_ind={};
    each_spot_3_frame_ind={};
    
    intra_cluster_cutoff=50;
    
    for i=1:5
        for j=1:5
            ROI_X=new_storm_table(:,3)>x_edge(j)& new_storm_table(:,3)<x_edge(j+1);
            ROI_Y=new_storm_table(:,4)>y_edge(i)& new_storm_table(:,4)<y_edge(i+1);
            ROI_XY=new_storm_table((ROI_X & ROI_Y),1:4);
            if numel(ROI_XY(:,1))<80
                j=j+1;
            else
                inter_molecule_dist=pdist([ROI_XY(:,3),ROI_XY(:,4)]);
                tree=linkage(inter_molecule_dist);
                
                cluster_idx=cluster(tree,'Cutoff',inter_cluster_cutoff,'Criterion','distance');
%                  gscatter(ROI_XY(:,3),ROI_XY(:,4),cluster_idx);
                 
                for ruler_ind=1:max(cluster_idx)
                    if numel(find(cluster_idx==ruler_ind))<70
                        ruler_ind=ruler_ind+1;
                    else
                        ruler_cent=median(ROI_XY(cluster_idx==ruler_ind,3:4));
                        ruler_X=new_storm_table(:,3)>ruler_cent(1)-150 & new_storm_table(:,3)<ruler_cent(1)+ 150;
                        ruler_Y=new_storm_table(:,4)>ruler_cent(2)-150 & new_storm_table(:,4)<ruler_cent(2)+150 ;
                        ruler_XY=new_storm_table((ruler_X & ruler_Y),1:4);
                        if numel(ruler_XY(:,1))<10
                             ruler_cent=median(ROI_XY(cluster_idx==ruler_ind,3:4));
                             ruler_X=new_storm_table(:,3)>ruler_cent(1)-250 & new_storm_table(:,3)<ruler_cent(1)+ 250;
                             ruler_Y=new_storm_table(:,4)>ruler_cent(2)-250 & new_storm_table(:,4)<ruler_cent(2)+250 ;
                             ruler_XY=new_storm_table((ruler_X & ruler_Y),1:4);
                        else
                        intra_ruler_dist=pdist(ruler_XY(:,3:4),'cityblock');
                        tree2=linkage(intra_ruler_dist,'single');
                        molecule_idx=cluster(tree2,'Cutoff',intra_cluster_cutoff,'Criterion','distance');
                        end
                        
                        %                     subplot (1,3,1)
                        %                     gscatter(ruler_XY(:,3),ruler_XY(:,4)),molecule_idx);
                        %                     title(sprintf('ruler index is: %d,number is: %d',ruler_ind,numel(molecule_idx)));
                        %                     axis square;
                        %
                        dens_spot_1=mode(molecule_idx);
                        if numel(find(molecule_idx==dens_spot_1))/numel(molecule_idx)>0.85...
                                || numel(find(molecule_idx==dens_spot_1))/numel(molecule_idx) <0.15
                            if numel(molecule_idx)>350
                                molecule_idx=cluster(tree2,'Cutoff',intra_cluster_cutoff/6,'Criterion','distance');
                            elseif numel(molecule_idx)>150 && numel(molecule_idx)<350
                                molecule_idx=cluster(tree2,'Cutoff',intra_cluster_cutoff/3.5,'Criterion','distance');
                            else
                                molecule_idx=cluster(tree2,'Cutoff',intra_cluster_cutoff/2.5,'Criterion','distance');
                            end
                            
                            %                         subplot (1,3,2)
                            %                         gscatter(ruler_XY(:,3),ruler_XY(:,4),molecule_idx);
                            %                         axis square;
                            dens_spot_1=mode(molecule_idx);
                            if numel(find(molecule_idx==dens_spot_1))/numel(molecule_idx)>0.75...
                                    || numel(find(molecule_idx==dens_spot_1))/numel(molecule_idx) <0.25
                                ruler_ind=ruler_ind+1;
                            else
                                molecule_idx(molecule_idx==dens_spot_1)=NaN;
                                dens_spot_2=mode(molecule_idx);
                                ratio=numel(find(molecule_idx==dens_spot_2))/numel(find(isnan(molecule_idx)));
                                summ=numel(find(molecule_idx==dens_spot_2))+numel(find(isnan(molecule_idx)));
                                if ratio>0.1 && summ/numel(molecule_idx)>0.6
                                    [molecule_idx,centroids]=kmeans(ruler_XY(:,3:4),2);
                                    all_spot_2_spot=cat(2,all_spot_2_spot,pdist(centroids));
                                    each_spot_1_frame_ind=cat(2,each_spot_1_frame_ind,ruler_XY(molecule_idx==1,2));
                                    each_spot_2_frame_ind=cat(2,each_spot_2_frame_ind,ruler_XY(molecule_idx==2,2));
                                else
                                    ruler_ind=ruler_ind+1;
                                    
                                end
                            end
                            
                        else
                            [molecule_idx,centroids]=kmeans(ruler_XY(:,3:4),2);
                            all_spot_2_spot=cat(2,all_spot_2_spot,pdist(centroids));
                            each_spot_1_frame_ind=cat(2,each_spot_1_frame_ind,ruler_XY(molecule_idx==1,2));
                            each_spot_2_frame_ind=cat(2,each_spot_2_frame_ind,ruler_XY(molecule_idx==2,2));
                        end
                        %                     subplot (1,3,3)
                        %                     scatter(ruler_XY(:,3),ruler_XY(:,4),10,'.');
                        %                     axis square;
                        %                     text(centroids(1,1),centroids(1,2),sprintf('%d',1),'FontSize',13);
                        %                     text(centroids(2,1),centroids(2,2),sprintf('%d',2),'FontSize',13);
                        %                     pause(2);
                    end
                    
                end
            end
            
        end
    end
    % ploting the spot to spot distance within each nanoruler with filter;
    if numel(all_spot_2_spot)>30
    [bin,bin_edge]=histcounts(all_spot_2_spot(all_spot_2_spot>1 & all_spot_2_spot<200),15);
    xx=linspace(bin_edge(1),bin_edge(end),15);
    f=fit(xx',bin','gauss1');
    else
      f.b1='not enough nanoruler found';
    end
    
    % alternative way of counting molecules, which looks at
    % the duration and gaps between molecules systmatically
    count_per_mol_in_spot_1=[];
    count_per_mol_in_spot_2=[];
    mol_in_1st_spot=[];
    mol_in_2nd_spot=[];
    for abc=1:numel(each_spot_1_frame_ind)
        % plot(each_spot_1_frame_ind{1,abc},1:numel(each_spot_1_frame_ind{1,abc}),'*');ave_mol_spot_1=numel(find(mol_in_1st_spot>50))/numel(each_spot_1_frame_ind);
        [~, count_per_mol_1,mol_counted_1]=count_molecules(each_spot_1_frame_ind{1,abc});
        [~,count_per_mol_2,mol_counted_2]=count_molecules(each_spot_2_frame_ind{1,abc});
        count_per_mol_in_spot_1=cat(2,count_per_mol_in_spot_1,count_per_mol_1);
        count_per_mol_in_spot_2=cat(2,count_per_mol_in_spot_2,count_per_mol_2);
        mol_in_1st_spot=cat(1,mol_in_1st_spot,mol_counted_1);
        mol_in_2nd_spot=cat(1,mol_in_2nd_spot,mol_counted_2);
    end
    ave_mol_in_spot_1=median(mol_in_1st_spot);
    ave_mol_in_spot_2=median(mol_in_2nd_spot);
    ave_count_per_mol_in_spot_1=median(count_per_mol_in_spot_1);
    ave_count_per_mol_in_spot_2=median(count_per_mol_in_spot_2);
    
    %  save the result
    result(file_number).file_name=filename;
    result(file_number).total_cluster_number=total_clusters;
    result(file_number).totoal_nano_ruler_number=total_nano_ruler;
    result(file_number).percentage_of_nanoruler_first_found=total_nano_ruler/total_clusters;
    result(file_number).percentage_of_nanoruler_finally_found=numel(mol_in_1st_spot)/total_clusters;
    result(file_number).average_molecule_in_1st_spot_of_nanoruler=ave_mol_in_spot_1;
    result(file_number).all_molecules_in_1st_spot_of_nanoruler=mol_in_1st_spot;
    result(file_number).average_molecule_in_2nd_spot_of_nanoruler=ave_mol_in_spot_2;
    result(file_number).all_molecules_in_2nd_spot_of_nanoruler=mol_in_2nd_spot;
    result(file_number).average_durations_of_frames_of_1st_spot_of_nanoruler=ave_count_per_mol_in_spot_1;
    result(file_number).all_durations_of_frames_of_1st_spot_of_nanoruler=count_per_mol_in_spot_1;
    result(file_number).average_durations_of_frames_of_2nd_spot_of_nanoruler=ave_count_per_mol_in_spot_2;
    result(file_number).all_durations_of_frames_of_2nd_spot_of_nanoruler=count_per_mol_in_spot_2;
    result(file_number).all_distance_between_two_spot_within_one_nanoruler=all_spot_2_spot;
    result(file_number).average_distance_between_two_spot_within_one_nanoruler=f.b1;
end
%
%  new_filename=erase(filename,'.txt');
%  new_filename= new_filename(find(~isspace(new_filename)));
save('all_result_combined_by_kmeans.mat','result');
clear all
close all
clc
