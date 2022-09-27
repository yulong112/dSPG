function nodes_to_add = Connect_pixel_btwsp2(nodes_to_super,superpix_img, distance, Height, Width, knn_distance_all)
N=size(distance,1);
nodes_to_add = false(N);
sp_num=size(nodes_to_super,1);
center_set=zeros(sp_num,1);

for i = 1:sp_num
    pixel_set = find(superpix_img==i);
    %% method 1: using hingle nodes
    tmp_distance=distance(pixel_set,pixel_set);
%     % Number of nearest neighbors
%     knn_param = 10;
%     % Calculating distances of k-nearest neighbors
%     knn_distance = zeros(length(pixel_set),1);
%     for kk = 1:length(pixel_set)
%         pixel_loc_tmp=pixel_set(kk);
%         % sort all possible neighbors according to distance
%         if strcmp(param2,'DIST')||strcmp(param2,'DIST2')||strcmp(param2,'DIST3')
%             temp = sort(distance(pixel_loc_tmp,:), 'ascend');
%         else
%             temp = sort(distance(pixel_loc_tmp,:), 'descend');
%         end
%         % select k-th neighbor: knn_param+1, as the node itself is considered
%         knn_distance(kk) = temp(knn_param + 1);
%     end

    knn_distance=knn_distance_all(pixel_set);
    
    nodes_to_knn = true(length(pixel_set));
    for kk = 1:length(pixel_set)
        pixel_loc_tmp=pixel_set(kk);
        nodes_to_knn(kk, tmp_distance(kk,:) < knn_distance(kk) ) = false;
        nodes_to_knn(kk,kk) = false; % diagonal should be zero
    end
    [~,center_point]=max(sum(nodes_to_knn));
    center_set(i)=pixel_set(center_point);

    %% method 2: centroidal nodes
%     YCoordSet=[];
%     XCoordSet=[];
%     for search_i=1:length(pixel_set)
%         centerloc=pixel_set(search_i);
%         [TempYCoord, TempXCoord] = f1DTo2DCoord([Height, Width],centerloc);
%         YCoordSet=[YCoordSet;TempYCoord];
%         XCoordSet=[XCoordSet;TempXCoord];
%     end
%     centroidal_Loc_YCoord=floor(mean(YCoordSet));
%     centroidal_Loc_XCoord=floor(mean(XCoordSet));
%     centroidal_Loc=(centroidal_Loc_XCoord-1)*Height+centroidal_Loc_YCoord;
%     center_set(i)=centroidal_Loc;
end

for i = 1:sp_num
    similar_set=find(nodes_to_super(i,:)>0);
    
    %% method1
%     point_idx1=center_set(i);
%     point_idx2=center_set(similar_set);
%     nodes_to_add(point_idx1,point_idx2)=true;
%     nodes_to_add(point_idx2,point_idx1)=true;

    %% method2
    similar_set=[similar_set i];
    
    if length(similar_set)>1
        for i_set=2:length(similar_set)
            for j_set=1:i_set-1
                point_idx1=center_set(similar_set(i_set));
                point_idx2=center_set(similar_set(j_set));
                nodes_to_add(point_idx1,point_idx2)=true;
                nodes_to_add(point_idx2,point_idx1)=true;
            end
        end
        
    else
        printtex=['³öÏÖ¹ÂÁ¢³¬ÏñËØ¿é,³¬ÏñËØ¿éÐòºÅ£º',num2str(i)]
    end
end
for i=1:N
    nodes_to_add(i,i)=false;
end
nodes_to_add( nodes_to_add ~= nodes_to_add' ) = true;
