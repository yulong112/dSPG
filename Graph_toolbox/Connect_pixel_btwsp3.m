function nodes_to_add = Connect_pixel_btwsp3(nodes_to_super,connect_representative_pixel,superpix_img, Height, Width)
N=Height*Width;
nodes_to_add = false(N);
sp_num=size(nodes_to_super,1);
center_set=zeros(sp_num,1);

% for i = 1:sp_num
%     pixel_set = find(superpix_img==i);
%     center_set(i)=pixel_set(1);
%     %% method 2: centroidal nodes
% %     YCoordSet=[];
% %     XCoordSet=[];
% %     for search_i=1:length(pixel_set)
% %         centerloc=pixel_set(search_i);
% %         [TempYCoord, TempXCoord] = f1DTo2DCoord([Height, Width],centerloc);
% %         YCoordSet=[YCoordSet;TempYCoord];
% %         XCoordSet=[XCoordSet;TempXCoord];
% %     end
% %     centroidal_Loc_YCoord=floor(mean(YCoordSet));
% %     centroidal_Loc_XCoord=floor(mean(XCoordSet));
% %     centroidal_Loc=(centroidal_Loc_XCoord-1)*Height+centroidal_Loc_YCoord;
% %     center_set(i)=centroidal_Loc;
% end

for i = 1:sp_num
    similar_set=find(nodes_to_super(i,:)>0);
    %% method1
%     point_idx1=center_set(i);
%     point_idx2=center_set(similar_set);
%     nodes_to_add(point_idx1,point_idx2)=true;
%     nodes_to_add(point_idx2,point_idx1)=true;

    %% method2
    representative_pixel1 = connect_representative_pixel(i,similar_set);
    representative_pixel2 = connect_representative_pixel(similar_set,i);
    for i_set=1:length(representative_pixel1)
        pixel1=representative_pixel1(i_set);
        pixel2=representative_pixel2(i_set);
        nodes_to_add(pixel1,pixel2)=true;
        nodes_to_add(pixel2,pixel1)=true;
    end
    
    %% method4
%     similar_set=[similar_set i];
%     if length(similar_set)>1
%         for i_set=1:length(similar_set)
%             superpixel_i=similar_set(i_set);
%             representative_pixel1 = connect_representative_pixel(superpixel_i,similar_set);
%             representative_pixel2 = connect_representative_pixel(similar_set,superpixel_i);
%             for j_set=1:length(representative_pixel1)
%                 pixel1=representative_pixel1(j_set);
%                 pixel2=representative_pixel2(j_set);
%                 if pixel1>0 && pixel2>0
%                     nodes_to_add(pixel1,pixel2)=true;
%                     nodes_to_add(pixel2,pixel1)=true;
%                 end
%             end
%         end
%     end
    
    %% method3
%     similar_set=[similar_set i];
%     
%     if length(similar_set)>1
%         for i_set=2:length(similar_set)
%             for j_set=1:i_set-1
%                 point_idx1=center_set(similar_set(i_set));
%                 point_idx2=center_set(similar_set(j_set));
%                 nodes_to_add(point_idx1,point_idx2)=true;
%                 nodes_to_add(point_idx2,point_idx1)=true;
%             end
%         end
%         
%     else
%         printtex=['³öÏÖ¹ÂÁ¢³¬ÏñËØ¿é,³¬ÏñËØ¿éĞòºÅ£º',num2str(i)]
%     end
end

for i=1:N
    nodes_to_add(i,i)=false;
end

nodes_to_add( nodes_to_add ~= nodes_to_add' ) = true;
