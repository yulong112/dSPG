
function superpix_img = superpixel_cut(HyperCube, numSuperpixels)
%% ��ԭͼ�������ҪPCA��ȡ���ɳɷ�������
[nl ns nb] = size(HyperCube);
input_pca = reshape(HyperCube, nl*ns, nb);

select_k=20;
% [PC1] = pca_row(input_pca,select_k);
[PC1] = fPCA_FinalVersionV2(HyperCube, select_k, 0, 0, 0);
img2 = reshape(PC1,nl,ns,select_k); % 145 X 145 X 20
%% �ú���������������PCA��ֱ���ú������������
% img2 = HyperCube; % 145 X 145 X 12

[nl ns nb] = size(img2);
x = img2;  % 610*340*20
x = reshape(x, nl*ns, nb);
x = x';
HIM = img2;
clear img2 PC;


%�����ݰ�����ת��һά������ʽ
input_img = zeros(1, nl * ns * nb);

startpos = 1;
for i = 1 : nl % ��
    for j = 1 : ns % ��
        input_img(startpos : startpos + nb - 1) = HIM(i, j, :); % ����
        startpos = startpos + nb;
    end
end


%% ���зָ����ָ���
if ~exist('numSuperpixels','var')
    numSuperpixels = 2000;  %ϣ���ָ������
end
compactness = 0.1; %compactness2 = 1-compactness, �����б�����Ϊ compactness*dxy+compactness2*dspectral
dist_type = 2; % 1:ŷ�Ͼ��룻2��SAD; 3:SID; 4:SAD-SID
seg_all = 1; % 1: ȫ����Ԫ��ɾ��࣬ 2������δ�������Ԫ
% labels:���������ȵķ�ʽ���ÿ����Ԫ�ķ���
% numlabels: ʵ�ʷ�����Ŀ
% seedx: ���ӵ��к������� seedy�����ӵ��к�����
[labels, numlabels, seedx, seedy] = hybridseg(input_img, nl, ns, nb, numSuperpixels, compactness, dist_type, seg_all);%numlabels is the same as number of superpixels
% note: "labels" is not in matlab format
clear input_img;

l = reshape(labels,ns,nl);
superpix_img = reshape(l',nl*ns,1);
% superpix_img=l';