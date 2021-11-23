function tract_clean = TractSeg(LabelMap,EndROI,tract)
% tract_clean = TractSeg(LabelMap,EndROI,tract)
% collect segments of tractography between two endpoints
% LabelMap: a Nifti file for a label map of brain
% EndROI=[ROI1 ROI2]: two endpoints of the fiber segment
% tract: tract structure of input tractography, tract.fibers: ras
% cooridinates of tractography, 
% tract_clean: cleaned tract structure

Nii=MRIread(LabelMap);
Label=Nii.vol;
[nx,ny,nz]=size(Nii.vol);
Vox2RAS=Nii.vox2ras;


M=length(tract.fibers);
IND=cell(1,M);
fl=zeros(1,M);% fiber length

fiber_keep = zeros(1,M);

for m=1:M
    lof = size(tract.fibers{m},2);
    ind_m = zeros(1,lof);
    IJK=(Vox2RAS)\[tract.fibers{m}; ones(1,size(tract.fibers{m},2))];
    IJK=ceil(IJK(1:3,:));
    IJK(IJK<=0)=1;
    IJK(3,IJK(3,:)>nz)=nz;
    ind=sub2ind([nx,ny,nz],IJK(2,:),IJK(1,:),IJK(3,:));
    fiber_label=Label(ind);
    if(sum(ismember(EndROI,fiber_label))==2)% if the fiber goes through both ROIs
        fiber_keep(m) = 1;
        id_1 = find(fiber_label==EndROI(1));
        id_2 = find(fiber_label==EndROI(2));
        if(median(id_1)<median(id_2))%EndROI(1)is the starting point
            points_id = [min(id_1) max(id_2)]; 
        else
            points_id = [min(id_2) max(id_1)]; 
        end
        tract.fibers{m}=tract.fibers{m}(:,points_id);
        ind_m(points_id) = 1;
        
        fl(m)=size(tract.fibers{m},2);
        IND{m}=ind_m;
    end
    ind = cell2mat(IND);   
    FN = fieldnames(tract); 
    for l = 1:length(FN)
        if(~strcmp(FN{l},'fibers'))
            eval(['tract.' FN{l} '=tract.' FN{l} '(ind>0);']);
        end
    end
end

tract_clean = tract;