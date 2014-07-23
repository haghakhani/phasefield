clc; close all; clear all;

directory_name = '/dinesh1/data/users/haghakha/TITAN2D/git phase field/titan-2.0.2/bin/colimafine/2proc';
files = dir(fullfile(directory_name, 'xdmf*.h5'));
format short;

% fileIndex = find(~[files.isdir]); Just select the files and delete the rest

name=zeros(length(files),2);

for i=1:length(files)
    filename = files(i).name;
    name(i,:) = sscanf(filename, ['xdmf' '%02d' '%08d' '.h5'])';
end

numProc=max(name(:,1));%+1;
numIter=length(name)/(numProc+1);

for timeIter=1:numIter
    
    for proc=1:numProc+1;
        timeIter
        
        filename=sprintf('xdmf%02d%08d.h5',(proc-1),name(timeIter,2));
        
        connec = hdf5read(filename,'/Mesh/Connections');
        points = hdf5read(filename,'/Mesh/Points');
        pileh = hdf5read(filename,'/Properties/PHI');
                
        nonZeroPile=find(pileh(1,:)>0.);%scale and norm the colormap based on non zero pile height

        %===========================================================
        %creating data required for the pileheight record
        %===========================================================
        nonZeroPile=find(pileh(1,:)>0.);%scale and norm the colormap based on non zero pile height
        if (proc==1)
            figure();
        end
        
        cmap=colormap;
        p = patch('Faces',(connec+1)','Vertices',points');
        minh=min(pileh(1,nonZeroPile));
        maxh=max(pileh(1,nonZeroPile));
        normh=round((pileh-minh)/(maxh-minh)*length(cmap));
        set(gca,'CLim',[0 length(cmap)]);
        set(p,'FaceColor','flat',...
            'FaceVertexCData',normh',...
            'CDataMapping','direct',...
            'edgecolor','interp');
        
        axis image;
        
        colorbar('YTickLabel',...
            {'Freezing','Cold','Cool','Neutral',...
            'Warm','Hot','Burning','Nuclear'});
        shading flat;
%         M(timeIter) = getframe(gcf);
        
%         pause(10);
    end
end