clc; close all; clear all;

directory_name = '/dinesh1/data/users/haghakha/TITAN2D/git phase field/titan-2.0.2/bin/colimafine';
files = dir(fullfile(directory_name, 'xdmf*.h5'));
format short;

% fileIndex = find(~[files.isdir]); Just select the files and delete the rest

name=zeros(length(files),2);

for i=1:length(files)
    filename = files(i).name;
    name(i,:) = sscanf(filename, ['xdmf' '%02d' '%08d' '.h5'])';
end

numProc=max(name(:,1));
numIter=length(name);

for timeIter=numIter:numIter
    for proc=1:numProc+1
        timeIter
        
        filename=sprintf('xdmf%02d%08d.h5',(proc-1),name(timeIter,2));
        
        connec = hdf5read(filename,'/Mesh/Connections');
        points = hdf5read(filename,'/Mesh/Points');
        pileh = hdf5read(filename,'/Properties/PILE_HEIGHT');
        %===========================================================
        %plotting the result
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
        
        %         pause(1);
        %         close all;
        %===========================================================
        %creating data required for the pileheight record
        %===========================================================
        rx=zeros(length(nonZeroPile),1);
        lx=rx;
        uy=rx;
        ly=rx;
        hsquare=rx;
        for k=1:length(nonZeroPile)
            %             k
            vertices=points(:,connec(:,nonZeroPile(k))+1);
            rx(k)=max(vertices(1,:));
            lx(k)=min(vertices(1,:));
            uy(k)=max(vertices(2,:));
            ly(k)=min(vertices(2,:));
            hsquare(k,1)=pileh(nonZeroPile(k));%reads the pile height for curent region
        end
        if (proc==1)
            nZH=hsquare;
            nZX(:,1)=lx;
            nZX(:,2)=rx;
            nZY(:,1)=ly;
            nZY(:,2)=uy;
        else
            
            newel=length(hsquare);
            currelem=length(nZH);
            
            nZH(currelem+1:currelem+newel,1)=hsquare;
            nZX(currelem+1:currelem+newel,1)=lx;
            nZX(currelem+1:currelem+newel,2)=rx;
            nZY(currelem+1:currelem+newel,1)=ly;
            nZY(currelem+1:currelem+newel,2)=uy;
        end;
        
        
    end
    
end

% for thresh=0:0%.025:.1
thresh=0.0;
tobeploted=find(nZH>=thresh);
numSquare=length(tobeploted);
connecSquare=zeros(numSquare,4);
verSquare=zeros(numSquare*4,2);


for square=1:numSquare
    ind=tobeploted(square);
    verNum=4*square;
    verSquare(verNum-3,:)=[nZX(ind,1),nZY(ind,1)];
    verSquare(verNum-2,:)=[nZX(ind,1),nZY(ind,2)];
    verSquare(verNum-1,:)=[nZX(ind,2),nZY(ind,2)];
    verSquare(verNum  ,:)=[nZX(ind,2),nZY(ind,1)];
    
    connecSquare(square,:)=[verNum-3,verNum-2,verNum-1,verNum];
    
end

fig=figure();
p = patch('Faces',connecSquare,'Vertices',verSquare);
clear maxh normh;
cmap=colormap;
% nonZeroPile=find(pileh(1,:)>0.);%scale and norm the colormap based on non zero pile height
minh=min(nZH(tobeploted));
maxh=max(nZH(tobeploted));
normh=round((nZH(tobeploted)-minh)/(maxh-minh)*length(cmap));
set(gca,'CLim',[0 length(cmap)]);
set(p,'FaceColor','flat',...
    'FaceVertexCData',normh,...
    'CDataMapping','direct',...
    'edgecolor','interp');

axis image;

colorbar('YTickLabel',...
    {'Freezing','Cold','Cool','Neutral',...
    'Warm','Hot','Burning','Nuclear'});
shading flat;
picname=sprintf('thresh_%f_iter_%d',thresh,numIter);
    print(fig,'-dpng',picname);

% end
