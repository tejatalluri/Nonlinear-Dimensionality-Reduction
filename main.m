tic
clear all;
close all;
clc;
%  rng default
load 'swiss_roll_data.mat'
 Y_data1=Y_data(1,1:3000);
R=Y_data1;

num_nodes = 3000; L = 100; max_seg_length = 30; ids = (1:num_nodes)';
    nodes = [ids L*R']; % nodes
    a = 25;
c = linspace(1,10,numel(nodes(:,2)));
figure
scatter(nodes(:,1),nodes(:,2),a,c,'filled')
     title('Input graph');
     figure
     N=10;
[X,Y,Z] = cylinder(nodes(:,2),N);
testsubject = surf(X,Y,Z); 
set(testsubject,'FaceAlpha',1)
set(testsubject,'FaceAlpha',0.5,'EdgeColor','red','EdgeAlpha',0.6,...
'DiffuseStrength',1,'AmbientStrength',1)
title('The high-dimensional data in R');
    i=1;
  j=1;
       while j~=3000
            w(i,j)=mod(R(i,j),R(i,j+1));
            j=j+1;
       end
       w=[w 0.325];
       w=w';
    num_segs = 0; segments = zeros(num_nodes,3);
    for i = 1:num_nodes-1 % create edges between some of the nodes
        for j = i+1:num_nodes
            d = sqrt(sum((nodes(i,2) - nodes(j,2)).^2));
            if and(d < max_seg_length,rand < 0.6)
                % add this link to the segments list
                num_segs = num_segs + 1;
                segments(num_segs,:) = [num_segs nodes(i,1) nodes(j,1)];
            end
        end
    end
    segments(num_segs+1:num_nodes,:) = [];
    omega={};
    pd=1;
    while (pd~=186)
    R=ceil(R);
    segments(num_segs+1:num_nodes,:) = [];
    start_id = abs((datasample(R,1))); 
    disp(['start id = ' num2str(start_id)]);
    finish_id = abs((R(:,end)-start_id )); disp(['finish id = ' num2str(finish_id)]);
    dijkstra
   path(cellfun(@(path) any(isnan(path)),path)) = [];
    newpath=path;
       [~,Id] = sort(cellfun(@length,newpath));
       newpath = newpath(Id);
       if size(newpath,2)>=20
           
       if (~isempty(newpath{1,20}))
  Lstar=newpath{1,20};
       
  omega{1,pd}=Lstar;
[row col]=size(R);
Lst=zeros(row,col);
for jj=1:numel(Lstar)
    Lst(1,jj)=Lstar(1,jj);
end
R=R-Lst;
       end
        end
  pd=pd+1;
     end
   
  Covernodes=cell2mat(omega);
  num_nodes1 = numel(Covernodes); L1 = 100; max_seg_length1 = 30; ids1 = (1:num_nodes1)';
    nodes1 = [ids1 L1*Covernodes']; % nodes
       a = 25;
c1= linspace(1,10,numel(nodes1(:,2)));
figure
scatter(nodes1(:,1),nodes1(:,2),a,c1,'filled')
     title('output graph');
toc 
figure
N=10;
[X,Y,Z] = cylinder(nodes1(:,2),N);
testsubject = surf(X,Y,Z); 
set(testsubject,'FaceAlpha',1)
set(testsubject,'FaceAlpha',0.5,'EdgeColor','red','EdgeAlpha',0.6,...
'DiffuseStrength',1,'AmbientStrength',1)
title('Result of the SSPC in R');

xx=nodes1(:,2);
yx=nodes1(:,1);
zx=xx.^2+yx.^2; % a function in this example but could be

tri=delaunay(xx,yx);
 tri=tri'; % reshape in 3 x ntri
 figure
 patch(xx(tri),yx(tri),zx(tri),zx(tri));
 title('unfolded manifold');
