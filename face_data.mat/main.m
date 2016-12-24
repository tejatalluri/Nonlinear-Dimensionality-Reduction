tic
clear all;
close all;
clc;
%  rng default
load 'face_data.mat'
 images=imresize(lights,[1 3000]);
R=images;

num_nodes = 3000; L = 100; max_seg_length = 30; ids = (1:num_nodes)';
    nodes = [ids L*R']; % nodes
    %     a = 25;
% c = linspace(1,10,numel(nodes(:,2)));
% figure
% scatter(nodes(:,1),nodes(:,2),a,c,'filled')
%      title('Input graph');
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
       if (~isempty(newpath))
  Lstar=newpath{1,20};
  omega{1,pd}=Lstar;
[row col]=size(R);
Lst=zeros(row,col);
for jj=1:numel(Lstar)
    Lst(1,jj)=Lstar(1,jj);
end
R=R-Lst;
        end
  pd=pd+1;
     end
   
  Covernodes=cell2mat(omega);
  num_nodes1 = numel(Covernodes); L1 = 100; max_seg_length1 = 30; ids1 = (1:num_nodes1)';
    nodes1 = [ids1 L1*Covernodes']; % nodes
      c1 = linspace(1,10,length( nodes1(:,2)));
      figure
scatter( nodes1(:,2),nodes1(:,1),[],c1)

     title('output graph');
toc 