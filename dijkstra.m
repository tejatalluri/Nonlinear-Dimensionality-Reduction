
node_ids = nodes(:,1);
    [num_map_pts,cols] = size(nodes);
    table = sparse(num_map_pts,2);
    shortest_distance = Inf(num_map_pts,1);
    settled = zeros(num_map_pts,1);
    path = num2cell(NaN(num_map_pts,1));
    col = 2;
    pidx = find(start_id == node_ids);
    shortest_distance(pidx) = 0;
    table(pidx,col) = 0;
    settled(pidx) = 1;
    path(pidx) = {start_id};
    if (nargin < 4) % compute shortest path for all nodes
        while_cmd = 'sum(~settled) > 0';
    else % terminate algorithm early
        while_cmd = 'settled(zz) == 0';
        zz = find(finish_id == node_ids);
    end
    while eval(while_cmd)
        % update the table
        table(:,col-1) = table(:,col);
        table(pidx,col) = 0;
        % find neighboring nodes in the segments list
        neighbor_ids = [segments(node_ids(pidx) == segments(:,2),3);
            segments(node_ids(pidx) == segments(:,3),2)];
        % calculate the distances to the neighboring nodes and keep track of the paths
        for k = 1:length(neighbor_ids)
            cidx = find(neighbor_ids(k) == node_ids);
            if ~settled(cidx)
                d = sqrt(sum((nodes(pidx,2:cols) - nodes(cidx,2:cols)).^2));
                if (table(cidx,col-1) == 0) || ...
                        (table(cidx,col-1) > (table(pidx,col-1) + d))
                    table(cidx,col) = table(pidx,col-1) + d;
                    tmp_path = path(pidx);
                    path(cidx) = {[tmp_path{1} neighbor_ids(k)]};
                else
                    table(cidx,col) = table(cidx,col-1);
                end
            end
        end
        % find the minimum non-zero value in the table and save it
        nidx = find(table(:,col));
        ndx = find(table(nidx,col) == min(table(nidx,col)));
        if isempty(ndx)
            break
        else
            pidx = nidx(ndx(1));
            shortest_distance(pidx) = table(pidx,col);
            settled(pidx) = 1;
        end
    end
    if (nargin < 4) % return the distance and path arrays for all of the nodes
        dist = shortest_distance';
        path = path';
    else % return the distance and path for the ending node
        dist = shortest_distance(zz);
        path = path(zz);
%         path = path{1};
    end
