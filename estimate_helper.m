network = 'kronecker-core-periphery-n1024-h10-r0_01-0_25-network.txt';
cascades = 'kronecker-core-periphery-n1024-h10-r0_01-0_25-1000-cascades.txt';
horizon = 10;
type_diffusion = 'exp';
num_nodes = 1024;
A = create_adj_matrix(network, num_nodes);
C = create_cascades(cascades, num_nodes);


num_cascades = zeros(1,num_nodes);
A_potential = sparse(zeros(size(A)));
A_bad = sparse(zeros(size(A)));
A_hat = sparse(zeros(size(A)));
total_obj = 0;

for c=1:size(C, 1),
    idx = find(C(c,:)~=-1); % used nodes
    [val, ord] = sort(C(c, idx));
    
    for i=2:length(val),
        num_cascades(idx(ord(i))) = num_cascades(idx(ord(i))) + 1;
        for j=1:i-1,
            if (strcmp(type_diffusion, 'exp'))
                A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+val(i)-val(j);
            elseif (strcmp(type_diffusion, 'pl') && (val(i)-val(j)) > 1)
                A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+log(val(i)-val(j));
            elseif (strcmp(type_diffusion, 'rayleigh'))
                A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+0.5*(val(i)-val(j))^2;
            end
        end
    end
    
    for j=1:num_nodes,
        if isempty(find(idx==j))
            for i=1:length(val),
                if (strcmp(type_diffusion, 'exp'))
                    A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + (horizon-val(i));
                elseif (strcmp(type_diffusion, 'pl') && (horizon-val(i)) > 1)
                    A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + log(horizon-val(i));
                elseif (strcmp(type_diffusion, 'rayleigh'))
                    A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + 0.5*(horizon-val(i))^2;
                end
            end
        end
    end
end

save('saveAp.mat','A_potential');
save('saveAb.mat','A_bad');
save('savecas.mat','num_cascades');