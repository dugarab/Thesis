function [A1] = sparsiofy(Aactual, technique, sparsity)
    N = length(Aactual);
    if ~exist('technique', 'var')
        technique = 'low';
    end
    
    if ~exist('sparsity', 'var')
        sparsity = (4*log2(N))/N;
    end
    temp = abs(Aactual);
    X = sort(temp')';
    
    
    if strcmpi(technique,'low')
        disp('Removing lower order terms');
        sort1 = X(:,N-ceil(N*sparsity));
        sort1 = sort1*ones(1,N);
        f1 = 1.0*(temp>sort1)  ;
        A1 = (f1.*Aactual + diag(diag(Aactual)));
        
    elseif strcmpi(technique, 'mid')
        sort1 = X(:,ceil(N*sparsity/2));
        sort2 = X(:,N-ceil(N*sparsity/2));
        sort1 = sort1*ones(1,N);
        sort2 = sort2*ones(1,N);
        f1 = 1.0*(temp<sort1) + 1.0*(temp>sort2) ;
        A1 = (f1.*Aactual + diag(diag(Aactual)));
        
    elseif strcmpi(technique, 'high')
        sort1 = X(:,floor(N*sparsity));
        sort1 = sort1*ones(1,N);
        f1 = 1.0*(temp<sort1)  ;
        A1 = (f1.*Aactual + diag(diag(Aactual)));
        
    else
        disp('Incorrect technique chosen');
        A1 = Aactual;
        
    end
        
end
