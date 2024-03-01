dataset_size = 5*no_inputs;
X= linspace(0,1,dataset_size);
X=X';
B=flipud(X);
X(:,2)=B;
X(:,1)=0;        % adjustment of bounds
Y = zeros(dataset_size,no_outputs);     % corresponding matrix of response data
for k = 1:size(X,1)
    X(k,:) = X(k,:).*(ub-lb)+lb;        % adjustment of bounds
    Y(k,:) = f(X(k,:));                 % calculation of response data
end
for k = 1:size(X,1)
    if Y(k,2)==inf
     ub(2)= X(k,2);
    else
       
    end
end
