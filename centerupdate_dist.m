% function
% center update of K-means

function W=centerupdate_dist(F,X)
N=size(F,1);
G=size(X,1);
num=sum(X,2);

W=zeros(N,G);

for g=1:G
    ind=find(X(g,:));
    A=zeros(N);
    for i=1:num(g)
        A=A+F(:,ind(i))*F(:,ind(i))'./num(g);
    end
    [V,D]=eig(A);
    W(:,g)=V(:,end).*(sqrt(N)/norm(V(:,end)));
end









