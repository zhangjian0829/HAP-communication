% function
% divide users into G groups and calculate WW

function [X,WW]=usergrouping_dist(alpha,rho,A,G)
[N,M]=size(A);

F=A*diag(sqrt(alpha).*rho);

% initialize WW
ind=zeros(1,G);
ind(1)=randi(M);
for g=2:G
    met=zeros(g-1,M);
    for k=1:M
        for i=1:g-1
%             met(i,k)=rho(k)^2*norm(A(:,k)'*A(:,ind(i)))^2-rho(k)^2*N;
            met(i,k)=norm(F(:,k)*F(:,k)'-F(:,ind(i))*F(:,ind(i))','fro')^2;
        end
    end
    if size(met,1)>1
%         met=max(met);
        met=min(met);
    end
%     ind(g)=find(met==min(met));
    ind(g)=find(met==max(met));
    
end

WW=A(:,ind);

eta=zeros(G,M);
flag=1;
while flag
    X=zeros(G,M);
    for k=1:M
        for g=1:G
%             eta(g,k)=rho(k)^2*norm(A(:,k)'*WW(:,g))^2-rho(k)^2*N;
            eta(g,k)=norm(F(:,k)*F(:,k)'-WW(:,g)*WW(:,g)','fro')^2;
        end
%         [val,lab]=max(eta(:,k));
        [val,lab]=min(eta(:,k));
        X(lab,k)=1;
    end
    
    % update WW
%     WW0=centerupdate(WW,rho,A,X);
    WW0=centerupdate_dist(F,X);
    
    % converge?
    if norm(WW-WW0,'fro')<1e-1
        flag=0;
    end
    
    WW=WW0;
end







