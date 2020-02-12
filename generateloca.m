% function
% generate user location

function [loca,lab]=generateloca(M,radius,type)

switch type
    case 1  % square
        loca=radius*(2*rand(2,M)-1); lab=zeros(1,M);
    case 2  % circle
        loca=zeros(2,M); lab=zeros(1,M);
        i=0;
        while i<M
            x=radius*(2*rand-1); y=radius*(2*rand-1);
            if norm([x y])<radius
                i=i+1;
                loca(:,i)=[x y]';
            end
        end
    case 3  % randomly generate 4 centers and locate users around them
        center=zeros(2,4);
        flag=1;
        while flag
            x=rand; y=rand;
            if norm([x y])<=1
                center(:,1)=radius*[x;y];
                flag=0;
            end
        end
        flag=1;
        while flag
            x=rand; y=rand;
            if norm([x y])<=1
                center(:,2)=radius*[-x;y];
                flag=0;
            end
        end
        flag=1;
        while flag
            x=rand; y=rand;
            if norm([x y])<=1
                center(:,3)=radius*[-x;-y];
                flag=0;
            end
        end
        flag=1;
        while flag
            x=rand; y=rand;
            if norm([x y])<=1
                center(:,4)=radius*[x;-y];
                flag=0;
            end
        end
        lab=randi(4,1,M);  % generate a 1-by-M vector, each element is an integer between 1 and 4
        loca=zeros(2,M);
        i=0;
        while i<M
            i=i+1;
            a=radius/30*rand;
            b=2*pi*rand;
            loca(:,i)=center(:,lab(i))+[a*cos(b) a*sin(b)]';
            if norm(loca(:,i))>radius
                i=i-1;
            end
        end    
end

