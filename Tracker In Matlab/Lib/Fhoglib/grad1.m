%仅计算一列的x和y梯度
% [Gx,Gy] = grad1(I+x*h+c*w*h,Gx+c*h4,Gy+c*h4,h,w,x);%仅计算一列的x和y梯度
function [Gx,Gy] = grad1(I,h4,h,w,x,c,d)
    Current = (x-1)*h+(c-1)*w*h+1;Previous = Current-h;Next = Current+h;r = 0.5;
    Gx = zeros(1,d*h4); Gy = zeros(1,d*h4);
    if x==1
        r=1;Previous = Previous+h;
    else if x==w
            r=1;Next = Next-h;
        end
        
    end
    %计算Gx
    if (h<4)||(mod(h,4)>0)
        for y=1:h
            Gx(y)=I(Next)-I(Previous);Next=Next+1;Previous=Previous+1;
        end
    else
        for y=1:h
            Gx(y)=(I(Next)-I(Previous))*r;Next=Next+1;Previous=Previous+1;
        end
    end
    %计算Gy
    Previous = Current;Next = Previous+1;
    for y=1:h
        if y==1
            Gy(y)=(I(Next)-I(Previous))*1;Next=Next+1;
        else if y==h
                Next=Next-1;Gy(y)=(I(Next)-I(Previous))*1;
        else
            Gy(y)=(I(Next)-I(Previous))*0.5;Next=Next+1;Previous=Previous+1;
            end
        end  
    end
end