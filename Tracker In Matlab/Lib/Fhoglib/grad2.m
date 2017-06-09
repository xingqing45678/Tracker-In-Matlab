%在每个位置计算x和y梯度
function [Gx,Gy]=grad2(I, Gx, Gy, h, w, d)
    a = w*h;
    for c=1:d
        for x=1:w
            griad1(I,h4,h,w,x,c,d);
        end
    end
end