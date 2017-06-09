%compute block normalization values 计算块归一化值
%N = hogNormMatrix(R2, nOrients, hb, wb, binSize)
%compute 2x2 block normalization values (padded by 1 pixel)HOG帮手：计算2x2块归一化值（用1像素）
function N = hogNormMatrix(R2, nOrients, hb, wb, binSize)
    hb1=hb+1;wb1=wb+1;ep=1e-4/4/binSize/binSize/binSize/binSize;%precise backward equality精确的反向等式
    N=zeros(1,hb1*wb1);
    current_N1=hb1;
   for o=1: nOrients
       for x=1:wb
           for y=1:hb
               N(current_N1+(x-1)*hb1+y)=N(current_N1+(x-1)*hb1+y)+R2((o-1)*wb*hb + (x-1)*hb + y)*R2((o-1)*wb*hb + (x-1)*hb + y);
           end
       end
   end
   for x=1:wb-1
       for y=1:hb-1
           n=current_N1+(x-1)*hb1+y;N(n)=1/(sqrt(N(n)+N(n+1)+N(n+hb1)+N(n+hb1+1)+ep));
       end
   end
   x=1;dx=1;dy=1;y=1;N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);
   x=1;dx=1;dy=0;for y=1:hb1,N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);end
   x=1;dx=1;dy=-1;y=hb1;N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);
   x=wb1;dx=-1;dy=1;y=1;N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);
   x=wb1;dx=-1;dy=0;for y=1:hb1,N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);end
   x=wb1;dx=-1;dy=-1;y=hb1;N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);
   dx=0;dy=1;y=1;for x=1:wb1,N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);end
   dx=0;dy=-1;y=hb;for x=1:wb1,N((x-1)*hb1+y) = N(((x-1) + dx)*hb1 + y + dy);
end