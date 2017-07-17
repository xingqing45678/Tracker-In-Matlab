%Qinwei----2017.3.27
%在每个位置计算梯度的大小和方向
% USAGE
%  [M,O] = gradientMag( I, [channel], [normRad], [normConst], [full] )
%
% INPUTS
%  I          - [hxwxk] input k channel single image
%  channel    - [0] if>0 color channel to use for gradient computation
%  normRad    - [0] normalization radius (no normalization if 0)
%  normConst  - [.005] normalization constant
%  full       - [0] if true compute angles in [0,2*pi) else in [0,pi)
%
% OUTPUTS
%  M          - [hxw] gradient magnitude at each location
%  O          - [hxw] approximate gradient orientation modulo PI

function [M,O] =  gradientMag(I, channel, normRad, normConst, full)
    
    if(nargin<1 || isempty(I)), M=single([]); O=M; return; end
    if(nargin<2 || isempty(channel)), channel=0; end
    if(nargin<3 || isempty(normRad)), normRad=0; end
    if(nargin<4 || isempty(normConst)), normConst=.005; end
    if(nargin<5 || isempty(full)), full=1; end
    
    %% 查表 build lookup table acosTable[] s.t. acosTable[x*n]~=acos(x) for x in [-1,1]
    n = 10000;b = 10;
    acosTable = zeros(1,n*2+b*2);
    for i=-n-b:-n-1 
        acosTable(n+b+i+1) = pi;
    end
    for i=-n:n-1     
        acosTable(n+b+i+1)=single(acos(i/single(n)));
    end
    for i=n:n+b-1    
        acosTable(n+b+i+1)=0;
    end
    for i=-n-b:n/10-1 
        if acosTable(n+b+i+1) > pi-(1e-6)
            acosTable(n+b+i+1)=pi-(1e-6);
        end 
    end
    %% 
    
    [h,w] = size(I); M = zeros(h,w); O = zeros(h,w);acMult = 10000;
    d = 1;%彩色图还是灰度图
    if mod(h,4)==0
        h4=h;
    else
        h4=h-mod(h,4)+4;
    end
    M2 = zeros(1,d*h4); Gx = zeros(1,d*h4); Gy = zeros(1,d*h4);
    for x=1:w
       for c=1:d
           [Gx,Gy] = grad1(I,h4,h,w,x,c,d);%仅计算一列的x和y梯度
           for y=1:h4
              y1 = h4*(c-1)+y;
              M2(y1) = (Gx(y1)*Gx(y1)+Gy(y1)*Gy(y1));%求Gx和Gy的平方和 M2   
              if c==1
                  continue;
              end
              if M2(y1) > M2(y)%取较大的值
                  M2(y) = M2(y1);
              end
              if Gx(y1) > Gx(y)
                  Gx(y) = Gx(y1);
              end
              if Gx(y1) > Gx(y)
                  Gx(y) = Gx(y1);
              end   
           end
       end
%        计算梯度幅值M和范数Gx
       for y=1:h4
          m = 1/sqrt(M2(y));
          if m>1e10,m = 1e10;end
          M2(y) = 1/m;%此处M2内容变成了sqrt（Gx2+Gy2）
          Gx(y) = (Gx(y)*m)*acMult;
          %
          if Gy(y)<0
              Gx(y)=Gx(y)*(-1);
          end
%         %按位与、按位异或作用未知
%           Gx(y) = bitxor(Gx(y),bitand(Gy(y),0.0));
       end
          M(:,x) = M2(1:h);%将M2的值复制到 M 中
          %计算和存储梯度方向（O）通过查找表
          %计算梯度方向时用到反正切，但是注意要将角度值映射到0~2pi
          for y=1:h
             try
                O((x-1)*h+y) = acosTable(round(Gx(y))+n+b);%+10011);%原方式查表法
             catch
                 disp('error');
             end
%                if Gx(y)==0
%                    if Gy(y)>=0,pv=pi/2;else pv=pi*3/2;end
%                else
%                    if Gy(y)==0
%                        if Gx>0,pv=0;else pv=pi;end
%                    else if (Gy(y)&&Gx(y))>0
%                        if Gy(y)>0,pv=atan(Gy(y)/Gx(y));else pv=atan(Gy(y)/Gx(y))+pi;end
%                    else if (Gy(y)&&Gx(y))<0
%                        if Gy(y)>0,pv=atan(Gy(y)/Gx(y))+pi;else pv=atan(Gy(y)/Gx(y))+2*pi;end
%                        end
%                        end
%                    end
%                end                                            
%                % pv = atan2(Gy(y),Gx(y));%atan2范围从-pi~pi，atan范围为-pi/2~pi/2;
%              O((x-1)*h+y) = pv;%O存储的是归一化后的梯度方向
          end
          if full%如果Gy小于0则将O相应位置方向值+180度，否则不作变换
              for y=1:h
                  if Gy(y)<0
                      O(y+(x-1)*h)=O(y+(x-1)*h)+pi;
                  end
              end
          end
    end
end