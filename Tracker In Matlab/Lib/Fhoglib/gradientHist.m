%% 计算梯度直方图
function R = gradientHist(M,O,h,w,binSize,nOrients,softBin,full)
%gradientHist(M, O, h, w, binSize, nOrients * 2, softBin, true);
    hb = floor(h/binSize);wb = floor(w/binSize);nb = hb*wb;w0 = wb*binSize;h0 = hb*binSize;
    sInv = 1/binSize;sInv2 = 1/binSize/binSize;R = zeros(1,wb*hb*nOrients);
    for x=1:w0
        [M0,M1,O0,O1] = gradQuantize(M,O,x,h,nb,h0,sInv2,nOrients,full,softBin>=0);%量化O和M
        if softBin<0&&mod(softBin,2)==0  %表明softBin为小于0的偶数，不采用插值
            current=floor((x-1)/binSize)*hb+1;
            if binSize==1
                for y=1:h0,R(current+O0(y))=R(current+O0(y))+M0(y);current=current+1;end
            else if binSize==2
                    for y=1:2:h0,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O0(y+1))=R(current+O0(y+1))+M0(y+1);current=current+1;end
                else if binSize==3
                       for y=1:3:h0,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O0(y+1))=R(current+O0(y+1))+M0(y+1);R(current+O0(y+2))=R(current+O0(y+2))+M0(y+2);current=current+1;end
                    else if binSize==4
                            for y=1:3:h0,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O0(y+1))=R(current+O0(y+1))+M0(y+1);R(current+O0(y+2))=R(current+O0(y+2))+M0(y+2);R(current+O0(y+3))=R(current+O0(y+3))+M0(y+3);current=current+1;end
                        else
                            for y=1:h0
                                for y1=1:binSize,R(current+O0(y))=R(current+O0(y))+M0(y);end
                            end
                        end
                    end
                end
            end
        else if mod(softBin,2)==0 ||binSize==1  %表明softBin为>=0的偶数，采用线性插值
             current=floor((x-1)/binSize)*hb+1;
             if binSize==1
                    for y=1:h0,R(current+O0(y))=R(current+O0(y))+M0(y);
                               R(current+O1(y))=R(current+O1(y))+M1(y);current=current+1;
                    end
                else if binSize==2
                        for y=1:2:h0,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O0(y+1))=R(current+O0(y+1))+M0(y+1);
                                     R(current+O1(y))=R(current+O1(y))+M1(y);R(current+O1(y+1))=R(current+O1(y+1))+M1(y+1);current=current+1;
                        end
                    else if binSize==3
                           for y=1:3:h0,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O0(y+1))=R(current+O0(y+1))+M0(y+1);R(current+O0(y+2))=R(current+O0(y+2))+M0(y+2);
                                        R(current+O1(y))=R(current+O1(y))+M1(y);R(current+O1(y+1))=R(current+O1(y+1))+M1(y+1);R(current+O1(y+2))=R(current+O1(y+2))+M1(y+2);current=current+1;
                           end
                        else if binSize==4
                                for y=1:3:h0,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O0(y+1))=R(current+O0(y+1))+M0(y+1);R(current+O0(y+2))=R(current+O0(y+2))+M0(y+2);R(current+O0(y+3))=R(current+O0(y+3))+M0(y+3);
                                             R(current+O1(y))=R(current+O1(y))+M1(y);R(current+O1(y+1))=R(current+O1(y+1))+M1(y+1);R(current+O1(y+2))=R(current+O1(y+2))+M1(y+2);R(current+O1(y+3))=R(current+O1(y+3))+M1(y+3);current=current+1;
                                end
                            else
                                for y=1:h0
                                    for y1=1:binSize,R(current+O0(y))=R(current+O0(y))+M0(y);R(current+O1(y))=R(current+O1(y))+M1(y);end
                                end
                            end
                        end
                    end
             end
            else   %表明softBin为奇数采用三线性插值
                if x==1 ,init=-0.5*sInv;xb=init;end%init=-0.5*sInv;%0.5*sInv-0.5
                if xb>=0,hasLf=1;else hasLf=0;end
                if hasLf,xb0=floor(xb);else xb0=-1;end
                if xb0<wb-1,hasRt=1;else hasRt=0;end
                xd=xb-xb0;xb=xb+sInv;yb=init;
                %宏定义
                expression='yd=yb-yb0;yb=yb+sInv;base=xb0*hb+yb0+1;xyd=xd*yd;ms(1)=1-xd-yd+xyd;ms(2)=yd-xyd;ms(3)=xd-xyd;ms(4)=xyd;';
                %leading rows,no top bin
                for y=1:floor(binSize/2)
                    yb0=-1;eval(expression);current=base;
                    if hasLf,R(current+O0(y)+1)=R(current+O0(y)+1)+ms(2)*M0(y);R(current+O1(y)+1)=R(current+O1(y)+1)+ms(2)*M1(y);end
                    if hasRt,R(current+O0(y)+hb+1)=R(current+O0(y)+hb+1)+ms(4)*M0(y);R(current+O1(y)+hb+1)=R(current+O1(y)+hb+1)+ms(4)*M1(y);end
                end
                %main rows,has top and bottom bins
                if softBin<0
                    while(1)
                        yb0=floor(yb);
                        if yb0>=hb-1,break;end
                        eval(expression);
                        if hasLf %源程序中sse加速是四位相乘，结果存入当前位置，未使用sse加速的考虑了后三个点情况
%                             current=base+O0(y);R(current)=R(current)+0*M0(y);R(current+1)=R(current+1)+0*M0(y);
%                             R(current+2)=R(current+2)+ms(2)*M0(y);R(current+3)=R(current+3)+ms(1)*M0(y);
                            current=base+O0(y);R(current)=R(current)+0*M0(y)+0*M0(y)+ms(2)*M0(y)+ms(1)*M0(y);
                        end
                        if hasRt
%                             current=base+O0(y)+hb;R(current)=R(current)+0*M0(y);R(current+1)=R(current+1)+0*M0(y);
%                             R(current+2)=R(current+2)+ms(4)*M0(y);R(current+3)=R(current+3)+ms(3)*M0(y);
                            current=base+O0(y)+hb;R(current)=R(current)+0*M0(y)+0*M0(y)+ms(4)*M0(y)+ms(3)*M0(y);
                        end
                        y=y+1;
                    end
                else
                    while(1)
                        yb=floor(yb);
                        if yb0>=hb-1,break;end
                        eval(expression);
                        if hasLf %同上分析，如果sse加速则当前位置值与后三个位置的无关，否则就要改为类似上述的表示
                            current=base+O0(y);R(current)=R(current)+0*M0(y);R(current+1)=R(current+1)+0*M0(y);
                            R(current+2)=R(current+2)+ms(2)*M0(y);R(current+3)=R(current+3)+ms(1)*M0(y);
                            current=base+O1(y);R(current)=R(current)+0*M1(y);R(current+1)=R(current+1)+0*M1(y);
                            R(current+2)=R(current+2)+ms(2)*M1(y);R(current+3)=R(current+3)+ms(1)*M1(y);
                        end
                        if hsaRt
                            current=base+O0(y)+hb;R(current)=R(current)+0*M0(y);R(current+1)=R(current+1)+0*M0(y);
                            R(current+2)=R(current+2)+ms(4)*M0(y);R(current+3)=R(current+3)+ms(3)*M0(y);
                            current=base+O1(y);R(current)=R(current)+0*M1(y);R(current+1)=R(current+1)+0*M1(y);
                            R(current+2)=R(current+2)+ms(4)*M1(y);R(current+3)=R(current+3)+ms(3)*M1(y);
                        end
                        y=y+1;
                    end
                end
                %final rows,no bottom bin
                for z=y:h0 %z相当于从y接着递增直到h0
                    yb0=floor(yb);eval(expression);
                    current=base;
                    if hasLf
                        R(current+O0(z))=R(current+O0(z))+ms(1)*M0(z);R(current+O1(z))=R(current+O1(z))+ms(1)*M1(z);
                    end
                    if hasRt
                        R(current+O0(z)+hb)=R(current+O0(z)+hb)+ms(3)*M0(z);R(current+O1(z)+hb)=R(current+O1(z)+hb)+ms(3)*M1(z);
                    end
                end
            end   %三线性插值的末尾
        end
    end%主循环x的末尾
    %normalize boundary bins which only get 7/8 of weight of interior bins
    %归一化边界bin，只有内部bin的7/8权重
    if mod(softBin,2)~=0
        for o=1:nOrients
            x=1;for y=1:hb,R((o-1)*nb+(x-1)*hb+y)=R((o-1)*nb+(x-1)*hb+y)*8.0/7.0;end
            y=1;for x=1:wb,R((o-1)*nb+(x-1)*hb+y)=R((o-1)*nb+(x-1)*hb+y)*8.0/7.0;end
            x=wb;for y=1:hb,R((o-1)*nb+(x-1)*hb+y)=R((o-1)*nb+(x-1)*hb+y)*8.0/7.0;end
            y=hb;for x=1:wb,R((o-1)*nb+(x-1)*hb+y)=R((o-1)*nb+(x-1)*hb+y)*8.0/7.0;end
        end       
    end
end