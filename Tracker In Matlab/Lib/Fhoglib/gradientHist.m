%% 计算梯度直方图
function R = gradientHist(M,O,h,w,binSize,nOrients,softBin,full)
%gradientHist(M, O, h, w, binSize, nOrients * 2, softBin, true);
    hb = floor(h/binSize);wb = floor(w/binSize);nb = hb*wb;w0 = wb*binSize;h0 = hb*binSize;
    sInv = 1/binSize;sInv2 = 1/binSize/binSize;R = zeros(1,wb*hb*nOrients*3/2);
    for x=1:w0
        [M0,M1,O0,O1] = gradQuantize(M,O,x,h,nb,h0,sInv2,nOrients,full,softBin>=0);%量化O和M
        if softBin==-1  %%表明softBin为奇数采用三线性插值           
            %出错，原因是c语言版本下标从0开始到216，而matlab下标从1开始，故会超界。
                if x==1 ,init=0.5*sInv-0.5;xb=init;end%此处有改动：init=-0.5*sInv;%源：0.5*sInv-0.5
                if xb>=0,hasLf=1;else hasLf=0;end
                if hasLf,xb0=floor(xb);else xb0=-1;end
                if xb0<wb-1,hasRt=1;else hasRt=0;end
                xd=xb-xb0;xb=xb+sInv;yb=init;base = 0;
                %宏定义
                expression='yd=yb-yb0;yb=yb+sInv;base=xb0*hb+yb0;xyd=xd*yd;ms(1)=1-xd-yd+xyd;ms(2)=yd-xyd;ms(3)=xd-xyd;ms(4)=xyd;';
                
                %leading rows,no top bin
                for y=1:floor(binSize/2)
                    yb0=-1;
                    eval(expression);
                    current=base+1;
                    if hasLf
                        R(current+O0(y)+1)=R(current+O0(y)+1)+ms(2)*M0(y);
                        R(current+O1(y)+1)=R(current+O1(y)+1)+ms(2)*M1(y);
                    end
                    if hasRt
                        R(current+O0(y)+hb+1)=R(current+O0(y)+hb+1)+ms(4)*M0(y);
                        R(current+O1(y)+hb+1)=R(current+O1(y)+hb+1)+ms(4)*M1(y);
                    end
                end
                y=y+1;
                %main rows,has top and bottom bins
                    while(1)
                        yb0=floor(yb);
                        if yb0>=hb-1,break;end
                        eval(expression);
                        try
                            if hasLf %源程序中sse加速是四位相乘，结果存入当前位置，未使用sse加速的考虑了后三个点情况
                                current=base+O0(y)+1;
%                                 if current>=213
%                                     current=213;
%                                 end
                                R(current)=R(current)+ms(1)*M0(y);
                                R(current+1)=R(current+1)+ms(2)*M0(y);
                                R(current+2)=R(current+2)+0*M0(y);
                                R(current+3)=R(current+3)+0*M0(y);%out of range
    %                             current=base+O0(y);R(current)=R(current)+0*M0(y)+0*M0(y)+ms(2)*M0(y)+ms(1)*M0(y);%by qw
                            end
                            if hasRt
                                current=base+O0(y)+hb+1;
%                                 if current>213
%                                     current=213;
%                                 end
                                R(current)=R(current)+ms(3)*M0(y);
                                R(current+1)=R(current+1)+ms(4)*M0(y);
                                R(current+2)=R(current+2)+0*M0(y);
                                R(current+3)=R(current+3)+0*M0(y);
    %                             current=base+O0(y)+hb;R(current)=R(current)+0*M0(y)+0*M0(y)+ms(4)*M0(y)+ms(3)*M0(y);%by qw
                            end
                        catch
                            disp('error');
                        end
                        y=y+1;
                    end
                    
                %final rows,no bottom bin
                for z=y:h0 %z相当于从y接着递增直到h0
                    yb0=floor(yb);eval(expression);
                    current=base;
                    if hasLf
                        R(current+O0(z))=R(current+O0(z))+ms(1)*M0(z);
                        R(current+O1(z))=R(current+O1(z))+ms(1)*M1(z);
                    end
                    if hasRt
                        R(current+O0(z)+hb)=R(current+O0(z)+hb)+ms(3)*M0(z);
                        R(current+O1(z)+hb)=R(current+O1(z)+hb)+ms(3)*M1(z);
                    end
                end
                
        end   %三线性插值的末尾
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