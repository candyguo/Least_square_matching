    run('/Users/c_c/Documents/MATLAB/correlation_coefficient_match.m');
    halfwindow=9;
    image_l=im2double(image_left);
    image_r=im2double(image_right);
    x1=y_orgin;
    y1=x_orgin;
    x2=y_target;
    y2=x_target;
    %初值
    h0=0;h1=1;
    a0=x1-x2;a1=1;a2=0;
    b0=y1-y2;b1=0;b2=1;
    % a0=0;a1=1;a2=0;
    % b0=0;b1=0;b2=1;
    for i=1:2*halfwindow+1
        for j=1:2*halfwindow+1
            moban1(i,j)=bilinc(image_r,x2-halfwindow+i-1,y2-halfwindow+j-1);
        end
    end
    moban=zeros(halfwindow*2+3,halfwindow*2+3);
    x=zeros(halfwindow*2+3,halfwindow*2+3);
    y=zeros(halfwindow*2+3,halfwindow*2+3);
    for i=1:halfwindow*2+3
        for j=1:halfwindow*2+3
            %下面的if语句用来防止超出图像边界而报错
            if(a0+a1*(x2-halfwindow+i-2)+a2*(y2-halfwindow+j-2)>1 && a0+a1*(x2-halfwindow+i-2)+a2*(y2-halfwindow+j-2)<size(image_l,1)-1 && b0+b1*(x2-halfwindow+i-2)+b2*(y2-halfwindow+j-2)>1 && b0+b1*(x2-halfwindow+i-2)+b2*(y2-halfwindow+j-2)<size(image_l,2)-1)
                moban(i,j)=h0+h1*bilinc(image_l,a0+a1*(x2-halfwindow+i-2)+a2*(y2-halfwindow+j-2),b0+b1*(x2-halfwindow+i-2)+b2*(y2-halfwindow+j-2));
            end
            x(i,:)=x2-halfwindow-1;
            y(:,j)=y2-halfwindow-1;
        end
    end
    moban2=moban(2:size(moban,1)-1,2:size(moban,2)-1);
    x=x(2:size(moban,1)-1,2:size(moban,2)-1);
    y=y(2:size(moban,1)-1,2:size(moban,2)-1);
    %%
    %列误差方程式
    C(:,1)=linspace(1,1,(2*halfwindow+1)^2)';
    C(:,2)=reshape(moban2',(2*halfwindow+1)^2,1)';
    chafenx=zeros(2*halfwindow+1,2*halfwindow+1);    %行方向差分
    chafeny=zeros(2*halfwindow+1,2*halfwindow+1);    %列方向差分
    for i=1:2*halfwindow+1
        chafenx(i,:)=(moban(i+2,2:2*halfwindow+2)-moban(i,2:2*halfwindow+2))./2;
        chafeny(:,i)=(moban(2:2*halfwindow+2,i+2)-moban(2:2*halfwindow+2,i))./2;
    end
    C(:,3)=reshape(chafenx',(2*halfwindow+1)^2,1)';
    C(:,4)=reshape((x.*chafenx)',(2*halfwindow+1)^2,1)';
    C(:,5)=reshape((y.*chafenx)',(2*halfwindow+1)^2,1)';
    C(:,6)=reshape(chafeny',(2*halfwindow+1)^2,1)';
    C(:,7)=reshape((x.*chafeny)',(2*halfwindow+1)^2,1)';
    C(:,8)=reshape((y.*chafeny)',(2*halfwindow+1)^2,1)';
    L=reshape((moban1-moban2)',(2*halfwindow+1)^2,1);   
%     X=(C'*C)\(C'*L);
    
    %%  截断奇异值法解决病态
    [U,S,VVV]=svd(C'*C);
    for ii=1:size(S,1)
        if(S(1,1)/S(ii,ii)>10^8)
            break;
        end
    end
    VVVt=VVV(:,1:ii-1);
    Ut=U(:,1:ii-1);
    S=S(1:ii-1,1:ii-1);
    X=VVVt*inv(S)*Ut'*(C'*L);
    
    
    r_old=corr2(moban1,moban2);
    %%
    %下面进行迭代计算
    canshu=zeros(6,6);    %迭代6次，计算6个相关系数，取最大的
    r=zeros(1,6);
    Iteration=1;
    while(1)
        canshu(Iteration,:)=[a0 a1 a2 b0 b1 b2];
        r(1,Iteration)=r_old;
        %参数改正
        a0=a0+X(3,1)+a0*X(4,1)+b0*X(5,1);
        a1=a1+a1*X(4,1)+b1*X(5,1);
        a2=a2+a2*X(4,1)+b2*X(5,1);
        b0=b0+X(6,1)+a0*X(7,1)+b0*X(8,1);
        b1=b1+a1*X(7,1)+b1*X(8,1);
        b2=b2+a2*X(7,1)+b2*X(8,1);
        h0=h0+X(1,1)+h0*X(2,1);
        h1=h1+h1*X(2,1);
        for i=1:2*halfwindow+1
            for j=1:2*halfwindow+1
                moban1(i,j)=bilinc(image_r,x2-halfwindow+i-1,y2-halfwindow+j-1);
            end
        end
        moban=zeros(halfwindow*2+3,halfwindow*2+3);
        x=zeros(halfwindow*2+3,halfwindow*2+3);
        y=zeros(halfwindow*2+3,halfwindow*2+3);
        for i=1:halfwindow*2+3
            for j=1:halfwindow*2+3
                %下面的if语句用来防止超出图像边界而报错
                if(a0+a1*(x2-halfwindow+i-2)+a2*(y2-halfwindow+j-2)>1 && a0+a1*(x2-halfwindow+i-2)+a2*(y2-halfwindow+j-2)<size(image_l,1)-1 && b0+b1*(x2-halfwindow+i-2)+b2*(y2-halfwindow+j-2)>1 && b0+b1*(x2-halfwindow+i-2)+b2*(y2-halfwindow+j-2)<size(image_l,2)-1)
                    moban(i,j)=h0+h1*bilinc(image_l,a0+a1*(x2-halfwindow+i-2)+a2*(y2-halfwindow+j-2),b0+b1*(x2-halfwindow+i-2)+b2*(y2-halfwindow+j-2));
                end
                x(i,:)=x2-halfwindow-1;
                y(:,j)=y2-halfwindow-1;
            end
        end
        moban2=moban(2:size(moban,1)-1,2:size(moban,2)-1);
        x=x(2:size(moban,1)-1,2:size(moban,2)-1);
        y=y(2:size(moban,1)-1,2:size(moban,2)-1);
        %列误差方程式
        C(:,1)=linspace(1,1,(2*halfwindow+1)^2)';
        C(:,2)=reshape(moban2',(2*halfwindow+1)^2,1)';
        chafenx=zeros(2*halfwindow+1,2*halfwindow+1);    %行方向差分
        chafeny=zeros(2*halfwindow+1,2*halfwindow+1);    %列方向差分
        for i=1:2*halfwindow+1
            chafenx(i,:)=(moban(i+2,2:2*halfwindow+2)-moban(i,2:2*halfwindow+2))./2;
            chafeny(:,i)=(moban(2:2*halfwindow+2,i+2)-moban(2:2*halfwindow+2,i))./2;
        end
        C(:,3)=reshape(chafenx',(2*halfwindow+1)^2,1)';
        C(:,4)=reshape((x.*chafenx)',(2*halfwindow+1)^2,1)';
        C(:,5)=reshape((y.*chafenx)',(2*halfwindow+1)^2,1)';
        C(:,6)=reshape(chafeny',(2*halfwindow+1)^2,1)';
        C(:,7)=reshape((x.*chafeny)',(2*halfwindow+1)^2,1)';
        C(:,8)=reshape((y.*chafeny)',(2*halfwindow+1)^2,1)';
        L=reshape((moban1-moban2)',(2*halfwindow+1)^2,1);   
%         X=(C'*C)\(C'*L);
        
        %%  截断奇异值法解决病态
        [U,S,VVV]=svd(C'*C);
        for ii=1:size(S,1)
            if(S(1,1)/S(ii,ii)>10^8)
                break;
            end
        end
        VVVt=VVV(:,1:ii-1);
        Ut=U(:,1:ii-1);
        S=S(1:ii-1,1:ii-1);
        X=VVVt*inv(S)*Ut'*(C'*L);
    
    
        r_new=corr2(moban1,moban2);
        
        if(r_new>r_old)
            r_old=r_new;
            Iteration=Iteration+1;
        else
            break
        end
    end

    %% 对坐标做加权平均
    xx=zeros(2*halfwindow+1,2*halfwindow+1);
    yy=zeros(2*halfwindow+1,2*halfwindow+1);
    for i=1:2*halfwindow+1
        xx(i,1:2*halfwindow+1)=x2-(halfwindow+1-i);
        yy(1:2*halfwindow+1,i)=y2-(halfwindow+1-i);
    end
    sumx=sum(sum(chafenx.^2));sumy=sum(sum(chafeny.^2));
    xt=sum(sum(xx.*chafenx.^2/sumx));
    yt=sum(sum(yy.*chafeny.^2/sumy));
    [~,lie]=find(r==max(r));
    
%     xt=x2;
%     yt=y2;
    
    XS=canshu(lie,1)+canshu(lie,2)*xt+canshu(lie,3)*yt;
    YS=canshu(lie,4)+canshu(lie,5)*xt+canshu(lie,6)*yt;
    
    if(isnan(XS) || isnan(YS) || isnan(xt) || isnan(yt))
        XS=x1;YS=y1;xt=x2;yt=y2;
    end
   %% 改正后的相关系数
    halfwindow=7;
    for i=1:2*halfwindow+1
        for j=1:2*halfwindow+1
            if(min([XS-halfwindow+i-1 YS-halfwindow+j-1])>1 && XS-halfwindow+i-1<size(image_l,1)-1 && YS-halfwindow+j-1<size(image_l,2)-1)
                moban_I1(i,j)=bilinc(image_l,XS-halfwindow+i-1,YS-halfwindow+j-1);
            end
            if(min([xt-halfwindow+i-1 yt-halfwindow+j-1])>1 && xt-halfwindow+i-1<size(image_r,1)-1 && yt-halfwindow+j-1<size(image_r,2)-1)
                moban_I2(i,j)=bilinc(image_r,xt-halfwindow+i-1,yt-halfwindow+j-1);
            end
        end
    end
    if(isnan(XS) || isnan(YS) || isnan(xt) || isnan(yt) || size(moban_I1,1)~=size(moban_I2,1) || size(moban_I1,2)~=size(moban_I2,2))
        r_new=0;
    else
        r_new=corr2(moban_I1,moban_I2);
    end
    plot(yt,xt,'g*');
    