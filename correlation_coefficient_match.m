%影像的纹理机构清晰，特征明显，信噪比较高
[filename_left,pathname_left]=uigetfile('.bmp','读取左像');%50为左像，49为右像
image_left=imread(strcat(pathname_left,filename_left));
[im_height_left,im_width_left,degree_left]=size(image_left);
[filename_right,pathname_right]=uigetfile('.bmp','读取右像');
image_right=imread(strcat(pathname_right,filename_right));
[im_height_right,im_width_right,degree_right]=size(image_right);
imshow(image_left);
%a=image_left(2,3)指的是两行三列 y=2,x=3
[x_orgin,y_orgin]=ginput;%x代表宽度，y代表高度
hold on;
plot(x_orgin,y_orgin,'r*');
%用5*5的目标窗口去搜索右像的大部分区域
correlation_p=0;
x_target=0;
y_target=0;
window_length=9;
N=window_length*window_length;
Sg=0;
Sgg=0;
for i=y_orgin-(window_length-1)/2:y_orgin+(window_length-1)/2
    for j=x_orgin-(window_length-1)/2:x_orgin+(window_length-1)/2
        pixel=double(image_left(round(i),round(j)));
        Sg=Sg+pixel;
        Sgg=Sgg+pixel^2;
    end
end
for x_right=1+(window_length-1)/2:0.75*im_width_right
    for y_right=1+(window_length-1)/2:im_height_right-(window_length-1)/2
        %计算目标窗口与搜索窗口的相关系数
        %先计算搜索窗口内的像素总和和平方和
        Sgt=0;
        Sgtgt=0;
        Sggt=0;%两窗口对应元素乘积之和
       for i=y_right-(window_length-1)/2:y_right+(window_length-1)/2
           for j=x_right-(window_length-1)/2:x_right+(window_length-1)/2
               pixel=double(image_right(round(i),round(j)));
               Sgt=Sgt+pixel;
               Sgtgt=Sgtgt+pixel^2;
               pixel_ggt=double(pixel)*double(image_left(round(y_orgin-(y_right-i)),round(x_orgin-(x_right-j))));
               Sggt=Sggt+pixel_ggt;
           end
       end
       %计算相关系数
       corre_p=(Sggt-Sg*Sgt/N)/(sqrt((Sgg-Sg^2/N)*(Sgtgt-Sgt^2/N)));
       if(corre_p>correlation_p)
           correlation_p=corre_p;
           x_target=x_right;
           y_target=y_right;
       end
    end
end
figure(2);
imshow(image_right);
hold on;
%计算x_target-1,y_target与x_target+1,y_target的相关系数，然后进行一个抛物线拟合
Sgt=0;
Sgtgt=0;
Sggt=0;%两窗口对应元素乘积之和
for i=y_target-(window_length-1)/2:y_target+(window_length-1)/2
     for j=x_target-1-(window_length-1)/2:x_target-1+(window_length-1)/2
           pixel=double(image_right(round(i),round(j)));
           Sgt=Sgt+pixel;
           Sgtgt=Sgtgt+pixel^2;
           pixel_ggt=double(pixel)*double(image_left(round(y_orgin-(y_target-i)),round(x_orgin-(x_target-1-j))));
           Sggt=Sggt+pixel_ggt;
     end
end
       %计算相关系数
corre_p_before=(Sggt-Sg*Sgt/N)/(sqrt((Sgg-Sg^2/N)*(Sgtgt-Sgt^2/N)));
Sgt=0;
Sgtgt=0;
Sggt=0;%两窗口对应元素乘积之和
 for i=y_target-(window_length-1)/2:y_target+(window_length-1)/2
    for j=x_target+1-(window_length-1)/2:x_target+1+(window_length-1)/2
          pixel=double(image_right(round(i),round(j)));
          Sgt=Sgt+pixel;
          Sgtgt=Sgtgt+pixel^2;
          pixel_ggt=double(pixel)*double(image_left(round(y_orgin-(y_target-i)),round(x_orgin-(x_target+1-j))));
          Sggt=Sggt+pixel_ggt;
    end
 end
corre_p_next=(Sggt-Sg*Sgt/N)/(sqrt((Sgg-Sg^2/N)*(Sgtgt-Sgt^2/N)));
%抛物线拟合的x方向坐标拟合，实验中选取了一点，其横向三点的相关系数分别为
%0.7583,0.9351,0.7436 由此进行抛物线拟合
x_target=x_target-(corre_p_next-corre_p_before)/(2*(corre_p_next+corre_p_before-2*correlation_p));
%抛物线修正后的像点坐标更加准确
plot(x_target,y_target,'b*');
hold on;



