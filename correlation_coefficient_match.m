%Ӱ�����������������������ԣ�����Ƚϸ�
[filename_left,pathname_left]=uigetfile('.bmp','��ȡ����');%50Ϊ����49Ϊ����
image_left=imread(strcat(pathname_left,filename_left));
[im_height_left,im_width_left,degree_left]=size(image_left);
[filename_right,pathname_right]=uigetfile('.bmp','��ȡ����');
image_right=imread(strcat(pathname_right,filename_right));
[im_height_right,im_width_right,degree_right]=size(image_right);
imshow(image_left);
%a=image_left(2,3)ָ������������ y=2,x=3
[x_orgin,y_orgin]=ginput;%x�����ȣ�y����߶�
hold on;
plot(x_orgin,y_orgin,'r*');
%��5*5��Ŀ�괰��ȥ��������Ĵ󲿷�����
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
        %����Ŀ�괰�����������ڵ����ϵ��
        %�ȼ������������ڵ������ܺͺ�ƽ����
        Sgt=0;
        Sgtgt=0;
        Sggt=0;%�����ڶ�ӦԪ�س˻�֮��
       for i=y_right-(window_length-1)/2:y_right+(window_length-1)/2
           for j=x_right-(window_length-1)/2:x_right+(window_length-1)/2
               pixel=double(image_right(round(i),round(j)));
               Sgt=Sgt+pixel;
               Sgtgt=Sgtgt+pixel^2;
               pixel_ggt=double(pixel)*double(image_left(round(y_orgin-(y_right-i)),round(x_orgin-(x_right-j))));
               Sggt=Sggt+pixel_ggt;
           end
       end
       %�������ϵ��
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
%����x_target-1,y_target��x_target+1,y_target�����ϵ����Ȼ�����һ�����������
Sgt=0;
Sgtgt=0;
Sggt=0;%�����ڶ�ӦԪ�س˻�֮��
for i=y_target-(window_length-1)/2:y_target+(window_length-1)/2
     for j=x_target-1-(window_length-1)/2:x_target-1+(window_length-1)/2
           pixel=double(image_right(round(i),round(j)));
           Sgt=Sgt+pixel;
           Sgtgt=Sgtgt+pixel^2;
           pixel_ggt=double(pixel)*double(image_left(round(y_orgin-(y_target-i)),round(x_orgin-(x_target-1-j))));
           Sggt=Sggt+pixel_ggt;
     end
end
       %�������ϵ��
corre_p_before=(Sggt-Sg*Sgt/N)/(sqrt((Sgg-Sg^2/N)*(Sgtgt-Sgt^2/N)));
Sgt=0;
Sgtgt=0;
Sggt=0;%�����ڶ�ӦԪ�س˻�֮��
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
%��������ϵ�x����������ϣ�ʵ����ѡȡ��һ�㣬�������������ϵ���ֱ�Ϊ
%0.7583,0.9351,0.7436 �ɴ˽������������
x_target=x_target-(corre_p_next-corre_p_before)/(2*(corre_p_next+corre_p_before-2*correlation_p));
%�����������������������׼ȷ
plot(x_target,y_target,'b*');
hold on;



