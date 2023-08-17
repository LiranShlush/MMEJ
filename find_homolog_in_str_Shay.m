function [final_res] = find_homolog_in_str_Shay(seq,min_homolog_size,window_size)
%for mismatch get a different input for the 2nd 'min_homolog_size' in my_seqdotplot

[~, Matrix] = my_seqdotplot(seq,seq,min_homolog_size,min_homolog_size);%change the last parameter for mismatches
seq_len=length(seq);
mat1=zeros(size(Matrix));
%only half a matrix is needed
[row,col]=find(Matrix);
col2=col(col>row);
row2=row(col>row);
%Remove homologs bigger than the window_size 
temp1=(col2-row2<window_size);
col2=col2(temp1);
row2=row2(temp1);
%write only homologies in the win. size in the half matrix to new matrix
mat1(sub2ind(size(mat1), row2,col2))=Matrix(sub2ind(size(Matrix), row2,col2));

mat2=mat1*100;
matLongHomolog=mat1(1:end-1,1:end-1)-mat2(2:end,2:end);

%Shay

%remove truncated homologs
%Find -99 in x=1 and y=1 and remove the end 1 of this -99 (on the diagonal)
[xh, yh]=find(matLongHomolog(1,2:end)==-99);
yh=yh+1;
if~isempty(yh)
    for i=1:length(yh)
        for j=1:(seq_len-yh)
            if(matLongHomolog(xh(i)+j,yh(i)+j)==1)
                matLongHomolog(xh(i)+j,yh(i)+j)=0;
                matLongHomolog(yh(i)+j,xh(i)+j)=0;
                break
            end
        end
    end
end
%remove homolog if 1 is on the end
h=matLongHomolog(1,:)==1;
matLongHomolog(1,h)=0;
matLongHomolog(h,1)=0;
%remove 1 from the end col and the corresponding -100



%for every '-100'count until found 1
[xh, yh]=find(matLongHomolog==-100);
[x2,y2]=find(matLongHomolog==1);%end of the homolog

%find where the x and y coor are of the same size and only then pick
%minimal distance
xDif=x2'-xh;
xDif(xDif<1)=-500;%remove negative unrelevant results
yDif=y2'-yh;
yDif(yDif<1)=-1500;%remove negative unrelevant results
[row,col]=find(xDif==yDif);

tempMat=zeros(size(xDif));
tempMat=tempMat+5000;
tempMat(sub2ind(size(tempMat), row,col))=xDif(sub2ind(size(tempMat), row,col));
tempMat=min(tempMat,[],2);
homolog_len=tempMat+min_homolog_size-1;
res= [xh+1, yh+1, homolog_len];

final_res=res;

%final_res=sortrows(res,1);

% 
% gg=sum(xyDif);
% for i=1:length(gg)%for each of the col. find the minimal value
%     if(gg(i)==1)
%         
%     end
% end
% 
% homolog_len=zeros(size(yh));
% homolog_len=homolog_len+min_homolog_size-1;
% 
% for i=1:length(yh)
%     y=yh(i);
%     x=xh(i);
%     counter=0;
%     while(matLongHomolog(x,y)~=1)
%         y=y+1;
%         x=x+1;
%         counter=counter+1;
%     end
%     homolog_len(i)=homolog_len(i)+counter;
% end
% 
% res= [xh+1, yh+1, homolog_len];
% 
% final_res=sortrows(res,1);
% 
% %Shay
% 
% 
% 
% %remove truncated homologs
% %Find -99 in x=1 and y=1 and remove the end 1 of this -99 (on the diagonal)
% [xh, yh]=find(matLongHomolog(1,2:end)==-99);
% yh=yh+1;
% if~isempty(yh)
%     for i=1:length(yh)
%         for j=1:(100-yh)%magic numbers change the 100 to lenght(seq)
%             if(matLongHomolog(xh(i)+j,yh(i)+j)==1)
%                 matLongHomolog(xh(i)+j,yh(i)+j)=0;
%                 matLongHomolog(yh(i)+j,xh(i)+j)=0;
%                 break
%             end
%         end
%     end
% end
% %remove homolog if 1 is on the end
% h=matLongHomolog(1,:)==1;
% matLongHomolog(1,h)=0;
% matLongHomolog(h,1)=0;
% %remove 1 from the end col and the corresponding -100
% 
% 
% 
% %find the start and end of each homolog
% [x1,y1]=find(matLongHomolog==-100);%start of the homolog
% [x2,y2]=find(matLongHomolog==1);%end of the homolog
% 
% 
% %find where the x and y coor are of the same size and only then pick
% %minimal distance
% xDif=x2'-x1;
% xDif(xDif<1)=-500;%remove negative unrelevant results
% yDif=y2'-y1;
% yDif(yDif<1)=-1500;%remove negative unrelevant results
% xyDif=xDif==yDif;
% 
% %change to a loop!!!! only for non standart cases
% xyDifSum=sum(xyDif)>1;
% if(any(xyDifSum))
%     xyDifSumCoorX=find(xyDifSum);%find which colomns have more than sinal value
%     for i=1:length(xyDifSumCoorX)%for each of the col. find the minimal value
%         tempCol=xDif(xyDif(:,xyDifSumCoorX(i)),xyDifSumCoorX(i));%take the distance values
%         tempCol=min(tempCol);%pick the minimal distance
%         xyDif(:,xyDifSumCoorX(i))=0;%put zero in ALL the places 
%         tempCol=find(xDif(:,xyDifSumCoorX(i))==tempCol, 1);
%         xyDif(tempCol,xyDifSumCoorX(i))=1;%return 1 to the place where the minimal value
%     end
% end
% 
% xDifMin=xDif(xyDif);
% %yDifMin=yDif(xyDif);
% 
% resgg= [x1+1, y1+1, min_homolog_size+xDifMin-1];
% 
% %final_res=sortrows(resgg,1);
% 
% final_res=resgg;
