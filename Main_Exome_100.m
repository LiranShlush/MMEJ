function Main_Exome_100_Shay(chr)

addpath('C:\Postdoc\MMEJ detection\noa\find homologies/')
path1='C:\Postdoc\MMEJ detection\noa\find homologies/';
if(isnumeric(chr))
    chr=num2str(chr);
end

pathEx=[path1, 'Hg19Exome/exome_', chr, '.fa']; 
pathGe=[path1, '\Hg19Genome/', chr, '.fa'];

FASTAData = fastaread(pathEx);
ref=fastaread(pathGe);

%make the win size double the max del. size + homolog = 240bp, and the step
%size the size of the max del. (100bp)
del_max = 101;
homolog_size_max = 50;
homolog_size_min = 5;
%mismatch=0;
window_size = 240;
step_size = 100;

file_temp = ['C:\Postdoc\MMEJ detection\noa\find homologies/chr_' chr, '_Exome_mmej_deleteMe.txt'];
file_wrong = ['C:\Postdoc\MMEJ detection\noa\find homologies/chr_' chr '_Exome_mmej_wrong.txt'];

fid = fopen(file_temp,'w'); 
fid2 = fopen(file_wrong,'w');

fprintf(fid,'%s\t%s\t%s\n','Pos_1','Pos_2','Length');
fprintf(fid2,'%s\t%s\t%s\t%s\t%s\n','seq1','seq2','Pos_1','Pos_2','Length');

final_res = [];%change the var. to have a size

for s=1:length(FASTAData)
    
    disp(s);
    seq = FASTAData(s).Sequence;
    seq_length = length(seq);
    
    for i=1:step_size:seq_length-window_size
        
        if(seq_length>i+window_size)
            cur_seq = seq(i:i+window_size-1);
        else
            cur_seq = seq(i:seq_length);
        end
        
        [res] = find_homolog_in_str_Shay(cur_seq,homolog_size_min,window_size);%change the input var. for mismatch detection
        
        if ~isempty(res)
            %find real cordinates
            %extract the exome coor. from sequence name
            temp = strsplit(FASTAData(s).Header,':');
            pos = strsplit(temp{2},'-');
            start_pos = str2num(pos{1});
            %convert homolog coor to real coor
            res(:,1) = res(:,1)+(i-1)+start_pos;
            res(:,2) = res(:,2)+(i-1)+start_pos;
            
            %Shay
            %for every line in final_res extract the sequence and check homologies.
            %If homologies dont match write the sequence i and s
            for r=1:size(res,1)
                tempseq1=ref.Sequence(res(r,1):res(r,1)+res(r,3)-1);
                tempseq2=ref.Sequence(res(r,2):res(r,2)+res(r,3)-1);
                if ~strcmp(tempseq1,tempseq2)
                    fprintf(fid2,'%s\t%s\t%d\t%d\t%d\n',tempseq1,tempseq2,res(r,1),res(r,2),res(r,3));
                    %remove the line from the final file
                end
            end
            %Shay
            
            final_res =  res;
%            final_res = unique(final_res,'rows');
            for j=1:size(final_res,1)
                fprintf(fid,'%d\t%d\t%d\n',final_res(j,1),final_res(j,2),final_res(j,3));
            end
            final_res = [];
        end
    end
    
    
end

fclose all;

mat_final = readmatrix(file_temp);
mat_final = unique(mat_final,'rows');%remove duplicates

mat_final(:,1:2) = mat_final(:,1:2)-1;%corect coor. to fit Aditee's script 
mat_final = sortrows(mat_final);%sort

temp1=mat_final(:,2)-mat_final(:,1);%remove del longer than del_max
mat_final=mat_final(temp1<del_max,:);

temp1=(mat_final(:,1)+mat_final(:,3))<mat_final(:,2);%remove overlapping homologies
mat_final=mat_final(temp1,:);

file_final= ['C:\Postdoc\MMEJ detection\noa\find homologies/chr_' chr, '_Exome_mmej.txt'];
writematrix(mat_final,file_final,'Delimiter','tab')




