clc();
% a = [1 1 1 3 5];
% b = [1 1 3 3 5];
% c = [1 1 3 4 5];
% d = [2 1 4 5 5];
% e = [2 2 1 4 5];
% 
% obj = MATLAB_DBA.DTW_distance(a,d,100);
% ans_d = obj.Similarity_answer
% obj = MATLAB_DBA.DTW_distance(a,e,100);
% ans_e = obj.Similarity_answer
% 
% sab = ((1+1+1+9+9+25)-0)/(1+1+1+9+9+25)
% sac = ((1+1+1+9+9+25)-1)/(1+1+1+9+9+25)
% sad = ((1+1+1+9+25+25)-2)/(1+1+1+9+25+25)
% sae = ((1+1+1+9+25)-3)/(1+1+1+9+25)

raw_sim = [];
raw_dis = [];
mean_sim = [];
mean_dis = [];
sd_sim = [];
sd_dis = [];

name = importdata('C:\Users\sura\Desktop\Dataset\UCR 47 dataset\47_Trainingset_Name.txt');
[number_of_file nc] = size(name);
Y = []; %% use in plot error graph
E = []; %% use in plot error graph
kendall_tau = [];
for round=1:number_of_file
    
file = {'C:\Users\sura\Desktop\Dataset\UCR 47 dataset\'};
file = strcat(file,name{round});
file=char(file)

Data = importdata(file);
[row,column] = size(Data);
ans_dtw = zeros(row,row);
ans_similar = zeros(row,row);
for i=1:row     %% i is candidate   &   j is query
    for j=1:row
        obj = MATLAB_DBA.DTW_distance(Data(j,2:column),Data(i,2:column),100);
        ans_dtw(i,j) = obj.DTW_answer();
        ans_similar(i,j) = obj.Similarity_answer();
    end
end
fdtw_name = strcat('C:\Users\sura\Desktop\BCI_workspace\New_Similarity\',strcat(name{round},'_TRAIN_DTW_dist_column_query_row_cand.txt'));
fdtwsim_name = strcat('C:\Users\sura\Desktop\BCI_workspace\New_Similarity\',strcat(name{round},'_TRAIN_DTW_sim_column_query_row_cand.txt'));
dlmwrite(fdtw_name,ans_dtw,'delimiter','\t');
dlmwrite(fdtwsim_name,ans_similar,'delimiter','\t');


% [B_d,I_d] = sort(ans_dtw);
% [B_s,I_s] = sort(ans_similar,'descend');
% 
% raw_dis = [raw_dis, B_d(1:4,1)'];
% mean_dis = [mean_dis, mean(B_d(1:4,1))];
% sd_dis = [sd_dis, std(B_d(1:4,1))];
% raw_sim = [raw_sim, B_s(1:4,1)'];
% mean_sim = [mean_sim, mean(B_s(1:4,1))];
% sd_sim = [sd_sim, std(B_s(1:4,1))];






% %%%%%%%%%%%%%%%%% Kandell tau rank correlation ......    using I_d as reference
% rank_a = (1:row)';
% rank_b = [];
% for i=1:row
%     rank_b = [rank_b;find(I_s == I_d(i))];
% end
% rank_b;
% concordant = 0;
% discordant = 0;
% for i=1:row-1
%     concordant = concordant + size(find(rank_b(i:end)>rank_b(i)),1);
%     discordant = discordant + size(find(rank_b(i:end)<rank_b(i)),1);
% end
% kendall_tau = [kendall_tau, (concordant-discordant) / (concordant+discordant)];
% fprintf('kendall tau = %f\n\n',kendall_tau(round));

end

% fprintf('mean_dis = %f , sd = %f\n',mean(raw_dis),std(raw_dis));
% fprintf('mean_sim = %f , sd = %f\n',mean(raw_sim),std(raw_sim));

% figure();
% errorbar([1:47],mean_dis,sd_dis,'x');
% ylabel('mean +- sd of the dtw for the first 4 rank that most similar');
% xlabel('dataset');
% xlim([0 47]);
% saveas(gcf,'C:\Users\sura\Desktop\BCI_workspace\รูป Threshold testing\dtw mean_sd of the dtw for the first 4 rank that most similar','epsc');

% figure();
% errorbar([1:47],mean_sim,sd_sim,'x');
% ylabel('mean +- sd of the similarity for the first 4 rank that most similar');
% xlabel('dataset');
% xlim([0 47]);
% ylim([0 1.0]);
% saveas(gcf,'C:\Users\sura\Desktop\BCI_workspace\รูป Threshold testing\similarity mean_sd of the dtw for the first 4 rank that most similar','epsc');


% figure;
% plot(kendall_tau,'x');
% ylim([0.0 1.0]);
% ylabel('Kendall tau correlation');
% xlabel('dataset');
% xlim([0 47]);
% saveas(gcf,'C:\Users\Tong\Desktop\BCI_workspace\รูป Threshold testing\Kendall tau correlation','epsc');
% fprintf('mean kendall tau = %f , sd kendall_tau = %f\n\n',mean(kendall_tau),std(kendall_tau));
