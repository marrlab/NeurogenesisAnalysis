function [] = plotPieCharts(w,model_str,lab1,lab2,opt_2,data_str)

ind_vec = 1:length(lab2);
pie_colors_short = [3 6 97; 28 60 164; 153 204 255;   %NSC colors
          114 155 155; 0 204 204; 164 246 246; %TAP colors
          118 62 99; 209 106 175; 247 203 233]./256; %NB colors
% pie_colors_long = [220 240 246; 127 181 222; 17 79 127; 10 45 71;  %NSC colors
%           201 233 243; 104 165 183; 54 103 118; 24 63 75; %TAP colors
%           236 206 247; 183 134 201; 107 75 118; 60 28 71]./256; %NB colors
pie_colors_long = [ 173 208 249; 255 236 188; 192 246 247; 249 205 216;
                    173 208 249; 255 236 188; 192 246 247; 249 205 216;
                    173 208 249; 255 236 188; 192 246 247; 249 205 216]/256;
% pie_colors_long = [ 31 88 152; 239 142 32; 13 146 151; 152 31 64;
%                     31 88 152; 239 142 32; 13 146 151; 152 31 64;
%                     31 88 152; 239 142 32; 13 146 151; 152 31 64]/256;
switch length(lab2)
  case 3
      pie_colors=pie_colors_short;
  case 4
      pie_colors=pie_colors_long;
end

P=zeros(length(lab1),length(lab2));
for i=1:length(lab1)
    for j=1:length(lab2)
        switch lab2{j}(1)
            case 'A'
                lab2_str = 'a';
            case 'S'
                lab2_str = 's';
            case 'C'
                lab2_str = 'd';
            case 'U'
                lab2_str = 'as';
            otherwise
                error('label name unknown!');
        end
        search_str = strcat('*',lab1(i),'_',lab2_str,'t*');
        idx = find(strwcmp(model_str,search_str));
        P(i,j) = sum(w(idx))/sum(w);
    end
    figure()
    set(gca,'fontsize',15)
    c=subplot(1,1,1);
    p=pie(P(i,:));
    for j=1:length(lab2)
        t=p(2*j);
        t.FontSize=14;
    end
    title(c,lab1{i})
    l=legend(lab2,'Location','southoutside','Orientation','vertical');
    set(l,'FontSize',14);
    colormap(pie_colors(length(lab2)*(i-1)+1:length(lab2)*(i),:));
    disp(lab1(i));
    saveFigs(opt_2,['DivStrategyPie_',data_str,'_',lab1{i}]);
end


end