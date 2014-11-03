close all
%i = YNan(5);

% figure(1);
% scatter3([X(IDX(i,:),1)],[X(IDX(i,:),2)], [X(IDX(i,:),3)],20,'fill')
% hold on
% scatter3(Y(i,1),Y(i,2), Y(i,3),20,'fill')
% scatter3(Y_3d(i,1),Y_3d(i,2), Y_3d(i,3),20,'fill')
% scatter3([proj_S1(i,1), proj_S2(i,1), proj_S3(i,1), proj_S4(i,1)],[proj_S1(i,2), proj_S2(i,2), proj_S3(i,2), proj_S4(i,2)],[proj_S1(i,3), proj_S2(i,3), proj_S3(i,3), proj_S4(i,3)],20)

% for j = 1:size(YNan,1)
%     i = YNan(j);
%     figure(2);
%     scatter(Y_2d(i,1),Y_2d(i,2),50,'fill')
%     hold on
%     scatter([proj_S1_2d(i,1), proj_S2_2d(i,1), proj_S3_2d(i,1), proj_S4_2d(i,1)],[proj_S1_2d(i,2), proj_S2_2d(i,2), proj_S3_2d(i,2), proj_S4_2d(i,2)],50,'fill')
%     
% %     figure(1);
% %     scatter3(X(IDX(i,:),1),X(IDX(i,:),2), X(IDX(i,:),3),50,'fill')
% %     hold on
% %     scatter3(Y(i,1),Y(i,2), Y(i,3),50,'fill')
%     
%     pause(1.5)
%     %clf(1);
%     clf(2);
% end


for i = 1:size(patches,1)
    if size(patches{i},2)==2
        figure(1)
        scatter(patches{i}(:,1),patches{i}(:,2),40,'fill')
        hold on
        scatter(cpX_tan_2d{i}(1),cpX_tan_2d{i}(2),40,'fill')
        pause(1.5)
        clf(1)
    end
end