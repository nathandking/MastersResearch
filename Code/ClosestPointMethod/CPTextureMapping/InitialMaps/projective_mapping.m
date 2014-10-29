function [Y_newxy] = projective_mapping(Y_2d, Ps1, Ps2, Ps3, Ps4, Pd1, Pd2, Pd3, Pd4)

% implement projective mapping for quadrilateral-quadrilateral as outlines
% in Heckbert 1989.



M = zeros(8,8);
Msd = zeros(3,3);
Y_newxy = zeros(size(Y_2d,1),2);
idxxx = zeros(size(Y_2d,1),1);
for i = 1:size(Y_2d,1)
   M(1:4, 1:2) = [Ps1(i,:); Ps2(i,:); Ps3(i,:); Ps4(i,:)];
   M(5:8, 4:5) = [Ps1(i,:); Ps2(i,:); Ps3(i,:); Ps4(i,:)];
   M(1:4, 3) = ones(4,1);
   M(5:8, 6) = ones(4,1);
   M(:,7) = -[Ps1(i,1)*Pd1(i,1); Ps2(i,1)*Pd2(i,1); Ps3(i,1)*Pd3(i,1);...
       Ps4(i,1)*Pd4(i,1); Ps1(i,1)*Pd1(i,2); Ps2(i,1)*Pd2(i,2);...
       Ps3(i,1)*Pd3(i,2); Ps4(i,1)*Pd4(i,2)];
   M(:,8) = -[Ps1(i,2)*Pd1(i,1); Ps2(i,2)*Pd2(i,1); Ps3(i,2)*Pd3(i,1);...
       Ps4(i,2)*Pd4(i,1); Ps1(i,2)*Pd1(i,2); Ps2(i,2)*Pd2(i,2);...
       Ps3(i,2)*Pd3(i,2); Ps4(i,2)*Pd4(i,2)];
   
   b = [Pd1(i,1); Pd2(i,1); Pd3(i,1); Pd4(i,1); Pd1(i,2); Pd2(i,2); Pd3(i,2); Pd4(i,2)];
   C = M\b;
   
   Msd(:,1) = C(1:3);
   Msd(:,2) = C(4:6);
   Msd(:,3) = [C(7:8); 1];

   Y_newxy_temp = [Y_2d(i,:), 1]*Msd;
   Y_newxy(i,:) = Y_newxy_temp(1:2);
   
    [warnmsg, msgid] = lastwarn;
    
    if strcmp(msgid,'MATLAB:singularMatrix')
        Y_newxy(i,:) = [NaN,NaN];
        lastwarn('');
        idxxx(i) = 1;
        M;
    end
    
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        Y_newxy(i,:) = [NaN,NaN];
        lastwarn('');
        idxxx(i) = 1;
        M;
    end
    
    if strcmp(msgid,'MATLAB:nearlySingularMatrix')
        Y_newxy(i,:) = [NaN,NaN];
        lastwarn('');
        idxxx(i) = 1;
        M;
    end
    
    if Y_newxy(i,1) <= min([Pd1(i,1), Pd2(i,1), Pd3(i,1), Pd4(i,1)])
       Y_newxy(i,:) = [NaN,NaN];
       
    end
    if Y_newxy(i,2) <= min([Pd1(i,2), Pd2(i,2), Pd3(i,2), Pd4(i,2)])
       Y_newxy(i,:) = [NaN,NaN];
       
    end
    
    if Y_newxy(i,1) >= max([Pd1(i,1), Pd2(i,1), Pd3(i,1), Pd4(i,1)])
       Y_newxy(i,:) = [NaN,NaN];
       
    end
    if Y_newxy(i,2) >= max([Pd1(i,2), Pd2(i,2), Pd3(i,2), Pd4(i,2)])
       Y_newxy(i,:) = [NaN,NaN];
       
    end
    

    
end


