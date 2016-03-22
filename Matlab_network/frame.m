clear all, close all,

%jani = csvread('RB1 orientation test.csv');

%quatRB1 = jani(:, 3:6);

permSequence = perms(1:3);

for pms = 1:6
    
    for aadd = -90:90:90
        for badd = -90:90:90
            for cadd = -90:90:90

for amult = -1:2:1
    for bmult = -1:2:1
        for cmult = -1:2:1
for i = 1:3
    for j = 1:3
        for k = 1:3
            
            if i == 1
                str1 = 'X';
            elseif i == 2
                str1 = 'Y';
            else
                str1 = 'Z';
            end
                
            if j == 1
                str2 = 'X';
            elseif j == 2
                str2 = 'Y';
            else
                str2 = 'Z';
            end
            
            if k == 1
                str3 = 'X';
            elseif k == 2
                str3 = 'Y';
            else
                str3 = 'Z';
            end
            
            angles = [-13 -79 161]/180*pi;
            angles = angles([permSequence(pms, :)]);
            angles(1) = angles(1) + aadd;
            angles(2) = angles(2) + badd;
            angles(3) = angles(3) + cadd;
            
            strSequence = [str1, str2, str3];
            multiSequence = num2str([amult, bmult, cmult]);
            
            try
                q = angle2quat(angles(1) *amult, angles(2) * bmult, angles(3)* cmult, strSequence); 
                q_ref = [0.63 0.01 -0.76 -0.05];
                qdiff = norm(q-q_ref);
                if qdiff < 0.1

                    disp(' ')
                    disp([strSequence, ', angleMultiplier: ', multiSequence, ', anglePermutation: ', num2str(permSequence(pms, :)), ', angle addition: ', num2str([aadd, badd, cadd])])
                    disp(['err: ', num2str(qdiff), ', quaternion:', num2str(q)])
                end
            catch msg
               
                %disp(['sequence: ', strSequence, ' not possible'])
            end
            
        end
    end
end
        end
    end
end
            end
        end
    end
end
              