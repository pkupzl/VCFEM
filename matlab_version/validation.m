validation_points=[[0,0.95];[2.5,0.95];[0,1.45];[2.5,1.45];[0,1.95];[2.5,1.95];[0,2.5];[2.5,2.5]];
format long
dy = cell(1,8);
for i=1:8
    dy{i}=getDisplacement(validation_points(i,1),validation_points(i,2),d_m);
    
    disp(dy{i})
end