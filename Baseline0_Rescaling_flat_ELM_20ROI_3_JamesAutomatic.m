% flat ELM line from _AS images and thickness data, plot ISez profile and
% save as .png file 
% 2022-07-16, Haohua Qian
% 2022-10-14, modified based on Human_ISez
% 2022-10-17, modified for mouse selection
% 2022-10-28, center IS band in plot
% 2024-09-02 Added rescaling MCP/AR James Soetedjo. 
% 2024-09-27 Set y-axis scale between -25 to 120 James Soetedjo.
% 2024-09-28 Baseline correction is working, now need to figure out the
% x-axis scale James Soetedjo
% 2024-10-01 Baseline set to 0 then index of max point determined. Scaling
% based on non-corrected baseline.

% 2024-12-30 Added a check to make sure the end point is higher than the
% start point. If it's not, the angle is reversed to make the baseline
% rotation correct. This is rrMCPAR

% 2025-07-21: Most recent version of rrMCPAR (confirmed)

%[oct_file, oct_path]= uigetfile ('Stack_AS--ALL.img','Select OCT image file');  % oct image file
%Hard coding the file path
oct_file = 'Stack_AS--ALL.img';
oct_path = 'C:\Users\Public\ALL OCT CODE\MCPAR\to be and analyzed'; %INSERT FILE DIRECTORY HERE

oct = readanalyze(fullfile(oct_path,oct_file));
oct=oct(:,:,1);    % first frame
moct=mean(oct,2);  % mean intensity
first=find(moct>0,1);  % first non zero intensity
last=find(moct>0,1,'last');  % last non zero intensity

oct=imrotate(oct,90);  % rotate 90 degree for Analyze format

sini=20;   % initial IS band starting point
eini=30;   % initial IS band ending point

%Plot the first Ioct retina map and get the boundaries to calculate the
%minimums on each side. Will use i = 1. Therefore, i < 10
i = 1;
if i>10                % Superior
    left=last-199+(i-11)*20;
    right=last-200+(i-10)*20;
else
    left=first+(i-1)*20;    % Inferior
    right=first-1+i*20;
end

low=121; %Originally 121
suffix=sprintf('%02d', i);

roct=oct;
roct(170,left:right)=255;          % ROI
roct(low,left:right)=255;
roct(low:170,left)=255;
roct(low:170,right)=255;

Ioct=oct(low:170,left:right);

Int=mean(Ioct,2);  % average intensity profile

f=figure(2);
ymin=min(Int)-10;
f.Position = [50 300 1240 328];
plot(Int)
set(gca,'xlim',[0 60])
xticks([0:2:60])
set(gca,'Ylim',[ymin ymin+140])
[x_left,~, ~] = ginput(1); %x_left refers to the left boundary where we want to find our first minimum
hold on
plot([x_left x_left],[ymin ymin+140])
s=round(x_left);

[x_right,~, ~] = ginput(1); %x_right refers to the right boundary where we want to find our second minimum
plot([x_right x_right],[ymin ymin+140])
e=round(x_right);

a_limit = 2; %Prevent program from searching for e and s too far out of bounds

for i=1:20
    if i>10                % Superior
        left=last-199+(i-11)*20;
        right=last-200+(i-10)*20;
    else
        left=first+(i-1)*20;    % Inferior
        right=first-1+i*20;
    end
    low=121; %originally 121
    suffix=sprintf('%02d', i);
    
    roct=oct;
    roct(170,left:right)=255;          % ROI
    roct(low,left:right)=255;
    roct(low:170,left)=255;
    roct(low:170,right)=255;

    imwrite(uint8(roct),fullfile(oct_path, [oct_file(1:end-4) '_ROI_' suffix '.jpg']));

    Ioct=oct(low:170,left:right);

    Int=mean(Ioct,2);  % average intensity profile

    f=figure(2);
    ymin=min(Int)-10;
    f.Position = [50 300 1240 328];
    plot(Int)
    set(gca,'xlim',[0 60])
    xticks([0:2:60])
    set(gca,'Ylim',[ymin ymin+140])
    hold on
    plot([x_left x_left],[ymin ymin+140])
    %s=round(x); s was already defined outside the loop

    %[x,~, ~] = ginput(1); user input was already taken outside the loop
    plot([x_right x_right],[ymin ymin+140])
    %e=round(x); e was already defined outside the loop


    s=round(x_left); %Must redfine each time?
    e=round(x_right);

    center=(s+e)/2;    % center position of IS band
    hold off
    saveas(gcf,fullfile(oct_path, [oct_file(1:end-4) '_Profile_' suffix '.png']));
    
    slop=(Int(e)-Int(s))/(e-s);

    while Int(s)+slop-Int(s+1)>0  && Int(s)-slop-Int(s-1)<0  && e < (x_right + a_limit) % adjust start point
         if e>s+1
            s=s+1;
            slop=(Int(e)-Int(s))/(e-s);
         else
             break
         end
    end

    while Int(s)+slop-Int(s+1)<0  && Int(s)-slop-Int(s-1)>0 && e < (x_right + a_limit)  % adjust start point
            s=s-1;
            slop=(Int(e)-Int(s))/(e-s);
     end

    while Int(e)+slop-Int(e+1)<0  && Int(e)-slop-Int(e-1)>0 && e < (x_right + a_limit)  % adjust start point
         if e>s+1
            e=e-1;
            slop=(Int(e)-Int(s))/(e-s);
         else
             break
         end
    end

    while Int(e)+slop-Int(e+1)>0  && Int(e)-slop-Int(e-1)<0 && e < (x_right + a_limit)  % adjust start point
            e=e+1;
            slop=(Int(e)-Int(s))/(e-s);
     end

    IS=Int(s:e);
    

    f=figure(1);
    %f.Position = [100 50 620 328];
    
    %Correct Baseline

    %Point at which graph will be rotated. We will be using the variable
    %"IS" as IS is the region of interest. We do not use the variable Int
    %as this variable contains the entire region of interest (ELM to RPE)

    x_rotation_point = (e + s)/2; %smidpoint of MCP/AR
    y_rotation_point = (Int(s) + Int(e)) / 2; %y-value at starting point

    %Rotation angle: we will be using cos(theta). So to calculate angle,
    %use inverse cosine
    adjacent_length = e - s;
    hypotenuse_length = sqrt( (e - s)^2 + (IS(end) - IS(1))^2 );
    if IS(end) >= IS(1) %if the end point is higher than the starting point, then the angle we are referring to is on the left side
        theta_radians = acos(adjacent_length/hypotenuse_length);
    elseif IS(end) < IS(1) %if the end point is lower than the starting point, then the angle we are referring to is on the right side
        theta_radians = (-1) * acos(adjacent_length/hypotenuse_length);
    end

    %Rotation matrix: clockwise
    R = [cos(-theta_radians) -sin(-theta_radians); 
         sin(-theta_radians) cos(-theta_radians)];
    
    %Shift points to origin
    number_of_points = max(size(Int(s:e)));
    x_points = linspace(s, e, number_of_points);
    y_points = Int(s:e)'; %originally (e - s) x 1 vector so need to make it 1 x (e - s)
    
    shifted_points = [x_points - x_rotation_point; y_points - y_rotation_point]; %Shift graph to origin

    %Rotate the graph
    rotated_points = R * shifted_points;

    %Shift points back to original center
    new_x = rotated_points(1, :) + x_rotation_point;
    new_y = rotated_points(2, :) + y_rotation_point;
    
    %Next step: given the new_y values, find the max value and
    %corresponding index and use this index on the original IS values to
    %get the max. Then rescale. (START HERE).
    
    [M,Index_max] = max(new_y);

    %Comparison to make sure rotated graph is truly rotated relative to
    %original IS. Can comment out this section once confirmed.
    % plot(new_x, new_y, 'b');
    % hold on
    % plot([new_x(1) new_x(end)], [new_y(1) new_y(end)], 'b');
    % plot(s:e, IS, 'r')
    % plot([s, e], [IS(1) IS(end)], 'r')
    

    %Rescaling. Once we get the maximum index, we use that index on the IS
    %region.
    Int_region = Int(s-3:e+5);
    max_Int_region = IS(Index_max);
    min_Int_region = min(IS);
    Int_region = (Int_region - min_Int_region ) / (max_Int_region  - min_Int_region);
    Int_region = Int_region * 100;
    
    center=(s+e)/2;
    plot((s-3:e+5),Int_region,'k')
    hold on
    plot([s,e],Int_region([4, end - 5]),'k')    
    set(gca,'xlim',[center-12.5 center+12.5])
    %set(gca,'Ylim',[-25 120])
    set(gca,'Ylim',[-20 160])
    
    saveas(gcf,fullfile(oct_path, [oct_file(1:end-4) '_ISez_' suffix '.png']));
    hold off
    
    
    %writematrix(Int(s:e), fullfile(oct_path, [oct_file(1:end-4) '_ISez_' suffix '.csv']))
end

close all

