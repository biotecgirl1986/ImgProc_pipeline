clear all; %clear the memory from all the variables
close all; %close all the open figures
clc;       %reset the writing space
%%  LOAD THE MULTI-TIFF FILE FOR THE NUP CHANNEL
%%% to segment the nuclei in the Nup channel we consider the maximal
%%% z-projection
mov_folder='E:\Carmina data analysis\Tanguy videos\Delta Zelda 3M\HisRFP-noNLSMCP\12142015\3\movie_input\MCP\';               %indicate the full path movies file
spot_mov='MCP1-479';            %title of spot movie
nuclei_mov='HisRFPMax1-479';    %title of nuclei movie
convert_tiff=1;                      % Convert from multi-slice movie to normal movie 1 if you wanna convert

%frame for segmentation file
seg_init = 1;
seg_end  = 35;
%select the frame range you want to analyze
it_start=1;
it_end=3;


%% *****WRITE HERE GLOBAL PARAMETERS OF THE MOVIE***** %%%
shift_left=518;     %length of the cat portion on the left
shift_rigth=646;    %length of the cat portion on the rigth
dt=13.05;           %time step in seconds that you read from the metafile
A_pole=1;           %anterior pole position 1=left
start_cycle=10;     %starting cycle of the movie
z_max=16;           %number of maximal z
x_resolution=0.197; %pixel size (micron/pixels)

%% IMPORTANT PARAMETERES FOR SPOT DETECTION
th1=45;                 % minimal threshold for spot detection
th2=20;                 % minimal threshold for median filter
th=[th1 th2];


averaging_radius=3;     % there is an average filter for the images. averaging_radius=1 means no average. About the spot size
voxels_min=10;          % minimal number of voxels for a spot
voxels_max=50;          % Maximum number of voxels for a spot
fact_r=1.3;             % tolerance radius (if fact_r=1) the distance is equal to the radius

%% Load basic information on movies

%rescale the images factor equal to rescale
%if the resolution is high (images of 2042x2042 pixels) the nup channel can
%be rescaled a lot, also by a factor 0.25
%the obtained coordinate will be rescaled back
%this operation saves memory and computing time
ms2_mov = [mov_folder spot_mov];
seg_mov = [mov_folder nuclei_mov];
rescale=1; %the nup images are rescaled by a factor equal to rescale
x_resolution=x_resolution/rescale;

info = imfinfo([seg_mov '.tif']);     %infos about the file
num_images = numel(info);    %number of frames of the multitiff file

count_rescaling=0;

I=imread([seg_mov '.tif'],1);         %read the frame 1
I = imresize(I, rescale);             %rescaling the image

Lx=size(I,1);
Ly=size(I,2);
%% Load the segmentation image
load([seg_mov '_Track_' num2str(seg_end) '.mat'],'nuc','BWL');

%% SPOTS INFORMATIONS (compulsory)
%load the movie and store some variables regarding the images
mkdir('output_images');
Is=imread([ms2_mov '.tif'],1);

Lx2=size(Is,1);
Ly2=size(Is,2);

max_frame = z_max*it_end;

% Convert the spot movie if needed
if convert_tiff
    ms2_mov = conv_tiff(ms2_mov,Lx2,Ly2,max_frame);
end

window_around_spot=40;  %pixels, this is the window around the spot to do the gaussian fit, 40 pixels should be fine for the usual resolution 2048x2048

if it_end>seg_end   % If spot movie is longer than segmentation movie 
    it_end=seg_end;
end

Ispot_new=zeros(Lx2,Ly2,it_end-it_start+1);
Ispot_raw=zeros(Lx2,Ly2,it_end-it_start+1);
it14=0;

time_offset=0;
tic
for it=it_start:it_end
    
    s_start=(it-1)*z_max+1;
    s_end=(it-1)*z_max+z_max;    
    
    fprintf(1,'time %d\n',it);
 
    [spot,Ispot,raw_spot]=find_ms2_spots_1spot_adjusted(ms2_mov,it,nuc,th,z_max,voxels_min,voxels_max,fact_r,window_around_spot,averaging_radius,time_offset);

    %figure,imshow(Ispot);
    
    fprintf(1,'\n');
    Ispot_new(:,:,it)=Ispot;
    raw_spot_t(:,:,it)=raw_spot;
    raw_nuclei_t(:,:,it)=BWL(:,:,it);
    infos_spot{it}=spot;  
    if it==it_start
        writemod='overwrite';
    else
        writemod='append';
    end
    imwrite(uint8(Ispot_new(:,:,it)),['output_images/img_seg.tif'],'tif','compression','none','writemode',writemod); %save the segmented images in the folder segmented_spots
    imwrite(uint8(raw_spot_t(:,:,it)),['output_images/img_spot.tif'],'tif','compression','none','writemode',writemod); %save the segmented images in the folder segmented_spots
    imwrite([double(raw_spot_t(:,:,it))/256;double(Ispot+0.5*(raw_nuclei_t(:,:,it)>=1))],['output_images/img_nuclei.tif'],'tif','compression','none','writemode',writemod); %save the segmented images in the folder segmented_spots
end
toc
%% connect and allocate the infos regarding nuclei and spots (compulsory)

id_max=max(nuc.ind(:));

spots_number=1;   % maximal number of spots allowed
n_spots_infos=11; % number of variables for each spot

infos_spot_nuc=zeros(id_max,it_end,n_spots_infos,spots_number);

for it=it_start:it_end
    if isempty(infos_spot{it})==0
        for in=1:size(nuc.frames,1)
            id=nuc.ind(in,it);
            count_spot=0;
            
            clear ind_spot;
            
            for iis=1:size(infos_spot{it}.id_n,2)  %%%loop on the spots
                if infos_spot{it}.id_n(iis)==id    %%%to identify the nuclei
                    count_spot=count_spot+1;
                    ind_spot(count_spot)=iis;
                end
            end
            
            if count_spot>0
                ord_spot=zeros(count_spot,1);
                ord_spot_I=sort([infos_spot{it}.I(ind_spot(:))],'descend');
                ord_spot_A=sort([infos_spot{it}.A(ind_spot(:))],'descend');
                
                for ic=1:count_spot
                    ord_spot(ic)=find(infos_spot{it}.A(:)==ord_spot_A(ic));
                end
                
                for i=1:count_spot
                    for j=1:count_spot
                        jj=ord_spot(i);
                        infos_spot_nuc(id,it,1,i)=infos_spot{it}.x(jj);
                        infos_spot_nuc(id,it,2,i)=infos_spot{it}.y(jj);
                        infos_spot_nuc(id,it,3,i)=infos_spot{it}.z(jj);
                        infos_spot_nuc(id,it,4,i)=infos_spot{it}.size(jj);
                        infos_spot_nuc(id,it,5,i)=infos_spot{it}.I(jj);
                        infos_spot_nuc(id,it,6,i)=infos_spot{it}.A(jj);
                        infos_spot_nuc(id,it,7,i)=infos_spot{it}.ssx(jj);
                        infos_spot_nuc(id,it,8,i)=infos_spot{it}.ssy(jj);
                        infos_spot_nuc(id,it,9,i)=infos_spot{it}.I2d(jj);
                        infos_spot_nuc(id,it,10,i)=infos_spot{it}.bkg(jj);
                        infos_spot_nuc(id,it,11,i)=infos_spot{it}.resid(jj);
                    end
                end
                
            end
        end
    end
end

%% print a table with all the infos for each nucleus in each time step (compulsory)
mkdir('table_summary');
filename = ['table_summary/th' num2str(th) '.txt'];
fp = fopen(filename,'w+');


id_max=max(nuc.ind(:));
% ir: index of nuclei
% id: id of nuclei
% it: time

for ir=1:size(nuc.frames,1)
    for it=it_start:it_end 
        id=nuc.ind(ir,it);
        if nuc.frames(ir,it)
            %if id<=size(v_frame(it).vs.ind,2) %in the case a nucleus was lost
                % Print the data
                daughter=nuc.daughter{ir,it};
                daughter_=zeros(1,2);
                for i=1:numel(daughter)
                    daughter_(i)=daughter(i);
                end
                fprintf(1,'%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',id,it,nuc.positions{ir,it}(1),nuc.positions{ir,it}(2),nuc.radius(ir,it),nuc.ecc(ir,it),nuc.ecc(ir,it),nuc.ind(ir,it),nuc.parent(ir,it),daughter_(1),daughter_(2),nuc.cycle(ir,it),dt,Lx2,Ly2,shift_left,shift_rigth,th,A_pole);
                fprintf(fp,'%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',id,it,nuc.positions{ir,it}(1),nuc.positions{ir,it}(2),nuc.radius(ir,it),nuc.ecc(ir,it),nuc.ecc(ir,it),nuc.ind(ir,it),nuc.parent(ir,it),daughter_(1),daughter_(2),nuc.cycle(ir,it),dt,Lx2,Ly2,shift_left,shift_rigth,th,A_pole);
                for j=1:spots_number
                    for i=1:n_spots_infos
                        fprintf(1,'s %g ',infos_spot_nuc(id,it,i,j));
                        fprintf(fp,'%g ',infos_spot_nuc(id,it,i,j));
                    end
                end
                fprintf(1,'\n');
                fprintf(fp,'\n');
            %end
        end
    end
end
fclose(fp);