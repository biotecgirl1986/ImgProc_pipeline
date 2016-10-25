function SEGNUP3
%% Author
% Mathieu Coppey
% Modified by Huy Tran
%% ideas
% http://fr.mathworks.com/help/images/ref/fspecial.html?refresh=true
% http://www.mathworks.com/matlabcentral/answers/57077-how-do-you-perform-a-difference-of-gaussian-filter-on-an-image
% http://www.activovision.com/octavi/doku.php?id=integral_images
% http://fr.mathworks.com/help/images/examples/marker-controlled-watershed-segmentation.html
% http://fr.mathworks.com/help/images/ref/watershed.html
%% Parameters for positions and colors of figures and some handles
% The position of the game figure will be [figX1,figX2,figY1,figY2]
figX1 = 0;
figX2 = 1300;
figY1 = 0;
figY2 = 800;

%% sets initial variables
FileName = '';
PathName = '';
FilterIndex = '';
frame = 1;
frameinit = 1;
framefin = 150;
Nf = 150;
showBW = 0;
Img = [];
BW = [];
BWL = [];
Nnuc = [];
nucActive = [];
open = 1;
seg = 0;

nuc.frames = [];
nuc.positions = [];
nuc.pixels = [];
nuc.cycle = [];
nuc.name = [];
nuc.box = [];
nuc.filiation = [];
        
%% Create and hide the GUI figure as it is being constructed.
segfigure = figure('Visible','on','Tag','segfigure','Position',[figX1,figY1,figX2,figY2]);
set ( gcf, 'Color', [0 0 0] )

% To capture mouse and keyboard
set(segfigure,'WindowKeyPressFcn',@keypress);
set(segfigure,'WindowScrollWheelFcn',@scroll);
set(segfigure,'WindowButtonUpFcn',@click);
%% Create buttons and others
hmessages = uicontrol('Style','text','String','Please load a NUP movie',...
    'Position',[50,730,350,50],'FontSize',14, ...
    'ForegroundColor','white','BackgroundColor','black'); % Display instructions

hquitbutton = uicontrol('Style','pushbutton',...
    'String','Quit','Callback',@hquitbutton_Callback,...
    'Position',[1250,750,50,50]); % Quit seg
%20.6628 -446.1809  144.5822   52.1283

hloadbutton = uicontrol('Style','pushbutton',...
    'String','Load','Callback',@hloadbutton_Callback,...
    'Position',[400,730,100,50]); % load file

hpath = uicontrol('Style','text','String','Path',...
    'Position',[510,730,400,50],'FontSize',10);

hframe = uicontrol('Style','slider',...
    'Min',1,'Max',Nf, 'value',frame,...
    'Position',[20,10,1200,20],...
    'SliderStep',[1/150 0.2],'Callback',@hframe_Callback);

hopen = uicontrol('Style','slider',...
    'Min',1,'Max',30, 'value',open,...
    'Position',[50,670,100,50],...
    'SliderStep',[1/30 0.2],'Callback',@hopen_Callback);

hchoice = uicontrol('Style', 'popup',...
           'String', 'soft|hard',...
           'Position', [50 720 60 30],...
           'Callback', @hchoice_Callback);   

hsegment = uicontrol('Style','pushbutton','String', 'Segment (S)',...
    'Position',[160,700,60,20],...
    'Callback',@hsegment_Callback);

hbwplus = uicontrol('Style','pushbutton','String', 'BW+',...
    'Position',[160,730,30,20],...
    'Callback',@hbwplus_Callback);

hbwminus = uicontrol('Style','pushbutton','String', 'BW-',...
    'Position',[190,730,30,20],...
    'Callback',@hbwminus_Callback);

hshow = uicontrol('Style','togglebutton','String', 'Show (X)',...
    'Position',[160,670,60,20],...
    'Callback',@hshow_Callback);

hremove = uicontrol('Style','pushbutton','String', 'Remove (R)',...
    'Position',[230,670,60,20],...
    'Callback',@hremove_Callback);

hadd = uicontrol('Style','togglebutton','String', 'Add (A)',...
    'Position',[230,700,60,20],...
    'Callback',@hadd_Callback);

hautosegfw = uicontrol('Style','pushbutton','String', 'Seg FW',...
    'Position',[300,670,60,20],...
    'Callback',@hautosegfw_Callback);

hautosegbk = uicontrol('Style','pushbutton','String', 'Seg BK',...
    'Position',[300,700,60,20],...
    'Callback',@hautosegbk_Callback);

hautosegbkf = uicontrol('Style','edit','String', '1',...
    'Position',[350,700,30,20],...
    'Callback',@hautosegbkf_Callback);

hautosegfwf = uicontrol('Style','edit','String', '150',...
    'Position',[350,670,30,20],...
    'Callback',@hautosegfwf_Callback);

hsavebw = uicontrol('Style','pushbutton','String', 'Save Mask',...
    'Position',[400,670,100,40],...
    'Callback',@hsavebw_Callback);

hloadbw = uicontrol('Style','pushbutton','String', 'Load Mask',...
    'Position',[520,670,100,40],...
    'Callback',@hloadbw_Callback);

htrack = uicontrol('Style','pushbutton','String', 'Track',...
    'Position',[640,670,100,40],...
    'Callback',@htrack_Callback);

%% Create axes and display image
ImgOp = zeros(100,100);
sI = size(ImgOp);
panelI_X1 = 40;
panelI_X2 = panelI_X1+sI(2);
panelI_Y1 = 40;
panelI_Y2 = panelI_Y1+sI(1);
ha = axes('Units','Pixels','Position',[panelI_X1,panelI_Y1,panelI_X2,panelI_Y2]);
imshow(ImgOp);

hbw = axes('Units','Pixels','Position',[panelI_X1,panelI_Y1+sI(1)+30,panelI_X2,panelI_Y2]);
BWop = zeros(100,100);
imshow(BWop);

% Change units to normalized so components resize automatically.
set([segfigure,ha,hmessages...
    hbw,hpath,hloadbutton,hquitbutton,...
    hframe,hopen,hsegment,hbwplus,hbwminus,hchoice...
    hshow,hadd,hremove,...
    hautosegfw,hautosegbk,hautosegfwf,hautosegbkf,...
    hsavebw,hloadbw,htrack],...
    'Units','normalized');
%% Final settings
% Assign the GUI a name to appear in the window title.
set(segfigure,'Name','Segment NUP V1')
% Move the GUI to the center of the screen.
movegui(segfigure,'center')
% Make the GUI visible.
set(segfigure,'Visible','on')
%% Callbacks
    function click(~,~)
        % Not used...
    end
    function hquitbutton_Callback(~,~)
        close(segfigure);
    end

    function hloadbutton_Callback(~,~)
        [FileName,PathName,FilterIndex] = uigetfile('*.tif','Select the tif file');
        set(hpath,'String',[PathName,FileName]);
        info = imfinfo([PathName,FileName]);
        Nf = length(info);
        Wd = info(1).Width;
        Hg = info(1).Height;
        sI = [Hg,Wd];
        Img = [];
        BW = zeros(Hg,Wd,Nf);
        BWL = zeros(Hg,Wd,Nf);
        Nnuc = zeros(Nf,1);
        for i=1:Nf
            Img = cat(3,Img,imread([PathName,FileName],i));
        end
        sI = size(Img(:,:,1));
        panelI_X1 = 40;
        panelI_X2 = panelI_X1+sI(2);
        panelI_Y1 = 40;
        panelI_Y2 = panelI_Y1+sI(1);
        ha = axes('Units','Pixels','Position',[panelI_X1,panelI_Y1,panelI_X2,panelI_Y2]);
        %imshow(ImgOp);
        
        hbw = axes('Units','Pixels','Position',[panelI_X1,panelI_Y1+sI(1)+30,panelI_X2,panelI_Y2]);
        BW(:,:,1) = im2bw(Img(:,:,1),graythresh(Img));
        %imshow(BWop);
        
        showI(Img(:,:,1),BW(:,:,1));
        set(hframe,'max',Nf,'SliderStep',[1/(Nf-1) 0.2])
        set(hmessages,'String',['Loaded. Frame ',num2str(frame)]);
        set(hautosegfwf,'String',num2str(Nf));
    end

    function hframe_Callback(~,~)
        frame = get(hframe,'Value');
        frame = round(frame);
        set(hframe,'Value',frame);
        showI(Img(:,:,frame),BW(:,:,frame));
        set(hmessages,'String',['Frame ',num2str(frame)]);
    end

    function scroll(src,evnt)
        deltaframe = evnt.VerticalScrollCount;
        frame = frame + deltaframe;
        if frame<1
            frame=1;
        end
        if frame>Nf
            frame = Nf;
        end
        set(hframe,'Value',frame);
        showI(Img(:,:,frame),BW(:,:,frame));
        set(hmessages,'String',['Frame ',num2str(frame)]);
    end

    function hchoice_Callback(~,~)
        choice = get(hchoice,'value');
        if choice ==1
            seg = 0;% soft
        elseif choice ==2
            seg = 1;% hard
        end
    end

    function hopen_Callback(~,~)
        open = get(hopen,'Value');
        open = round(open);
        set(hopen,'Value',open);
        switch seg
            case 0 %soft
                f1 = fspecial('Gaussian', open, open/3);
                background = imopen(Img(:,:,frame),strel('disk',20));
                BW(:,:,frame) = imfilter(Img(:,:,frame)-background,f1,'replicate');
                %f2 = fspecial('Gaussian', open, open/2);
                %df = f1-2.*f2;
                %BW(:,:,frame) = conv2(double(Img(:,:,frame)), df, 'same');
                %BW(:,:,frame) = imcomplement(BW(:,:,frame));
            case 1 %hard
                f1 = fspecial('Gaussian', open, open/3);
                background = imopen(Img(:,:,frame),strel('disk',20));
                BW(:,:,frame) = imfilter(Img(:,:,frame)-background,f1,'replicate');
        end
        showop(BW(:,:,frame));
        set(hmessages,'String',['Frame ',num2str(frame)]);
    end
    
    % Show button - add mask over the image in Panel 2.
    function hshow_Callback(~,~)
        show_state = get(hshow,'Value');
        if show_state == get(hshow,'Max')
            showBW = 1;
        end
        if show_state == get(hshow,'Min')
            showBW = 0;
        end
        showI(Img(:,:,frame),BW(:,:,frame));
    end

    % Remove a mask - Modification on Panel 1.
    function hremove_Callback(~,~)
        [x,y] = getpts(hbw);
        BW2 = bwselect(BW(:,:,frame),x,y, 8);
        BW(:,:,frame) = BW(:,:,frame)-BW2;
        showbw(BW(:,:,frame));
        showI(Img(:,:,frame),BW(:,:,frame));
    end

    function hadd_Callback(~,~)
        r = getrect(ha);
        tempfig = figure('Name','please draw the nucleus');        
        locI = imcrop(Img(:,:,frame),r);
        imshow(locI);
        % Resize the frame 2x for better draw
        framesize=get(gca,'Position');zoomlv=2;
        set(gca,'Position',[0 0 framesize(3:4)*zoomlv]);
        % Draw a nucleus
        addingreg = imfreehand;
        BWtemp = createMask(addingreg);
        close(tempfig)
        r = round(r);
        rs = size(BWtemp);
        r(3) = rs(2)-1;
        r(4) = rs(1)-1;
        BW(r(2):r(2)+r(4),r(1):r(1)+r(3),frame) = BW(r(2):r(2)+r(4),r(1):r(1)+r(3),frame)|BWtemp;
        [idx,ixy] = find(BW(:,:,frame)==2);
        BW(idx,ixy,frame)=ones;
        showbw(BW(:,:,frame));
        showI(Img(:,:,frame),BW(:,:,frame));
    end

    function hsegment_Callback(~,~)
        switch seg
            case 0% soft
                BW(:,:,frame) = segment(BW(:,:,frame));
                BW(:,:,frame) = imclearborder(BW(:,:,frame));
            case 1% hard
                BW(:,:,frame) = segment(BW(:,:,frame));
                %BW(:,:,frame) = imclose(BW(:,:,frame),strel('disk',1));
                %BW(:,:,frame) = imfill(BW(:,:,frame),'holes');
                BW(:,:,frame) = imclearborder(BW(:,:,frame));
        end
        showbw(BW(:,:,frame));
        showI(Img(:,:,frame),BW(:,:,frame));
    end

    % Expanding the size of the mask
    function hbwplus_Callback(~,~)
        if numel(unique(BW(:,:,frame)))<3 % Only work if it is already a mask
            h=[1 1;1 1];
            BW(:,:,frame)=double(conv2(BW(:,:,frame),h,'same')>=1);            
            showbw(BW(:,:,frame));
            showI(Img(:,:,frame),BW(:,:,frame));
        end
    end

    function hbwminus_Callback(~,~)
        if numel(unique(BW(:,:,frame)))<3 % Only work if it is already a mask
            h=[1 1;1 1];
            BW(:,:,frame)=1-double(conv2(1-BW(:,:,frame),h,'same')>=1);
            showbw(BW(:,:,frame));
            showI(Img(:,:,frame),BW(:,:,frame));
        end
    end

    function hautosegfw_Callback(~,~)
        framefin = str2num(get(hautosegfwf,'String'));
        f1 = fspecial('Gaussian', open, open/3);
        switch seg
            case 0 %soft
                for i=frame+1:framefin
                    BW(:,:,i) = imfilter(Img(:,:,i),f1,'replicate');
                    BW(:,:,i) = segment(BW(:,:,i));
                    BW(:,:,i) = imclearborder(BW(:,:,i));
                    set(hmessages,'String',['Frame ',num2str(i)]);
                    showI(Img(:,:,i),BW(:,:,i));
                    drawnow;
                end
            case 1 %hard
                [centers,R,num] = getdisks(BW(:,:,frame));
                BW(:,:,frame) = makedisk(R,centers);

                % select best size
                R = bestsize(Img(:,:,frame),centers,R);
                BW(:,:,frame) = makedisk(R,centers);
                
               % SE = getpatron(Img(:,:,frame),centers,R,num);

%                 % centering disks
%                 centers = centerdisks(centers,R,num,Img(:,:,frame),SE);
                showI(Img(:,:,frame),BW(:,:,frame)); 
                axes(hbw)
                imshowpair(Img(:,:,frame),BW(:,:,frame),'blend')
                
                for i=frame+1:framefin
                    centers = centerdisks(centers,R,num,Img(:,:,i-1),Img(:,:,i));
                    R = bestsize(Img(:,:,i),centers,R);
                    BW(:,:,i) = makedisk(R,centers);
                    showI(Img(:,:,i),BW(:,:,i));
                    axes(hbw)
                    imshowpair(Img(:,:,i),BW(:,:,i),'blend')
                    set(hmessages,'String',['Frame ',num2str(i)]);
                    drawnow;
                end
                frame = framefin;
                set(hframe,'Value',frame);
        end
    end

   function hautosegbk_Callback(~,~)
       frameinit = str2num(get(hautosegbkf,'String'));
       f1 = fspecial('Gaussian', open, open/3);
       
       
       for i=frame-1:-1:frameinit
           BW(:,:,i) = imfilter(Img(:,:,i),f1,'replicate');
           BW(:,:,i) = segment(BW(:,:,i));
           BW(:,:,i) = imclearborder(BW(:,:,i));
           set(hmessages,'String',['Frame ',num2str(i)]);
           showI(Img(:,:,i),BW(:,:,i));
           drawnow;
       end
       frame = frameinit;
       set(hframe,'Value',frame);
   end

    function hautosegfwf_Callback(~,~)
        framefin = str2num(get(hautosegfwf,'String'));
    end

    function hautosegbkf_Callback(~,~)
        frameinit = str2num(get(hautosegbkf,'String'));
    end

    function hsavebw_Callback(~,~)
        k = strfind(FileName,'.');
        save([PathName,FileName(1:k-1),'_Mask.mat'],'BW');
    end

    function hloadbw_Callback(~,~)
        k = strfind(FileName,'.');
        temp = load([PathName,FileName(1:k-1),'_Mask.mat']);
        BW = temp.BW;
        showbw(BW(:,:,frame))
    end

    % Track nuclei identity
    function htrack_Callback(~,~)
        % Count cell:
        count_cell();
        
        % Create a panel to vísualize the process
        tempfig = figure('Name','tracking...');
        locframe = 1;       % Identifier for currently processed frame
        
        

        CC = bwconncomp(BW(:,:,locframe),8);
        stats = regionprops(CC,'Centroid','BoundingBox');
        Nnuc(locframe) = CC.NumObjects;
        countnuc = 1;
        bwtemp = zeros(sI);
        % Allocating the nuclei array
        for j=1:Nnuc(locframe)
            nuc.pixels{countnuc,locframe} = CC.PixelIdxList{j};
            nuc.frames(countnuc,locframe) = 1;
            nuc.positions{countnuc,locframe} = stats(j).Centroid;
            nuc.box{countnuc,locframe} = stats(j).BoundingBox;
            nuc.name{countnuc,locframe} = num2str(countnuc);
            nuc.filiation(countnuc,locframe) = countnuc;
            bwtemp(nuc.pixels{countnuc,locframe}) =  nuc.filiation(countnuc,locframe);
            countnuc = countnuc+1;
        end
        BWL(:,:,locframe) =  bwtemp;
        nucActive{locframe} = find(nuc.frames(:,locframe)==1);
        
        RGB = label2rgb(BWL(:,:,locframe),'lines','k');
        imshow(RGB,[]);
        hold on
        texts = cat(1,{nuc.name{nucActive{locframe},locframe}});
        posal = cat(1,nuc.positions{nucActive{locframe},locframe});
        text(posal(:,1),posal(:,2),texts,'FontSize',14,'FontWeight','bold','Color','w','HorizontalAlignment','center')
        text(-10,-10,strcat('frame ',num2str(locframe)),'FontSize',14,'FontWeight','bold','Color','k')
        
        % Processing following frames
        for i=2:Nf
            locframe = i;
            CC = bwconncomp(BW(:,:,locframe),8);
            stats = regionprops(CC,'Centroid','BoundingBox');
            Nnuc(locframe) = CC.NumObjects;
            
            % Check whether there is division or not:
            if Nnuc(locframe)==Nnuc(locframe-1) % no division
                bwtemp = zeros(sI);
                for j=1:Nnuc(locframe)
                    bwtemp2 = zeros(sI);
                    bwtemp2(CC.PixelIdxList{j}) = 1;
                    bwtemp2 = bwtemp2.*BWL(:,:,locframe-1);
                    idx = find(bwtemp2,1);
                    ident = bwtemp2(idx);
                    % Only process nuclei that exists in previous frame
                    if numel(ident)
                        nuc.pixels{ident,locframe} = CC.PixelIdxList{j};
                        nuc.frames(ident,locframe) = 1;
                        nuc.positions{ident,locframe} = stats(j).Centroid;
                        nuc.box{ident,locframe} = stats(j).BoundingBox;
                        nuc.name{ident,locframe} =  nuc.name{ident,locframe-1};
                        nuc.filiation(ident,locframe) = ident;
                        bwtemp(nuc.pixels{ident,locframe}) =  nuc.filiation(ident,locframe);
                    end
                end
            else % division or gain of nuclei
                hits = zeros(Nnuc(locframe-1),2);
                hits(:,1) = nucActive{locframe-1};
                for j=1:Nnuc(locframe)
                    bwtemp2 = zeros(sI);
                    bwtemp2(CC.PixelIdxList{j}) = 1;
                    bwtemp2 = bwtemp2.*BWL(:,:,locframe-1);
                    idx = find(bwtemp2,1);
                    ident = bwtemp2(idx);
                    hits(ident,2) = hits(ident,2)+1;
                end
                idxdouble = find(hits(:,2)==2);
                bwtemp = zeros(sI);
                dc = 0;
                for j=1:Nnuc(locframe)
                    bwtemp2 = zeros(sI);
                    bwtemp2(CC.PixelIdxList{j}) = 1;
                    bwtemp2 = bwtemp2.*BWL(:,:,locframe-1);
                    idx = find(bwtemp2,1);
                    if isempty(idx)
                        Nnuc(locframe)=Nnuc(locframe)-1
                    else
                        ident = bwtemp2(idx);
                        if (sum(ident==idxdouble)>0)
                            nuc.pixels{countnuc,locframe} = CC.PixelIdxList{j};
                            nuc.frames(countnuc,locframe) = 1;
                            nuc.positions{countnuc,locframe} = stats(j).Centroid;
                            nuc.box{countnuc,locframe} = stats(j).BoundingBox;
                            if dc==0
                                nuc.name{countnuc,locframe} = strcat(nuc.name{ident,locframe-1},'A');
                                dc = 1;
                            elseif dc==1;
                                nuc.name{countnuc,locframe} = strcat(nuc.name{ident,locframe-1},'B');
                                dc = 0;
                            end
                            nuc.filiation(countnuc,locframe) = ident;
                            bwtemp(nuc.pixels{countnuc,locframe}) =  nuc.filiation(countnuc,locframe);
                            countnuc = countnuc+1;
                        else
                            nuc.pixels{ident,locframe} = CC.PixelIdxList{j};
                            nuc.frames(ident,locframe) = 1;
                            nuc.positions{ident,locframe} = stats(j).Centroid;
                            nuc.box{ident,locframe} = stats(j).BoundingBox;
                            nuc.name{ident,locframe} = nuc.name{ident,locframe-1};
                            nuc.filiation(ident,locframe) = ident;
                            bwtemp(nuc.pixels{ident,locframe}) =  nuc.filiation(ident,locframe);
                        end
                    end
                end
            end
            BWL(:,:,locframe) =  bwtemp;
            nucActive{locframe} = find(nuc.frames(:,locframe)==1);
            RGB = label2rgb(BWL(:,:,locframe),'lines','k');
            clf
            imshow(RGB,[]);
            hold on
            texts = cat(1,{nuc.name{nucActive{locframe},locframe}});
            posal = cat(1,nuc.positions{nucActive{locframe},locframe});
            text(posal(:,1),posal(:,2),texts,'FontSize',14,'FontWeight','bold','Color','w','HorizontalAlignment','center')
            text(-10,-10,strcat('frame ',num2str(locframe)),'FontSize',14,'FontWeight','bold','Color','k')
            drawnow
        end

%         for i=frame+1:Nf
%             CC = bwconncomp(BW(:,:,i),8); 
%              Nnuc(i) = CC.NumObjects;
%              for j=1:Nnuc(i)
%                  bwtemp = zeros(sI);
%                  bwtemp(CC2.PixelIdxList{j}) = 1;
%                  overlap = bwtemp.*BW(:,:,i);
%                  idx = find(overlap==1,'first');
%                  
%                  [r,c] = ind2sub(sI,CC2.PixelIdxList{j});
%                  BW2 = bwselect(BW,c,r,n);
%              end
%         L = labelmatrix(CC);;
%         end
    end

%% Functions
    function showI(image,bwi)
        axes(ha);
        if showBW == 0
            imshow(image);
        end
        if showBW == 1
            locsize = size(image);
            ImBW = zeros(locsize(1),locsize(2),3);
            ImBW(:,:,1)=double(image);
            ImBW(:,:,2)=double(image);
            ImBW(:,:,3)=double(image);
            ImBW(:,:,2)=ImBW(:,:,2).*imcomplement(bwi);
            ImBW(:,:,3)=ImBW(:,:,3).*imcomplement(bwi);
            imshow(uint8(ImBW));
        end
        axes(hbw);
        imshow(bwi,[]);
    end

    function showbw(bwimage)
        axes(hbw);
        imshow(bwimage,[]);
    end

    function showop(opimage)
        axes(hbw);
        imshow(opimage,[]);
    end

   function showbwl(bwlimage)
        axes(hbw);
        RGB = label2rgb(bwlimage,'lines','k');
        imshow(RGB,[]);
        hold on
        texts = cat(1,{nuc.name{nucActive{frame},frame}});
        posal = cat(1,nuc.positions{nucActive{frame},frame});
        text(posal(:,1),posal(:,2),texts,'FontSize',14,'FontWeight','bold','Color','w','HorizontalAlignment','center')
   end

    function bwo = segment(bw)
        l = graythresh(uint8(bw));
        bw = im2bw(uint8(bw),l);
        D = bwdist(~bw);
        l = graythresh(uint8(D));
        bw = im2bw(uint8(D),l);
        bw = imdilate(bw,strel('diamond',3));
        bw = imclose(bw,strel('disk',5));
        bw = imfill(bw,'holes');
        bwl = bwconncomp(bw);
        S = regionprops(bwl,'Area');
        Sm = mean([S.Area]);
        bw = bwareaopen(bw,round(Sm/3));
        D = -bwdist(~bw);
        mask = imextendedmin(D,2);
        D2 = imimposemin(D,mask);
        LD = watershed(D2);
        bw(LD == 0) = 0;
        bwo = imopen(bw,strel('diamond',3));
    end

    function bwo = makedisk(R,centers)
        SE = fspecial('disk',R);
        padsize = ceil(R+R/2);
        bwe = zeros(sI+2*padsize);      
        num = size(centers);
        num = num(1);
        for i=1:num
            bwe(centers(i,2)+padsize,centers(i,1)+padsize) = 1;
        end
        bwd = conv2(bwe,SE,'same');
        bwf = im2bw(bwd,graythresh(bwd));
        D = -bwdist(~bwf);
        mask = imextendedmin(D,2);
        D2 = imimposemin(D,mask);
        LD = watershed(D2);
        bwf(LD == 0) = 0;
        bwo = bwf(padsize+1:end-padsize,padsize+1:end-padsize);
    end

    function [centers,R,num] = getdisks(bw)
        [BL,num] = bwlabel(bw);
        %rgb = label2rgb(BL,'jet', 'k', 'shuffle');
        %showbw(rgb);
        s  = regionprops(BL,'centroid','Area','BoundingBox');
        area = cat(1, s.Area);
        centers = cat(1, s.Centroid);
        centers = round(centers);
        ma = mean(area);
        R = ceil(sqrt(ma/3.1415));
    end

    function newcenters = centerdisks(centers,R,num,Im2,Im1)
        bwe = ones(sI);
        % padding image
        padsize = ceil(R+R/2);
        pIm1 = padarray(Im1,[padsize padsize],'symmetric');
        pIm2 = padarray(Im2,[padsize padsize],'symmetric');
        figure()
        imshow(pIm2,[])
        hold on
%          w=hann(3*R+1);
%          [maskr,maskc]=meshgrid(w,w);
%          w=maskr.*maskc;
% w = fspecial('disk',3*R/2);
% w = w./max(max(w));
% size(w)
       f1 = fspecial('Gaussian', open, open/3);
       f2 = fspecial('Gaussian', 15, 5);
        for i=1:num
            rect = [centers(i,1),centers(i,2),3*R,3*R];
            cropi1 = imcrop(pIm1,rect);
            cropi2 = imcrop(pIm2,rect);
            %size(cropi1)


            %cropi1 = double(cropi1).*w;
            %cropi2 = double(cropi2).*w;
            %cropi1 = imfilter(cropi1,f1,'replicate');
            %cropi2 = imfilter(cropi1,f1,'replicate');
            C = normxcorr2(cropi1,cropi2);
            % figure(10)
            % imshow(C,[])
           %C = imfilter(C,f1,'replicate');
%             figure(3)
%             surfc(C)
            [max_c, imax] = max(C(:));
            [ypeak, xpeak] = ind2sub(size(C),imax(1));
            xpeak = xpeak-size(cropi1,2);%+R;
            ypeak = ypeak-size(cropi1,1);%+R;
            
            %newcenters(i,:) = round([xpeak+rect(1)-padsize,ypeak+rect(2)-padsize]);
            newcenters(i,:) = centers(i,:)-[xpeak,ypeak];
            newrect = [newcenters(i,1),newcenters(i,2),3*R,3*R];
            rectangle('Position',newrect,'EdgeColor','w')
           scatter(newcenters(i,1)+padsize,newcenters(i,2)+padsize,10,'MarkerEdgeColor','r','MarkerFaceColor','r')
        end
    end

    function Rb = bestsize(Im,centers,R)
        r = [];
        for radi = R-1:R+1
            bwtemp = makedisk(radi,centers);
            r = [r,corr2(Im,bwtemp)];
        end
        ro = find(r==max(r));
        switch ro
            case 1
                Rb = R-1;
            case 2
                Rb = R;
            case 3
                Rb = R+1;
        end
    end

%     function SE = getpatron(Im,centers,R,num);
%         all = zeros(3*R,3*R,num);
%          padsize = ceil(R+R/2);
%         pIm = padarray(Im,[padsize padsize],'symmetric');
%         for i=1:num
%             rect = [centers(i,1),centers(i,2),3*R,3*R];
%             all(:,:,i) = imcrop(pIm,rect);
%         end
%         SE = mean(double(all),3);
%         f1 = fspecial('disk',R);
%         f1 = padarray(f1,[R R]);
%         SE = SE.*f1;
%         SE = SE./max(max(SE));
%         figure()
%     end

    function propagate_next()
        deltaframe = 1;
        frame = frame + deltaframe;
        if frame<1
            frame=1;
        end
        if frame>Nf
            frame = Nf;
        end
        set(hframe,'Value',frame);
        hopen_Callback();
        hsegment_Callback();
        set(hmessages,'String',['Frame ',num2str(frame)]);
        
    end
    
    function count_cell()
        for locframe=1:Nf
            % Remove small errors in the first frame
            BW(:,:,locframe) = imfill(BW(:,:,locframe),'holes');
            BW(:,:,locframe) = imerode(BW(:,:,locframe),strel('diamond',3));
            BW(:,:,locframe) = bwareaopen(BW(:,:,locframe),50);
            BW(:,:,locframe) = imdilate(BW(:,:,locframe),strel('diamond',3));

            CC = bwconncomp(BW(:,:,locframe),8);
            stats = regionprops(CC,'Centroid','BoundingBox');
            Nnuc(locframe) = CC.NumObjects;
        end
        figure;
        plot(Nnuc);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOTKEY HOTKEY HOTKEY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function keypress(~, eventdata, ~)
        switch eventdata.Key
            case 'a'
                hadd_Callback();
            case 'r'
                hremove_Callback();
            case 's'
                hsegment_Callback();
            case 'x'
                hshow_Callback();
            case 'space'
                propagate_next();
            case '['
                hbwplus_Callback();
            case ']'
                hbwminus_Callback();
            case 'f1'
                show_help();
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELP HELP HELP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function show_help()
        Msg={'Welcome to Nuclei Segmentation of NUP(or Histone v.1.?)',...
            '',...
            'Brief procedure:',...
            '      -Load the maximum projection movies',...
            '      -Select a proper filter size for the image',...
            '      -Segment the frame: ',...
            '            + Use Automatic segmentation',...
            '            + Add or Remove the nuclei mask if needed',...
            '            + Adjust the size of the mask (with BW+, BW-) if needed',...
            '            + Move on to the next frame',...
            '      -Track nuclei and lineages',...
            '',...
            'Hotkeys:',...
            '      S: automatic segmentation of the current frame',...
            '      A: add a mask (using panel 2)',...
            '      R: remove a mask (using panel 1)',...
            '      [: increase mask size by 1',...
            '      ]: decrease mask size by 1',...
            '',...
            'Hints',...
            '      Always save the mask before Track'...
            '      Backup the mask often (?)'...
            };
        msgbox(Msg,'Help','help')
    end
end
