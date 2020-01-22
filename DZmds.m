function varargout = DZmds(varargin)

% Last Modified 3-Dec-2016 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DZmds_OpeningFcn, ...
                   'gui_OutputFcn',  @DZmds_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DZmds is made visible.
function DZmds_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DZmds (see VARARGIN)
%imshow('uhlogo_red.png', 'Parent', handles.axes1);
% Choose default command line output for DZmds
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DZmds wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DZmds_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%-----------Section to select data file--------%%%%%%%%%%%%%%%
% --- Executes on button press in Browser.
function Browser_Callback(hObject, eventdata, handles)
% hObject    handle to Browser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*'},'File Selector');
fullpathname = strcat(pathname, filename);
text1 = fileread(fullpathname);
set(handles.Filepath, 'String', fullpathname); %show path name
d1 = [fullpathname];
[numbers, text1, data_tmp] = xlsread(d1);

data = cell2mat(data_tmp(2:end,:));
[dataR,dataC]=size(data);
nsamples=(dataC/2);
set(handles.numsamples, 'String', nsamples);
text1=data_tmp(1,:);


handles.data=data;
handles.nsamples=nsamples;
handles.text1=text1;
guidata(hObject,handles);

%%%%%%%%%%%%---------Find Dimensions Button----------%%%%%%%%%%%%%%
% --- Executes on button press in dimensions.
function dimensions_Callback(hObject, eventdata, handles)
% hObject    handle to dimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;
nsamples=handles.nsamples;
text1=handles.text1;
for i=1:nsamples
    headers(1,i)=text1(1,2*i-1);
end

rad_on=get(handles.sigma,'selectedobject');
switch rad_on
    case handles.sigma1
        %
    case handles.sigma2
        for i=1:nsamples
            data(:,2*i)=data(:,2*i)./2;
        end
    otherwise
    set(handles.edit_radioselect,'string','');
end        

a = 0;
b = 4000;
c = 1;
x = a:c:b;
x = transpose(x);

tin=linspace(1,length(x),length(x));
pdp_out = zeros(length(x),nsamples+1);
pdp_cdf_out = zeros(length(x),nsamples+1);
bandwidth_out=ceil(transpose(1:1:nsamples));

%select distribution
dist_on=get(handles.comparison,'selectedobject');
switch dist_on
    case handles.PDP
        for i = 1:nsamples;
            mi = data(:,i*2-1);
            mi =mi(isfinite(mi(:,1)),:);
            si = data(:,i*2);
            si =si(isfinite(si(:,1)),:);
            pdpi = pdp5(mi, si, a, b, c);
            pdpi = transpose(pdpi);
            pdp_out(:,i+1) = pdpi;
            pdp_out(:,1) = x;
            pdp_cdfi = transpose(pdpi);
            pdp_normi = pdp_cdfi/sum(pdp_cdfi);
            cumsumi = cumsum(pdp_normi);
            pdp_cdf_out(:,i+1) = (cumsumi);
            pdp_cdf_out(:,1) = x;
        end
        pdp_sample = pdp_out(:,2:end);
        assignin('base','pdp_sample',pdp_sample)
        
    case handles.KDE
        pdp_out=zeros(size(x,1),nsamples+1);
        bandwidth = str2num(get(handles.kde_bandwidth,'String'))
        for i = 1:nsamples;
            mi = data(:,i*2-1);
            mi = mi(isfinite(mi(:,1)),:);
            test=isfinite(bandwidth)
            if test==1;
                si = ones(length(mi)).*bandwidth;
                kdeAi = pdp5(mi, si, a, b, c);
                xmesh1i=x;
                bandwidthi=bandwidth;
            else
                [bandwidthi,kdeAi,xmesh1i,cdfi]=kde(mi,length(tin),a,b);
            end
            pdpi=transpose(interp1(xmesh1i, kdeAi, tin));
            %assignin('base','pdpi', pdpi);
            bandwidth_out(i,2) = bandwidthi;
            pdp_out(:,i+1) = pdpi(1:4001,1);
            pdp_out(:,1) = x;
            pdp_cdfi = transpose(pdpi);
            pdp_normi = pdp_cdfi/sum(pdp_cdfi);
            cumsumi = cumsum(pdp_normi);
            pdp_cdf_out(:,i+1) = (cumsumi);
            pdp_cdf_out(:,1) = x;
        end
        pdp_out(4001, 2:end)=0;
        pdp_sample = pdp_out(:,2:end);
        F = max(pdp_out(:,2:nsamples+1));
        F=max(F);
        
        
    otherwise
        %
end
assignin('base','pdp_sample',pdp_sample);

% Select comparison
crit_on=get(handles.criterion,'selectedobject');
switch crit_on
    case handles.metric
        crit='metricstress';
    case handles.metric_squared
        crit='metricsstress';
    otherwise
        crit='stress';
end 

rad_on=get(handles.comparison_type,'selectedobject');
switch rad_on
    case handles.R2
        %calculate Cross-correlation coefficients for all samples
        for j=1:nsamples;
            for i=1:nsamples;
                [R2(i,j)] = r2(pdp_sample(:,j),pdp_sample(:,i));
            end
        end
        diss=ones(size(R2));
        diss=(diss-R2);
        
        i=1
        [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        
        while stress(1,i)>0.05 && i<nsamples+1;
            i=i+1;
            [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        end
    case handles.Likeness
        %calculate Likeness/Mismatch coefficients for all samples
        for j=1:nsamples;
            for i=1:nsamples;
                [Like(i,j)] = sum(abs(pdp_sample(:,j)-pdp_sample(:,i)))/2;
            end
        end
        assignin('base','Like',Like);
        diss=ones(size(Like));
        diss=Like;
        i=1
        [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        while stress(1,i)>0.05 || i<3;
            i=i+1;
            [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        end   
    case handles.KSv
        [dataR,dataC]=size(data);
       
        for j=1:nsamples;
            for i=1:nsamples;
                ages1=data(:,2*j-1);
                ages2=data(:,2*i-1);
                [decision(i,j), pKup(i,j), dKup(i,j)] = kstest2(ages1,ages2);
            end
        end
        diss=dKup;
        i=1
        [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        while stress(1,i)>0.05 || i<3;
            i=i+1;
            [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        end
   case handles.Kuiper_D
        [dataR,dataC]=size(data);
       
        for j=1:nsamples;
            for i=1:nsamples;
                ages1=data(:,2*j-1);
                ages2=data(:,2*i-1);
                [pKup(i,j), dKup(i,j)] = kuipertest2c_nan(ages1,ages2);
            end
        end
        diss=dKup;
        i=1
        [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        while stress(1,i)>0.05 || i<3;
            i=i+1;
            [XY,stress(1,i),disparities] = mdscale_new(diss,i,'Criterion',crit);
        end
   otherwise
        set(handles.edit_radioselect,'string','');
end
axes(handles.Scree);
plot (1:i,stress);
ylim([0 max(stress)]);
xlim([1 i]);
set(handles.Scree,'Xtick',0:1:i);
xlabel('Dimensions')
ylabel('Stress')
title('Select dimension at "elbow" in Scree plot')

%%%%%%%%%%%%%%---------------Plot Button-----------%%%%%%%%%%%%%%%%%
%Section to compute MDS and plot
% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.MDSplot,'reset');
cla(handles.Shepard,'reset');

dim_in = str2num(get(handles.dimensions_in,'String'));
data=handles.data;
nsamples=handles.nsamples;
text1=handles.text1;
for i=1:nsamples
    headers(1,i)=text1(1,2*i-1);
end

dx=0.01;

% Select one or two sigma uncertainty inputs
rad_on=get(handles.sigma,'selectedobject');
switch rad_on
    case handles.sigma1
        %
    case handles.sigma2
        for i=1:nsamples
            data(:,2*i)=data(:,2*i)./2;
        end
    otherwise
    set(handles.edit_radioselect,'string','');
end 

% Preallocate PDP matrix
a = 0;
b = 4000;
c = 1;
x = a:c:b;
x = transpose(x);
tin=linspace(1,length(x),length(x));
pdp_out = zeros(length(x),nsamples+1);
pdp_cdf_out = zeros(length(x),nsamples+1);
bandwidth_out=ceil(transpose(1:1:nsamples));

%select distribution
dist_on=get(handles.comparison,'selectedobject');
switch dist_on
    case handles.PDP
        for i = 1:nsamples;
            mi = data(:,i*2-1);
            mi =mi(isfinite(mi(:,1)),:);
            si = data(:,i*2);
            si =si(isfinite(si(:,1)),:);
            pdpi = pdp5(mi, si, a, b, c);
            pdpi = transpose(pdpi);
            pdp_out(:,i+1) = pdpi;
            pdp_out(:,1) = x;
            pdp_cdfi = transpose(pdpi);
            pdp_normi = pdp_cdfi/sum(pdp_cdfi);
            cumsumi = cumsum(pdp_normi);
            pdp_cdf_out(:,i+1) = (cumsumi);
            pdp_cdf_out(:,1) = x;
        end
        pdp_sample = pdp_out(:,2:end);
        assignin('base','pdp_sample',pdp_sample)
        
    case handles.KDE
        pdp_out=zeros(size(x,1),nsamples+1);
        bandwidth = str2num(get(handles.kde_bandwidth,'String'))
        for i = 1:nsamples;
            mi = data(:,i*2-1);
            mi = mi(isfinite(mi(:,1)),:);
            test=isfinite(bandwidth)
            if test==1;
                si = ones(length(mi)).*bandwidth;
                kdeAi = pdp5(mi, si, a, b, c);
                xmesh1i=x;
                bandwidthi=bandwidth;
            else
                [bandwidthi,kdeAi,xmesh1i,cdfi]=kde(mi,length(tin),a,b);
            end
            pdpi=transpose(interp1(xmesh1i, kdeAi, tin));
            %assignin('base','pdpi', pdpi);
            bandwidth_out(i,2) = bandwidthi;
            pdp_out(:,i+1) = pdpi(1:4001,1);
            pdp_out(:,1) = x;
            pdp_cdfi = transpose(pdpi);
            pdp_normi = pdp_cdfi/sum(pdp_cdfi);
            cumsumi = cumsum(pdp_normi);
            pdp_cdf_out(:,i+1) = (cumsumi);
            pdp_cdf_out(:,1) = x;
        end
        pdp_out(4001, 2:end)=0;
        pdp_sample = pdp_out(:,2:end);
        F = max(pdp_out(:,2:nsamples+1));
        F=max(F);
        
        
    otherwise
        %
end

% Select comparison
crit_on=get(handles.criterion,'selectedobject');
switch crit_on
    case handles.metric
        crit='metricstress';
    case handles.metric_squared
        crit='metricsstress';
    otherwise
        crit='stress';
end 

rad_on=get(handles.comparison_type,'selectedobject');
switch rad_on
    case handles.R2
        %calculate Cross-correlation coefficients for all samples
        for j=1:nsamples;
            for i=1:nsamples;
                [R2(i,j)] = r2(pdp_sample(:,j),pdp_sample(:,i));
            end
        end
        diss=ones(size(R2));
        diss=(diss-R2);
                assignin('base','diss',diss)
        %calculate and plot MDS plot for cross-correlation coefficient (R2)
        [XY,stress,disparities] = mdscale_new(diss,dim_in,'Criterion',crit);
        handles.XY=XY;
        
        % set figure parameters
        axes(handles.MDSplot);
        set(gca,'FontUnits','points');
        set(gca, 'OuterPosition',[0.449 -.016 0.556 .952]); % [xLeft, yBottom, width, height]
        colours = colormap(jet((nsamples)));
        % plot 2D or 3D MDS
        if dim_in==2
            X=XY(:,1);
            Y=XY(:,2);
            assignin('base','X',X)
            assignin('base','Y',Y)
            %set(gca, 'Position',[0.48 .075 0.5 .799]); % [xLeft, yBottom, width, height]
            hold on; 
            for i=1:nsamples
                scatter(X(i),Y(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black'); 
            end
            [rubbish,i] = sort(R2,1,'descend');
            YX=XY(i(2,:),:);
            YZ=XY(i(3,:),:);
            arrow3(XY,YZ,':l1',0.4);
            arrow3(XY,YX,'k1',0.5);
            text(X+dx,Y+dx,headers(1,:));
            title('2D MDS plot');
            hold off
            
        else 
            X=XY(:,1);
            Y=XY(:,2);
            Z=XY(:,3);
            assignin('base','X',X)
            assignin('base','Y',Y)
            assignin('base','Z',Z)
            hold on; 
            set(gca, 'OuterPosition',[0.4 0.15 0.65 0.65]); % [xLeft, yBottom, width, height]
            %set(gca, 'Position',[0.6 .2 0.3 .5]); % [xLeft, yBottom, width, height]
            for i=1:nsamples
                scatter3(X(i),Y(i),Z(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black');
            end
            
            [rubbish,i] = sort(R2,1,'descend');
            xlim([min(X) max(X)]);
            ylim([min(Y) max(Y)]);
            zlim([min(Z) max(Z)]);
            
            YX=XY(i(2,:),1:3);
            YZ=XY(i(3,:),1:3);
            XZ=XY(i(1,:),1:3);
            arrow3(XZ,YZ,':l1',0.4);
            arrow3(XZ,YX,'k1',0.5);
            text(X+dx,Y+dx,Z+dx,headers(1,:));
            if dim_in>3;
                title('Caution: multidimensional MDS plot shown in 3 dimensions');
            else
                title('3D MDS plot')
            end
            grid on
            hold off
            view([45 45]);
            
        end
        % add legend to MDS plot   
        %legend(headers,'FontSize', 10,'FontName','Arial','Location','Best')% add legend to scatter plot
        %box on; 
        %axis fill;
               
    case handles.Likeness
        %calculate Likeness/Mismatch coefficients for all samples
        for j=1:nsamples;
            for i=1:nsamples;
                [Like(i,j)] = sum(abs(pdp_sample(:,j)-pdp_sample(:,i)))/2;
            end
        end
        assignin('base','Like',Like);
        diss=ones(size(Like));
        diss=Like;
                
        %calculate and plot MDS plot for cross-correlation coefficient (R2)
        [XY,stress,disparities] = mdscale_new(diss,dim_in,'Criterion',crit);
        
        % figure parameters
        axes(handles.MDSplot);
        set(gca,'FontUnits','points');
        set(gca, 'OuterPosition',[0.449 -.016 0.556 .952]); % [xLeft, yBottom, width, height]
        colours = colormap(jet((nsamples)));
        
        % Plot 2D or 3D MDS
        if dim_in==2;
            X=XY(:,1);
            Y=XY(:,2);
            hold on; 
            for i=1:nsamples
                scatter(X(i),Y(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black'); 
            end
            [rubbish,i] = sort(Like,1,'ascend');
            YX=XY(i(2,:),:);
            YZ=XY(i(3,:),:);
            arrow3(XY,YZ,':l1',0.4);
            arrow3(XY,YX,'k1',0.5);
            text(X+dx,Y+dx,headers(1,:));
            title('2D MDS plot');
            hold off            
        else 
            X=XY(:,1);
            Y=XY(:,2);
            Z=XY(:,3);
            hold on; 
            set(gca, 'OuterPosition',[0.4 0.15 0.65 0.65]); % [xLeft, yBottom, width, height]
            for i=1:nsamples
                scatter3(X(i),Y(i),Z(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black');
            end
            
            [rubbish,i] = sort(Like,1,'ascend');
            xlim([min(X) max(X)]);
            ylim([min(Y) max(Y)]);
            zlim([min(Z) max(Z)]);
            YX=XY(i(2,:),1:3);
            YZ=XY(i(3,:),1:3);
            XZ=XY(i(1,:),1:3);
            arrow3(XZ,YZ,':l1',0.4);
            arrow3(XZ,YX,'k1',0.5);
            text(X+dx,Y+dx,Z+dx,headers(1,:));
            if dim_in>3;
                title('Caution: multidimensional MDS plot shown in 3 dimensions');
            else
                title('3D MDS plot')
            end
            grid on
            hold off
            view([45 45]);
            
        end
        
    case handles.KSv
        
        [dataR,dataC]=size(data);
        
       % Calculate D and p values 
        for j=1:nsamples;
            for i=1:nsamples;
                ages1=data(:,2*j-1);
                ages2=data(:,2*i-1);
                [decision(i,j),pKup(i,j), dKup(i,j)] = kstest2(ages1,ages2);
            end
        end
        
        [XY,stress,disparities] = mdscale_new(dKup,dim_in,'Criterion',crit);
        
        [rubbish,sort_index] = sort(dKup,1,'ascend');
        
        % set figure parameters
        axes(handles.MDSplot);
        set(gca,'FontUnits','points');
        set(gca, 'OuterPosition',[0.449 -.016 0.556 .952]); % [xLeft, yBottom, width, height]
        colours = colormap(jet((nsamples)));
        
        % plot 2D and 3D MDS
        if dim_in==2;
            X=XY(:,1);
            Y=XY(:,2);
            hold on; 
            for i=1:nsamples
                scatter(X(i),Y(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black'); 
            end
            
            YX=XY(sort_index(2,:),:);
            YZ=XY(sort_index(3,:),:);
            arrow3(XY,YZ,':l1',0.4);
            arrow3(XY,YX,'k1',0.5);
            text(X+dx,Y+dx,headers(1,:));
            title('2D MDS plot');
            hold off
        else 
            X=XY(:,1);
            Y=XY(:,2);
            Z=XY(:,3);
            hold on; 
            set(gca, 'OuterPosition',[0.4 0.15 0.65 0.65]); % [xLeft, yBottom, width, height]
            for i=1:nsamples
                scatter3(X(i),Y(i),Z(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black'); 
            end
            
            xlim([min(X) max(X)]);
            ylim([min(Y) max(Y)]);
            zlim([min(Z) max(Z)]);
            YX=XY(sort_index(2,:),1:3);
            YZ=XY(sort_index(3,:),1:3);
            XZ=XY(sort_index(1,:),1:3);
            arrow3(XZ,YZ,':l1',0.4);
            arrow3(XZ,YX,'k1',0.5);
            text(X+dx,Y+dx,Z+dx,headers(1,:));
            if dim_in>3;
                title('Caution: multidimensional MDS plot shown in 3 dimensions');
            else
                title('3D MDS plot')
            end
            grid on
            hold off
            view([45 45]);            
        end
        diss=dKup;
        
    case handles.Kuiper_D
        
        [dataR,dataC]=size(data);
       
        % calculate V and p values 
        for j=1:nsamples;
            for i=1:nsamples;
                ages1=data(:,2*j-1);
                ages2=data(:,2*i-1);
                [pKup(i,j), dKup(i,j)] = kuipertest2c_nan(ages1,ages2);
            end
        end
        
        [XY,stress,disparities] = mdscale_new(dKup,dim_in,'Criterion',crit);
        [rubbish,sort_index] = sort(dKup,1,'ascend');
        
        %set figure parameters
        axes(handles.MDSplot);
        set(gca,'FontUnits','points');
        set(gca, 'OuterPosition',[0.449 -.016 0.556 .952]); % [xLeft, yBottom, width, height]
        colours = colormap(jet((nsamples)));
        
        % plot 2D and 3D MDS
        if dim_in==2;
            X=XY(:,1);
            Y=XY(:,2);
            hold on; 
            for i=1:nsamples
                scatter(X(i),Y(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black'); 
            end
            
            YX=XY(sort_index(2,:),:);
            YZ=XY(sort_index(3,:),:);
            arrow3(XY,YZ,':l1',0.4);
            arrow3(XY,YX,'k1',0.5);
            text(X+dx,Y+dx,headers(1,:));
            title('2D MDS plot');
            hold off
        else
            X=XY(:,1);
            Y=XY(:,2);
            Z=XY(:,3);
            hold on; 
            set(gca, 'OuterPosition',[0.4 0.15 0.65 0.65]); % [xLeft, yBottom, width, height]
            for i=1:nsamples
                scatter3(X(i),Y(i),Z(i),60, 'MarkerFaceColor',colours(i, :)...
                    ,'MarkerEdgeColor','black'); 
            end
            
            xlim([min(X) max(X)]);
            ylim([min(Y) max(Y)]);
            zlim([min(Z) max(Z)]);
            YX=XY(sort_index(2,:),1:3);
            YZ=XY(sort_index(3,:),1:3);
            XZ=XY(sort_index(1,:),1:3);
            arrow3(XZ,YZ,':l1',0.4);
            arrow3(XZ,YX,'k1',0.5);
            text(X+dx,Y+dx,Z+dx,headers(1,:));
            if dim_in>3;
                title('Caution: multidimensional MDS plot shown in 3 dimensions');
            else
                title('3D MDS plot')
            end
            grid on
            hold off
            view([45 45]);
            
            
        end
        diss=dKup;
        
    otherwise
        set(handles.edit_radioselect,'string','');
end 
           
%%%%Plot Shepard plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
axes(handles.Shepard);
set(gca,'FontUnits','normalized');
fsize = get(gca,'FontSize')/4;
set(gca,'FontUnits','points');
colours = colormap(jet((nsamples)));
        
        dissimilarities=squareform(diss);
        disp=squareform(disparities);
        distances=pdist(XY);
        high=max(distances);
        
        [dum,ord] = sortrows([disp(:) dissimilarities(:)]);
        
        hold on
        xlim([0 high+0.2]);
        ylim([0 high+0.2]);
        plot(dissimilarities, distances, 'bo')
        plot(dissimilarities(ord),disp(ord),'r.-');
        plot([0 high+0.2],[0 high+0.2],'k:');
        xlabel('Dissimilarities');
        ylabel('Distances/Disparities');
        title(strcat('Shepard Plot, stress = ',num2str(stress)));
        legend(handles.Shepard,{'Distances' 'Disparities' '1:1'}, 'Location','NorthWest');
        hold off
        
handles.diss=diss;
handles.dips=disp;
handles.high=high;
handles.dissimilarities=dissimilarities;
handles.stress=stress;
handles.disparities=disparities;
handles.XY=XY;
handles.headers=headers;
guidata(hObject,handles);
        
        
%%%%%%%%%%%%%%%%---------Clear Plot Button-----------%%%%%%%%%%%%%%%%
%Section to clear plots
% --- Executes on button press in clrplot.
function clrplot_Callback(hObject, eventdata, handles)
% hObject    handle to clrplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.MDSplot,'reset');

cla(handles.Shepard,'reset');



% --- Executes when selected object is changed in comparison_type.
function comparison_type_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in comparison_type 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when selected object is changed in sigma.
function sigma_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in sigma 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%--------------Export Figures Button---------%%%%%%%%%%%
% --- Executes on button press in ExportFig.
function ExportFig_Callback(hObject, eventdata, handles)
% hObject    handle to ExportFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nsamples=handles.nsamples;
text1=handles.text1;
diss=handles.diss;
disp=handles.dips;
high= handles.high;
dissimilarities=handles.dissimilarities;
disparities=handles.disparities;
stress=handles.stress;
XY=handles.XY;
dim_in = str2num(get(handles.dimensions_in,'String'));

for i=1:nsamples
    headers(1,i)=text1(1,2*i-1);
end

f = figure;
MDSplot=handles.MDSplot;
fig_axes=copyobj(MDSplot,f);
if dim_in>2
    set(gca,'Units', 'normalized', 'outerposition', [0.15 0.15 0.65 0.65]);
else
    set(gca,'Units', 'normalized', 'outerposition', [0 0 1 1]);
end


g=figure;
set(gca,'FontUnits','normalized');
fsize = get(gca,'FontSize')/4;
set(gca,'FontUnits','points');
colours = colormap(jet((nsamples)));
        
        dissimilarities=squareform(diss);
        disp=squareform(disparities);
        distances=pdist(XY);
        high=max(distances);
        
        [dum,ord] = sortrows([disp(:) dissimilarities(:)]);
        
        hold on
        xlim([0 high+0.2]);
        ylim([0 high+0.2]);
        plot(dissimilarities, distances, 'bo')
        plot(dissimilarities(ord),disp(ord),'r.-');
        plot([0 high+0.2],[0 high+0.2],'k:');
        xlabel('Dissimilarities');
        ylabel('Distances/Disparities');
        title(strcat('Shepard Plot, stress = ',num2str(stress)));
        legend({'Distances' 'Disparities' '1:1'}, 'Location','NorthWest');
        hold off
        






% --- Executes on button press in ExportMDS.
function ExportMDS_Callback(hObject, eventdata, handles)
% hObject    handle to ExportMDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
XY=handles.XY;
headers=handles.headers;

XY=num2cell(XY)
headers=transpose(headers);
out=horzcat(headers, XY);
[file,path] = uiputfile('*.xls','Save file');
xlswrite([path file], out);

function dimensions_in_Callback(hObject, eventdata, handles)
% hObject    handle to dimensions_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
function dimensions_in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dimensions_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Example_Data_Set.
function Example_Data_Set_Callback(hObject, eventdata, handles)
% hObject    handle to Example_Data_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Example_Data_Set;

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function kde_bandwidth_Callback(hObject, eventdata, handles)
% hObject    handle to kde_bandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kde_bandwidth as text
%        str2double(get(hObject,'String')) returns contents of kde_bandwidth as a double


% --- Executes during object creation, after setting all properties.
function kde_bandwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kde_bandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
