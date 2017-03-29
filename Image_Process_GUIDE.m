function varargout = Image_Process_GUIDE(varargin)
% IMAGE_PROCESS_GUIDE MATLAB code for Image_Process_GUIDE.fig
%      IMAGE_PROCESS_GUIDE, by itself, creates a new IMAGE_PROCESS_GUIDE or raises the existing
%      singleton*.
%
%      H = IMAGE_PROCESS_GUIDE returns the handle to a new IMAGE_PROCESS_GUIDE or the handle to
%      the existing singleton*.
%
%      IMAGE_PROCESS_GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_PROCESS_GUIDE.M with the given input arguments.
%
%      IMAGE_PROCESS_GUIDE('Property','Value',...) creates a new IMAGE_PROCESS_GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Image_Process_GUIDE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Image_Process_GUIDE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Image_Process_GUIDE

% Last Modified by GUIDE v2.5 16-Jul-2016 19:47:49

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Image_Process_GUIDE_OpeningFcn, ...
                   'gui_OutputFcn',  @Image_Process_GUIDE_OutputFcn, ...
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
end

% --- Executes just before Image_Process_GUIDE is made visible.
function Image_Process_GUIDE_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Image_Process_GUIDE (see VARARGIN)

% UIWAIT makes Image_Process_GUIDE wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Choose default command line output for Image_Process_GUIDE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% make input structure available from every function in the GUI
% hObject.UserData.vargin = varargin;
handles.InputStruct = varargin;

% convert inputStruct into a structure format
handles.InputStruct = cell2mat(handles.InputStruct);

% obtain the list of image names to measure the length of the image stack
% this is used to scale the slider
ListOfImageNames = handles.InputStruct(2).ImageLoc;

% set the maximum value of the slider as the total number of images
set(handles.slider1, 'Max', length(ListOfImageNames));
set(handles.slider1, 'Value',1);
set(handles.slider1, 'Min', 1);

% make the swtich callback variable false to start out, and store it in
% handles
handles.SwitchCallBack = false;
% obtain first original image from stack
FirstIM = handles.InputStruct(2).ImageSet(:,:,1);

% display Original image under "Original Image" (second axes)
axes(handles.axes2);
handles.IM2 = imagesc(handles.axes2,FirstIM); colormap(handles.axes2,gray); axis image; axis off;

% obtain first segmented image from stack
FirstSegIM = handles.InputStruct(3).ImageSet(1).SegFile;

% display segmented image under "Segmented Image"
axes(handles.axes1);
handles.IM1 = imagesc(handles.axes1,FirstSegIM); axis image; axis off;

% create a mouse click function to retrieve lication of click, and remove
% object
set(handles.IM1,'ButtonDownFcn', {@ImageClickCallBack,handles,1});

% update the recently changed handles
guidata(hObject,handles);

uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = Image_Process_GUIDE_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles2 = guidata(findobj('Tag','figure1'));
% Get default command line output from handles structure
varargout{1} = handles.InputStruct(1).ImageSet;

delete(handles.figure1);
end

% --- Executes on button press in Done.
function Done_Callback(~,~, handles)
% hObject    handle to Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

ListOfImageNames = handles.InputStruct(2).ImageLoc;

% retrieve image paths of imoverlay images from input structure
SegFileNames = handles.InputStruct(3).ImageLoc;

% save all segmented images to file
for i = 1:length(ListOfImageNames)
    imwrite(handles.InputStruct(3).ImageSet(i).SegFile,SegFileNames{i});
end

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');

% call on uiresume so output function happend
uiresume(handles.figure1);
% 
% % close the GUI
% Image_Process_GUIDE_OutputFcn(handles);

end

% --- Executes on slider movement.
function slider1_Callback(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% if handles.SwitchCallBack is false - Original images should be displayed,
% if handles.SwitchCallBack is true - imfuse images should be displayed.

% call on uiresume to let execution take place
% uiresume(handles.figure1);

% get slider position
Value = round(get(hObject,'Value')); 

% make the slice number available to all functions
handles.SliderVal = Value;

% update slice number on GUI
set(handles.SliceNum,'String',num2str(Value));

% if handles.SwitchCallBack is false - Original images should be displayed,
% if handles.SwitchCallBack is true - imfuse images should be displayed.

if handles.SwitchCallBack == false
    
    % obtain first original image from stack
    origIM = handles.InputStruct(2).ImageSet(:,:,Value);

    % display Original image under "Original Image" (second axes)
    axes(handles.axes2);
    handles.IM2 = imagesc(handles.axes2,origIM); colormap(handles.axes2,gray); axis image; axis off;
    
elseif handles.SwitchCallBack == true
   
    % extract image corresponding to the position of the slider
    TransIM = handles.InputStruct(4).ImageSet(:,:,Value);
    
    % display Original image under "Transparent Segmented Image" (second axes)
    axes(handles.axes2);
    handles.IM2 = imagesc(handles.axes2,TransIM); axis image; axis off
end


% display Segmented image under "Processed Image"
Image = handles.InputStruct(3).ImageSet(Value).SegFile;

axes(handles.axes1);
handles.IM1 = imagesc(handles.axes1,Image); axis image; axis off;

set(handles.IM1,'ButtonDownFcn', {@ImageClickCallBack,handles,Value});

guidata(hObject,handles);

% uiwait(handles.figure1);
end

function ImageClickCallBack(~,eventData, handles, Value)

% call on uiresume to let execution take place
% uiresume(handles.figure1);

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

handles_lastUpdated = guidata(handles.axes1);
coordinates = get(handles_lastUpdated.axes1,'CurrentPoint'); 
x_coord = round(coordinates(1,1));
y_coord = round(coordinates(1,2));
% disp(['coordinates: ', num2str([x_coord y_coord])]);
% disp(['Slice: ',num2str(Value)]);

LabeldIM = handles.InputStruct(1).ImageSet;

objVal = LabeldIM(y_coord-6,x_coord-3,Value);

% disp(['Object Label: ', num2str(objVal)]);

if objVal ~= 0
    % store previous state
    handles.OldLabeldIM = LabeldIM;
    
    % 3D matrix of original images
    OrigMat = handles.InputStruct(2).ImageSet;
    
    % identify the slices where the object existed
    LinInds = find(LabeldIM == objVal); % finds linear indices
    
    [~, ~, dim] = ind2sub(size(LabeldIM),LinInds); % convert linear indices to subscripts (row, col, dim)
    
    objSliceRange = [min(dim) max(dim)];
    
    handles.ErasedObjRange = objSliceRange;
    
    % disp(['Slice Range of Object: ', num2str(objSliceRange)]);
    
    % set all pixels of selected object to 0 (black)
    LabeldIM(LabeldIM == objVal) = 0;
    
    % store new LabeldIM in handles
    handles.InputStruct(1).ImageSet = LabeldIM;
    
    % update handles structure across GUI
    guidata(handles.figure1,handles);
    
    % convert back to binary (threshold of 0.1 --> convert all numbered objects
    % to 1)
    BIM = LabeldIM>0.1;
    % re-draw segmented image superimposed over original
    for i = objSliceRange(1):1:objSliceRange(2)
        
        % use imoverlay to super impose segmented image as red over
        % original
        handles.InputStruct(3).ImageSet(i).SegFile = imoverlay(OrigMat(:,:,i),BIM(:,:,i),[1 0 0]);
        
        % use imfuse to super impose segmented image as transparent over
        % original
        handles.InputStruct(4).ImageSet(:,:,i) = imfuse(OrigMat(:,:,i),BIM(:,:,i),'diff');
        
    end
    
    slider1_Callback(findobj('Tag','slider1'),eventData, handles);
    
%     % extract image corresponding to the position of the slider
%     TransIM = handles.InputStruct(3).ImageSet(:,:,Value);
%     
%     % display Original image under "Transparent Segmented Image" (second axes)
%     axes(handles.axes2);
%     handles.IM2 = imagesc(handles.axes2,TransIM); axis image; axis off

else
    disp('No Object Detected');
end

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');

% uiwait(handles.figure1);
end
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on button press in imtool1.
function imtool1_Callback(~, ~, handles)
% hObject    handle to imtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the original image currently displayed
IM = getimage(handles.axes1);

% use the imtoo of the image processing toolbox
imtool(IM);
end

% --- Executes on button press in imtool2.
function imtool2_Callback(~, ~, handles)
% hObject    handle to imtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IM = getimage(handles.axes2);
imtool(IM);
end


% --- Executes on button press in isosurface_render.
function isosurface_render_Callback(~, ~, handles)
% hObject    handle to isosurface_render (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Display a 3D rendering of the image stack using the isosurface command

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% update handles structure across GUI
Newhandles = guidata(handles.figure1);

L3 = Newhandles.InputStruct(1).ImageSet;

figure('Name','3D Lacunae Rendering'); hold on;

% the isovalue separates stuff greater then and less than that vlaue. I
% believe isosruface uses it for a gradient.
fv = isosurface(L3,0.5);

% note the color of the faces is given in RGB, can be changed per user
% preference. currently displaying some shade of green
patch(fv,'FaceColor',[61/255 163/255 47/255], 'EdgeColor','none'); hold on;

% this is a shortcut for view(-37.5, 30) which defines the view point in
% polar coordinates from 0 view(azimuith, elevation)
view(3);

% lighiting
camlight;
lighting gouraud; % interpolates lighting

% perspective
camproj perspective; % parallel lines drawn from perspective. more closely correlates with human vision.
                     % displays parallel lines as converging

VoxSize = handles.InputStruct(1).Properties;

% makes the visualization of the data true to the image dimensions
% Ex: if each voxel size is of size (0.447,0.447, 1) microns,
% daspect([0.447, 0.447, 1] will make 0.447 units in x
daspect([1/VoxSize(1) 1/VoxSize(2) 1/VoxSize(3)]);
 
% transparency
% alpha (0.6) % this option allows to see the surface transparent. useful
              % when encountering objects within objects and such
hold off

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Executes on button press in Switch_View.
function Switch_View_Callback(hObject, ~, handles)
% hObject    handle to Switch_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function is supposed to switch the Original images displayed to be
% the images created with imfuse, so the user may see a trace of the
% segmented image, transparenly layed over the original image.

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% this is useful because it is hard to evaluate if the red in the segmented
% image goes beyond the lacuna trace

% if handles.SwitchCallBack is false - Original images should be displayed,
% if handles.SwitchCallBack is true - imfuse images should be displayed.

% recal when this function executes, the variable changes, so when in
% inside this function, if the switchcallback is false, we change it to
% true, and display the imfuse images. if it is true, we change it to false
% and display original images.

% display imfuse images 
if handles.SwitchCallBack == false
   
    % extract image corresponding to the position of the slider
    Image = handles.InputStruct(4).ImageSet(:,:,handles.SliderVal);
    
    % display Original image under "Transparent Segmented Image" (second axes)
    axes(handles.axes2);
    handles.IM2 = imagesc(handles.axes2,Image); axis image; axis off
    
    % Change the title from "Original Image" to "Transparent Segmented Image"
    set(handles.ax2_title, 'String', 'Transparent Segmented Image');
    
    handles.SwitchCallBack = true;
    
elseif handles.SwitchCallBack == true
    
    % extract image corresponding to the position of the slider
    Image = handles.InputStruct(2).ImageSet(:,:,handles.SliderVal);

    % display Original image under "Original Image" (second axes)
    axes(handles.axes2);
    handles.IM2 = imagesc(handles.axes2,Image); colormap(handles.axes2,gray); axis image; axis off;
    % Change the title from "Transparent Segmented Image" to "Original Image"
    set(handles.ax2_title, 'String', 'Original Image');
    
    handles.SwitchCallBack = false;
end

guidata(hObject,handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Executes on button press in Back.
function Back_Callback(~, ~, handles)
% hObject    handle to Back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;
% This function should be like Ctrl + Z. It should store the Image set
% before the last button click

if isfield(handles,'OldLabeldIM')
    
    % obtain last LabeldIM
    LabeldIM = handles.OldLabeldIM;
    
    % obtain the range of slices of last erased object
    objSliceRange = handles.ErasedObjRange;
    
    % store the previous LabeldIM in handles
    handles.InputStruct(1).ImageSet = LabeldIM;
    
    % 3D matrix of original images
    OrigMat = handles.InputStruct(2).ImageSet;
    
    % convert back to binary (threshold of 0.1 --> convert all numbered objects
    % to 1)
    BIM = LabeldIM>0.1;
    
    % re-draw segmented image superimposed over original
    for i = objSliceRange(1):1:objSliceRange(2)
        
        % use imoverlay to super impose segmented image as red over
        % original
        handles.InputStruct(3).ImageSet(i).SegFile = imoverlay(OrigMat(:,:,i),BIM(:,:,i),[1 0 0]);
        
        % use imfuse to super impose segmented image as transparent over
        % original
        handles.InputStruct(4).ImageSet(:,:,i) = imfuse(OrigMat(:,:,i),BIM(:,:,i),'diff');
        
    end

    Image = handles.InputStruct(3).ImageSet(handles.SliderVal).SegFile;
    
    % re-plot Image of current slice
    axes(handles.axes1);
    handles.IM1 = imagesc(handles.axes1,Image); axis image; axis off;

    % Since Image is redrawn, must re-apply buttondownfunction 
    set(handles.IM1,'ButtonDownFcn', {@ImageClickCallBack,handles,handles.SliderVal});

else
    disp('No Previous object deleted');
end

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');

end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(~, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% execute the Done call back function - save images and assign output
% Done_Callback('Place Holder',eventData, handles);
Done_Callback(handles);

end
