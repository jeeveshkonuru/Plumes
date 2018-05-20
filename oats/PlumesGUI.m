function varargout = PlumesGUI(varargin)
% PLUMESGUI MATLAB code for PlumesGUI.fig
%      PLUMESGUI, by itself, creates a new PLUMESGUI or raises the existing
%      singleton*.
%
%      H = PLUMESGUI returns the handle to a new PLUMESGUI or the handle to
%      the existing singleton*.
%
%      PLUMESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUMESGUI.M with the given input arguments.
%
%      PLUMESGUI('Property','Value',...) creates a new PLUMESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlumesGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlumesGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlumesGUI

% Last Modified by GUIDE v2.5 06-Dec-2017 18:10:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlumesGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PlumesGUI_OutputFcn, ...
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


% --- Executes just before PlumesGUI is made visible.
function PlumesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlumesGUI (see VARARGIN)

% load data from ParamEntry.m

h = findobj('Tag','ParamEntry');

% get handles and other user-defined data associated to Gui1
if ~isempty(h)
    ParamEntry_data = guidata(h);
    handles.S = ParamEntry_data.S;
    handles.z1 = ParamEntry_data.z1;
    handles.drho1 = ParamEntry_data.drho1;
    handles.vel1 = ParamEntry_data.vel1;
    handles.rad1 = ParamEntry_data.rad1;    
    handles.sout = ParamEntry_data.sout;
    handles.rho0 = ParamEntry_data.sout.rho0;    
    handles.Ua = ParamEntry_data.sout.Ua;
    handles.u_para = ParamEntry_data.sout.U_para
    handles.u_parb = ParamEntry_data.sout.U_parb
    handles.rho = ParamEntry_data.rhoa
    %handles.X = ParamEntry_data.S{1}.X;
    handles.X = ParamEntry_data.X
    
    
end

% set old data
handles.sout_old = handles.sout;

% plot current data on axis1

% access popupmenu info
h = findobj('Tag','popupmenu1');
str = get(h, 'String');
val = get(h,'Value');

% Set current data to the selected data set.
switch str{val};
case 'Density' % User selects density.
   handles.xax1 = handles.drho1;
   handles.xunits1 = '(kg/m^3)';
case 'rhoa' % User selects density.
   junk = 1
   handles.xax1 = handles.rho(handles.z1);
   handles.xunits1 = '(kg/m^3)';
case 'Ua' % User selects density.
   junk = 1
   handles.xax1 = handles.Ua(handles.z1,handles.u_para,handles.u_parb);
   handles.xunits1 = '(kg/m^3)';
case 'Velocity' % User selects velocity.
   handles.xax1 = handles.vel1;
   handles.xunits1 = '(m/s)';
case 'Width' % User selects width.
   handles.xax1 = handles.rad1;
   handles.xunits1 = '(m)';
end
axes(handles.axes1)
plot(handles.xax1,handles.z1)
set(gca, 'Ydir', 'reverse')
xlabel(strcat(str{val},{' '},handles.xunits1))
ylabel('z (m)')

% populate edit text boxes with current values
h = findobj('Tag','edit1');
set(h,'String',num2str(handles.sout.zr));

h = findobj('Tag','edit2');
set(h,'String',num2str(handles.sout.bi));

h = findobj('Tag','edit3');
set(h,'String',num2str(handles.sout.flowrate));

% Choose default command line output for PlumesGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlumesGUI wait for user response (see UIRESUME)
% uiwait(handles.PlumesGUI);


% --- Outputs from this function are returned to the command line.
function varargout = PlumesGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% Re-Define Zr
handles.sout.zr = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% Re-Define bi
handles.sout.bi = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

% Re-Define floware
handles.sout.flowrate = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% Determine the selected data set.
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val};
case 'Density' % User selects density.
   handles.xax1 = handles.drho1;
   handles.xunits1 = '(kg/m^3)';
case 'rhoa' % User selects density. 
   class(handles.X)
   class(handles.rad1)
   handles.xax1 = handles.rho(handles.z1);
   %handles.xax1 = handles.rho({1,2,3,4,5,6,7});
   handles.xunits1 = '(kg/m^3)';
case 'x' % User selects x. 
   handles.xax1 = handles.X;
   handles.xunits1 = 'm';
case 'Ua' % User selects density.
   handles.xax1 = handles.Ua(handles.z1,handles.u_para,handles.u_parb);
   handles.xunits1 = '(kg/m^3)';
case 'Velocity' % User selects velocity.
   handles.xax1 = handles.vel1;
   handles.xunits1 = '(m/s)';
case 'Width' % User selects width.
   handles.xax1 = handles.rad1;
   handles.xunits1 = '(m)';
case 'x' % User selects width.
   handles.xax1 = handles.X;
   handles.xunits1 = '(m)';
end

% plot data
axes(handles.axes1)
plot(handles.xax1,handles.z1)
set(gca, 'Ydir', 'reverse')
xlabel(strcat(str{val},{' '},handles.xunits1))
ylabel('z (m)')

% reset loglog button
log = findobj('Tag','togglebutton2');
set(log,'Value',0);

% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set fullprof = 1
handles.sout.fullprof = 1;
% set title to plume profile
handles.plumeprofile = 1;
guidata(hObject, handles);
run CurrentParams.m

% save old copy of handles
handles.sout_old = handles.sout;

% run plumes 
S = plumes_main(handles.sout);
handles.S = S;
handles.z1 = S{1}.Z; 

%Changed by Andrew, Exact density-------------

%handles.drho1 = S{1}.rho0./S{1}.g.*S{1}.F./(S{1}.Q);
handles.drho1 = S{1}.rho0./S{1}.g.*S{1}.F./(S{1}.Q)+S{1}.rhoa(S{1}.Z);

%End changes------------------
handles.vel1 = S{1}.M./(S{1}.Q);
handles.rad1 = S{1}.Q./sqrt(S{1}.M);

% plot
% Determine the selected data set.
h = findobj('Tag','popupmenu1');
str = get(h, 'String');
val = get(h,'Value');
% Set current data to the selected data set.
switch str{val};
case 'Density' % User selects density.
   handles.xax1 = handles.drho1;
   handles.xunits1 = '(kg/m^3)';
case 'rhoa' % User selects density.
   junk = 3
   handles.xax1 = handles.rho(handles.z1);
   handles.xunits1 = '(kg/m^3)';
case 'Ua' % User selects density.
   junk = 3
   handles.xax1 = handles.Ua(handles.z1,handles.u_para,handles.u_parb);
   handles.xunits1 = '(kg/m^3)';
case 'Velocity' % User selects velocity.
   handles.xax1 = handles.vel1;
   handles.xunits1 = '(m/s)';
case 'Width' % User selects width.
   handles.xax1 = handles.rad1;
   handles.xunits1 = '(m)';
end

% plot data
axes(handles.axes1)
plot(handles.xax1,handles.z1)
set(gca, 'Ydir', 'reverse')
xlabel(strcat(str{val},{' '},handles.xunits1))
ylabel('z (m)')

guidata(hObject, handles);



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

% check if sweep has been calculated yet
if isfield(handles,'S2') == 1
    
    % Determine the selected data set.
    str = get(hObject, 'String');
    val = get(hObject,'Value');
    % Set current data to the selected data set.
    switch str{val};
    case 'Final Density' % User selects density.
       handles.yax2 = handles.drho2;
       handles.yunits2 = '(kg/m^3)';
    case 'Final rhoa' % User selects density.
       junk = 4
       handles.yax2 = handles.rho(handles.z1);
       handles.yunits2 = '(kg/m^3)';
    case 'Final Ua' % User selects density.
       junk = 4
       handles.yax2 = handles.Ua(handles.z1,handles.u_para,handles.u_parb);
       handles.yunits2 = '(kg/m^3)';
    case 'Final Width' % User selects width.
       handles.yax2 = handles.rad2;
       handles.yunits2 = '(m)';
    case 'Vertical Extent' % User selects vertical extent.
       handles.yax2 = handles.z2;
       handles.yunits2 = '(m)';
    %11/2
    case 'X' % User selects vertical extent.
       handles.yax2 = handles.X2;
       handles.yunits2 = '(m)';
    
    end

    % plot data
    h = findobj('Tag','popupmenu3');
    strx = get(h,'String');
    valx = get(h,'Value');
    axes(handles.axes2)
    plot(handles.xax2,handles.yax2)
    %set(gca, 'Ydir', 'reverse')
    ylabel(strcat(str{val},{' '},handles.yunits2))
    xlabel(strcat(strx{valx},{' '},handles.xunits2))

    % reset loglog button
    log = findobj('Tag','togglebutton4');
    set(log,'Value',0);
end

% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

% Reset data set so that multiple params don't have multiple values
handles.sout = handles.sout_old;

% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

% set min of desired param
handles.min = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

% set step size of desired varied param
handles.stepsize = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% set max of desired varied param
handles.max = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set full profile to zero
handles.sout.fullprof = 0;

% set title to parameter sweep
handles.plumeprofile = 0;

% set swept parameter to min:stepsize:max
h = findobj('Tag','popupmenu3');
valx = get(h,'Value');
strx = get(h,'String');
switch strx{valx};
case 'Reference Density' 
    handles.sout.rho0 = handles.min:handles.stepsize:handles.max;
case 'Entrainment Coefficient' 
    handles.sout.entr = handles.min:handles.stepsize:handles.max;
case 'Flow Rate' 
    handles.sout.flowrate = handles.min:handles.stepsize:handles.max;
case 'Release Depth' 
    handles.sout.zr = handles.min:handles.stepsize:handles.max;  
case 'Initial Radius' 
    handles.sout.bi = handles.min:handles.stepsize:handles.max;  
case 'Density of Waste' 
    handles.sout.rho_waste = handles.min:handles.stepsize:handles.max; 
case 'Density of Fines' 
    handles.sout.rho_f = handles.min:handles.stepsize:handles.max;
case 'Mass Flow Rate of Fines' 
    handles.sout.m_f = handles.min:handles.stepsize:handles.max;
case 'Density of Sediment' 
    handles.sout.rho_s = handles.min:handles.stepsize:handles.max;
case 'Mass Flow Rate of Sediment' 
    handles.sout.m_s = handles.min:handles.stepsize:handles.max;  
case 'Initial Water Density' 
    handles.sout.rho_w = handles.min:handles.stepsize:handles.max;
case 'ParamA' 
    handles.sout.U_para = handles.min:handles.stepsize:handles.max;
end

% run plumes for param. sweep
S2 = plumes_main(handles.sout);

handles.S2 = S2;

% calculate output params

%Plume Length and release Depth
handles.z2 = [];
for i = 1:length(S2)
    %if S2{i}.Z < 350
        handles.z2 = [handles.z2 S2{i}.Z];
    %end    
end

%11/2
handles.X2 = [];
for i = 1:length(S2)
    %if S2{i}.Z < 350
        handles.X2 = [handles.X2 S2{i}.X];
    %end    
end

%Radius
handles.rad2 = [];
for i = 1:length(S2)
    handles.rad2 = [handles.rad2 S2{i}.Q(end)./sqrt(S2{i}.M(end))];
end

%mean turbulentVelocity
handles.vel2 = [];
for i = 1:length(S2)
    handles.vel2 = [handles.vel2 S2{i}.M(end)./(S2{i}.Q(end))];
end

%Potential Density anomily
handles.drho2 = [];
for i = 1:length(S2)
    %Changed by Andrew, absolute Denisty---------------------
    %handles.drho2 = [handles.drho2 S2{i}.rho0(end)./S2{i}.g.*S2{i}.F(end)./(S2{i}.Q(end))];
    handles.drho2 = [handles.drho2 S2{i}.rho0(end)./S2{i}.g.*S2{i}.F(end)./(S2{i}.Q(end))+S2{i}.rhoa(S2{i}.Z(end))];
    %End changes-------------------
end

% put default plot on axes2 (depth)

% Determine x variable
h = findobj('Tag','popupmenu3');
valx = get(h,'Value');
strx = get(h,'String');
switch strx{valx};
case 'Reference Density' 
    handles.xax2 = []
    for i = 1:length(S2)
        handles.xax2 = [handles.xax S2{i}.rho0];
    end    
    handles.xunits2 = '(kg/m^3)';
case 'Entrainment Coefficient' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.entr];
    end 
    handles.xunits2 = '[/]';
case 'Flow Rate' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.flowrate];
    end 
    handles.xunits2 = '(m^3/s)';
case 'Release Depth' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.zr];
    end 
    handles.xunits2 = '(m)';
case 'Initial Radius' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.bi];
    end 
    handles.xunits2 = '(m)';    
case 'Density of Waste' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.rho_waste];
    end 
    handles.xunits2 = '(kg/m^3)';    
case 'Density of Fines' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.rho_f];
    end 
    handles.xunits2 = '(kg/m^3)';
case 'Mass Flow Rate of Fines' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.m_f];
    end 
    handles.xunits2 = '(kg/s)';
case 'Density of Sediment' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.rho_s];
    end 
    handles.xunits2 = '(kg/m^3)';    
case 'Mass Flow Rate of Sediment' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.m_s];
    end 
    handles.xunits2 = '(kg/s)';    
case 'Initial Water Density' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.rho_w];
    end 
    handles.xunits2 = '(kg/m^3)';
case 'ParamA' 
    handles.xax2 = [];
    for i = 1:length(S2)
        handles.xax2 = [handles.xax2 S2{i}.U_para];
    end 
    handles.xunits2 = '(kg/m^3)';
%Unsure
end

% Determine the y variable
h1 = findobj('Tag','popupmenu2');  
str = get(h1, 'String');
val = get(h1,'Value');
% Set current data to the selected data set.
switch str{val};
case 'Final Density' % User selects density.
   handles.yax2 = handles.drho2;
   handles.yunits2 = '(kg/m^3)';

case 'Final Width' % User selects width.
   handles.yax2 = handles.rad2;
   handles.yunits2 = '(m)';
case 'Vertical Extent' % User selects width.
   handles.yax2 = handles.z2;
   handles.yunits2 = '(m)';
case 'X' % User selects width.
   handles.yax2 = handles.X2;
   handles.yunits2 = '(m)';
end

% plot data
axes(handles.axes2)
plot(handles.xax2,handles.yax2)
ylabel(strcat(str{val},{' '},handles.yunits2))
xlabel(strcat(strx{valx},{' '},handles.xunits2))

% run current params
guidata(hObject, handles);
run CurrentParams.m

%save handles data
guidata(hObject, handles);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% popout axes1
figure;
axes;

% determine if loglog and plot
h = findobj('Tag','togglebutton2');
h1 = findobj('Tag','popupmenu1');
str = get(h1,'String');
val = get(h1,'Value');

if get(h,'Value') == 0
    plot(handles.xax1,handles.z1)
    set(gca, 'Ydir', 'reverse')
    xlabel(strcat(str{val},{' '},handles.xunits1))
    ylabel('z (m)')
else
    loglog(handles.xax1,handles.z1)
    set(gca, 'Ydir', 'reverse')
    xlabel(strcat(str{val},{' '},handles.xunits1))
    ylabel('z (m)')
end
   

guidata(hObject, handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% popout axes2
figure;
axes;

hx = findobj('Tag','popupmenu3');
hy = findobj('Tag','popupmenu2');
strx = get(hx,'String');
valx = get(hx,'Value');
stry = get(hy,'String');
valy = get(hy,'Value');

% determine if loglog and plot
h = findobj('Tag','togglebutton4');

if get(h,'Value') == 0
    plot(handles.xax2,handles.yax2)
    ylabel(strcat(stry{valy},{' '},handles.yunits2))
    xlabel(strcat(strx{valx},{' '},handles.xunits2))
else
    loglog(handles.xax2,handles.yax2)
    ylabel(strcat(stry{valy},{' '},handles.yunits2))
    xlabel(strcat(strx{valx},{' '},handles.xunits2))
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Return to Param Entry GUI
run ParamEntry.m

% Close PlumesGUI 
close(handles.PlumesGUI)
h = findobj('Tag','CurrentParams');
handles1 = guidata(h);
close(handles1.CurrentParams)


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2
log = get(hObject,'Value');
h = findobj('Tag','popupmenu1');
str = get(h,'String');
val = get(h, 'Value');
if log == 1
    % plot data
    axes(handles.axes1)
    loglog(handles.xax1,handles.z1)
    set(gca, 'Ydir', 'reverse')
    xlabel(strcat(str{val},{' '},handles.xunits1))
    ylabel('z (m)')
else
    % plot data
    axes(handles.axes1)
    plot(handles.xax1,handles.z1)
    set(gca, 'Ydir', 'reverse')
    xlabel(strcat(str{val},{' '},handles.xunits1))
    ylabel('z (m)')
end

% Save the handles structure.
guidata(hObject,handles)

% --- Executes on button press in togglebutton4.
function togglebutton4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton4
log = get(hObject,'Value');
hx = findobj('Tag','popupmenu3');
strx = get(hx,'String');
valx = get(hx, 'Value');
hy = findobj('Tag','popupmenu2');
stry = get(hy,'String');
valy = get(hy, 'Value');
if log == 1
    % plot data
    axes(handles.axes2)
    loglog(handles.xax2,handles.yax2)
    ylabel(strcat(stry{valy},{' '},handles.yunits2))
    xlabel(strcat(strx{valx},{' '},handles.xunits2))
else
    % plot data
    axes(handles.axes2)
    plot(handles.xax2,handles.yax2)
    ylabel(strcat(stry{valy},{' '},handles.yunits2))
    xlabel(strcat(strx{valx},{' '},handles.xunits2))
end

% Save the handles structure.
guidata(hObject,handles)



function para_Callback(hObject, eventdata, handles)
% hObject    handle to para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of para as text
%        str2double(get(hObject,'String')) returns contents of para as a double
handles.sout.para = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function para_CreateFcn(hObject, eventdata, handles)
% hObject    handle to para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
