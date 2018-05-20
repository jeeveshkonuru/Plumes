function varargout = CurrentParams(varargin)
% CURRENTPARAMS MATLAB code for CurrentParams.fig
%      CURRENTPARAMS, by itself, creates a new CURRENTPARAMS or raises the existing
%      singleton*.
%
%      H = CURRENTPARAMS returns the handle to a new CURRENTPARAMS or the handle to
%      the existing singleton*.
%
%      CURRENTPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CURRENTPARAMS.M with the given input arguments.
%
%      CURRENTPARAMS('Property','Value',...) creates a new CURRENTPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CurrentParams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CurrentParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CurrentParams

% Last Modified by GUIDE v2.5 27-Apr-2017 23:20:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CurrentParams_OpeningFcn, ...
                   'gui_OutputFcn',  @CurrentParams_OutputFcn, ...
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


% --- Executes just before CurrentParams is made visible.
function CurrentParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CurrentParams (see VARARGIN)

% Choose default command line output for CurrentParams

h = findobj('Tag','PlumesGUI');
if ~isempty(h)
    Plumes_data = guidata(h);
    handles.sout = Plumes_data.sout;
    handles.plumeprofile = Plumes_data.plumeprofile;
else
    h = findobj('Tag','ParamEntry');
    ParamEntry_data = guidata(h);
    handles.sout = ParamEntry_data.sout;
    handles.plumeprofile = ParamEntry_data.plumeprofile;
end

%set title
h = findobj('Tag','text116');
if handles.plumeprofile == 1
   set(h,'String','Plume Profile');
else
    set(h,'String','Parameter Sweep');
end    

h1 = findobj('Tag','text31');
set(h1,'String',num2str(handles.sout.rho0));
h1 = findobj('Tag','text32');
set(h1,'String',num2str(handles.sout.g));

h1 = findobj('Tag','text33');
if length(handles.sout.entr) > 1
    set(h1,'String',strcat(num2str(handles.sout.entr(1)),{':'},num2str(handles.sout.entr(end))));
else   
    set(h1,'String',num2str(handles.sout.entr));
end

h1 = findobj('Tag','text34');
if length(handles.sout.flowrate) > 1
    set(h1,'String',strcat(num2str(handles.sout.flowrate(1)),{':'},num2str(handles.sout.flowrate(end))));
else   
    set(h1,'String',num2str(handles.sout.flowrate));
end    

h1 = findobj('Tag','text35');
if length(handles.sout.zr) > 1
    set(h1,'String',strcat(num2str(handles.sout.zr(1)),{':'},num2str(handles.sout.zr(end))));
else   
    set(h1,'String',num2str(handles.sout.zr));
end

h1 = findobj('Tag','text36');
if length(handles.sout.bi) > 1
    set(h1,'String',strcat(num2str(handles.sout.bi(1)),{':'},num2str(handles.sout.bi(end))));
else   
    set(h1,'String',num2str(handles.sout.bi));
end

h1 = findobj('Tag','text37');
set(h1,'String',num2str(handles.sout.rho_waste));

h1 = findobj('Tag','text38');
if length(handles.sout.rho_f) > 1
    set(h1,'String',strcat(num2str(handles.sout.rho_f(1)),{':'},num2str(handles.sout.rho_f(end))));
else   
    set(h1,'String',num2str(handles.sout.rho_f));
end    

h1 = findobj('Tag','text39');
if length(handles.sout.m_f) > 1
    set(h1,'String',strcat(num2str(handles.sout.m_f(1)),{':'},num2str(handles.sout.m_f(end))));
else   
    set(h1,'String',num2str(handles.sout.m_f));
end

h1 = findobj('Tag','text40');
if length(handles.sout.rho_s) > 1
    set(h1,'String',strcat(num2str(handles.sout.rho_s(1)),{':'},num2str(handles.sout.rho_s(end))));
else   
    set(h1,'String',num2str(handles.sout.rho_s));
end

h1 = findobj('Tag','text41');
if length(handles.sout.m_s) > 1
    set(h1,'String',strcat(num2str(handles.sout.m_s(1)),{':'},num2str(handles.sout.m_s(end))));
else   
    set(h1,'String',num2str(handles.sout.m_s));
end

h1 = findobj('Tag','text42');
if length(handles.sout.rho_w) > 1
    set(h1,'String',strcat(num2str(handles.sout.rho_w(1)),{':'},num2str(handles.sout.rho_w(end))));
else   
    set(h1,'String',num2str(handles.sout.rho_w));
end

h1 = findobj('Tag','text43');
set(h1,'String',num2str(handles.sout.mu));
h1 = findobj('Tag','text48');
set(h1,'String',num2str(handles.sout.kw));
h1 = findobj('Tag','text49');
set(h1,'String',num2str(handles.sout.cp));
h1 = findobj('Tag','text50');
set(h1,'String',num2str(handles.sout.finalfrac));
h1 = findobj('Tag','text51');
set(h1,'String',num2str(handles.sout.fullprof));
h1 = findobj('Tag','text52');
set(h1,'String',handles.sout.profile);
h1 = findobj('Tag','text53');
set(h1,'String',num2str(handles.sout.constprof));
h1 = findobj('Tag','text54');
set(h1,'String','Param Entry');     % edit this part
h1 = findobj('Tag','text55');
set(h1,'String',num2str(handles.sout.c));
h1 = findobj('Tag','text56');
set(h1,'String',num2str(handles.sout.zmax));
h1 = findobj('Tag','text57');
set(h1,'String',num2str(handles.sout.thermal));
h1 = findobj('Tag','text58');
set(h1,'String',num2str(handles.sout.nondim));
h1 = findobj('Tag','text59');
set(h1,'String',num2str(handles.sout.climatology));
h1 = findobj('Tag','text60');
set(h1,'String',num2str(handles.sout.surfacetemp));
h1 = findobj('Tag','text60');
set(h1,'String',num2str(handles.sout.surfacetemp));
h1 = findobj('Tag','ent');
set(h1,'String',num2str(handles.sout.entr_num));
h1 = findobj('Tag','horiz');
set(h1,'String',num2str(handles.sout.entrb));
h1 = findobj('Tag','para');
set(h1,'String',num2str(handles.sout.U_para));
h1 = findobj('Tag','parb');
set(h1,'String',num2str(handles.sout.U_parb));

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CurrentParams wait for user response (see UIRESUME)
% uiwait(handles.CurrentParams);


% --- Outputs from this function are returned to the command line.
function varargout = CurrentParams_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close GUI]
handles;
close(handles.CurrentParams)
