function varargout = btsvq_gui_04(varargin)

%   
%    BTSVQ_GUI Application M-file for btsvq_gui_04.fig
%    FIG = BTSVQ_GUI_01 launch btsvq_gui_04 GUI.
%    BTSVQ_GUI_04('callback_name', ...) invoke the named callback.
%    
%    Written by Mujahid Sultan; 
%    mujahid.sultan@ryerson.ca; m.sultan@utoronto.ca
% 

warning off
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).

%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% --------------------------------------------------------------------
function varargout = file_open_Callback(h, eventdata, handles, varargin)
try
    [filename, pathname] = uigetfile({'*.xls';'*.txt'; '*.mat';'*.csv'},'Pick a file');
cd (pathname)
end
handles.filename = filename;
handles.pathname = pathname;
guidata(h,handles) % Save the handles structure

%cd ('handles.pathname')
%if strcmp(get(handles.figure1,'SelectionType'),'open')
%figure (handles.figure1);
     [path,name,ext] = fileparts(handles.filename);    
      % save the name portion of the file to save mat file in save_mat_file call back (next two lines) 
      handles.name=name;
      handles.ext = ext;
      guidata(h,handles);

      switch ext
       case '.fig'
		   guide (filename)
	    case '.xls'
         try
            % load wait bar for all cases of the file types 
	    	h1 = waitbar(.25,'Loading XLS File...'); 
   
            %open(filename)
            [data textdata] = xlsread (filename);
            [r c] = size (textdata);
            [m n] = size(data);
            
            % cnames will be equal to the data in the file: get it by adding 1 to the differnce of c and n 
            cnames = textdata (1,((c-n)+1):c)';
            % if (r-m) => 1, it means there is lables column in data, otherwise there are no labels for rows 
            if  (r-m) >= 1
            labels = textdata (((r-m)+1):r,1);
            else
            labels = '';
            end

            % make an artificial loading sequence            
            %h1 = waitbar(.25,'Loading File...');
            for i = .5:.01:1
            waitbar (i,h1)
            end
            close(h1);

            % update handles
            handles.data = data;
            handles.raw_data = data;
            handles.cnames = cnames;
            handles.labels = labels;
            guidata(h,handles);
        catch
			errordlg(lasterr,'File Type Error','modal')
		 end
        case '.csv'
         try
        	prompt={'Total cols';'Text Columns:';'Header Lines:'};
            def={'13';'2';'1'};
            dlgTitle='Labels and specimens';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            dim = str2num (char(answer(1)));
            textcols = str2num (char(answer(2)));
            headlines = str2num (char(answer(3)));
            
            %  initialize Load Wait bar
            h1 = waitbar(0,'Loading CSV File...');
            
            %open(filename)
            data = csvread (filename,headlines,textcols); % all should be numeric values
            [m n] = size(data);

            if textcols <= 0
                labels = '';
                fsts = {''};
            else
                % read first column as labels
                [l] = textread (filename,'%s%*[^\n]','delimiter',',','headerlines',headlines);
                labels = l;
                % prepare the format string
                fsts = {''};
                for i = 1:textcols
                fsts = strcat({'%*s'},fsts);
                end
            end
            
            fstrs = strcat(fsts,{'%s%*[^\n]'});
            D.fstrs(1,:)= fstrs;    
            for i= 2:(dim-textcols)
             % this is to just show the load status bar
             waitbar (2/(dim-textcols),h1);    
              
             fstrs = strcat({'%*s'},fstrs);
             D.fstrs(i) = fstrs;
            end
            %close the wait bar
            close (h1)
            for i = 1:(dim-textcols)
                [s] = textread (filename,[char(D.fstrs(i))],1,'delimiter',',');
                D.cnames(i,:)= s;
            end
            
            % usual check
            if n == length(D.cnames)
                handles.data = data;
                handles.raw_data = data;
                handles.cnames = D.cnames;
                handles.labels = labels;
                guidata(h,handles);
            else
            errordlg(lasterr,'File Type Error','modal')
            return
            end

		  catch
			errordlg(lasterr,'File Type Error','modal')
		  end
        case '.txt'
         try
            prompt={'Total cols';'Text Columns:';'Header Lines:'};
            def={'13';'2';'1'};
            dlgTitle='Labels and specimens';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            dim = str2num (char(answer(1)));
            textcols = str2num (char(answer(2)));
            headlines = str2num (char(answer(3)));
            
            % mD now has the load waitbar
            mD = big_textread (filename,dim,textcols,headlines);
            if isempty (mD.labels)
            handles.labels = '';    
            else
            handles.labels = mD.labels;            
            end
            data = mD.data;
            cnames = mD.cnames;
            [m n] = size(data);

            % usual check
            if n == length(cnames)
                handles.data = data;
                handles.raw_data = data;
                handles.cnames = cnames;
                guidata(h,handles);
            else
                errordlg('File Type Error','modal')
                return
            end
            
         catch
			errordlg(lasterr,'File Type Error','modal')
		 end
        otherwise %case
            try
                load (name);
                %handles.data = data;
                %handles.raw_data = data;
%                 try
%                     handles.labels= labels;
%                     handles.cnames= cnames;
%                 end
                try
                    handles.data = sD.data;
                    handles.cnames = sD.comp_names;
                    handles.labels = sD.labels;    
                    handles.sD = sD;    
                catch
                    errordlg('Data struct not found not found....  ','modal')
                end
                try
                handles.sM = sM;
                catch
                    errordlg('Map struct not found not found....  ','modal')
                end
                try % some times some old mat files might not have ptree or gtree
                    handles.ptree = ptree;
                catch
                    errordlg('ptree not found....  ','modal')
                end
                try % some times some old mat files might not have ptree or gtree
                    handles.ptree_classify = ptree_classify;
                catch
                    errordlg('ptree_classify not found....  ','modal')
                end
                try
                    handles.gtree = gtree;
                catch
                    errordlg('gtree not found....  ','modal')
                end
                
                if isfield(handles,'ptree')==1
                    Handle = findobj('Tag','btsvq_push_button');
                    set(Handle,'Enable','on');
                end     
                if isfield(handles,'ptree_classify')==1
                    Handle = findobj('Tag','btsvq_push_button');
                    set(Handle,'Enable','on');
                end     
                if isfield(handles,'sM')==1
                    Handle = findobj(gcbf,'Tag','som_push_button');
                    set(Handle,'Enable','on');    
                end
                
                handles.normalization_type= 'Previously normalized data';
                Handle = findobj(gcbf,'Tag','som_push_button');
                set(Handle,'Enable','on');
                
                guidata(h,handles);
            catch
			    errordlg(lasterr,'All variables not found in mat file ','modal')
		    end
            guidata(h,handles);
         end
       set(handles.Text2,'String',pathname)  
       set(handles.Text1,'String',filename) % Display current directory    
       Handle = findobj(gcbf,'Tag','som_push_button');
       set(Handle,'Enable','off');
       Handle = findobj(gcbf,'Tag','topology_push_button');
       set(Handle,'Enable','off');
       Handle = findobj(gcbf,'Tag','btsvq_push_button');
       set(Handle,'Enable','off');
       off = [handles.Specimens_radio_button];
       mutual_exclude (off);
       off = [handles.genes_radio_button];
       mutual_exclude (off);
       %axes (handles.axes1);
       %cla;
 
       Handle = findobj(gcbf,'Tag','plot_data_push_button');
       set(Handle,'Enable','on');
       Handle = findobj(gcbf,'Tag','normalization_pop_up_button');
       set(Handle,'Enable','on');
       Handle = findobj(gcbf,'Tag','Specimens_radio_button');
       set(Handle,'Enable','on');
       Handle = findobj(gcbf,'Tag','genes_radio_button');
       set(Handle,'Enable','on');
       
       handles.ploted = 'no';
       %end %if	
       %handles.filename = filename;
       guidata (h,handles);
%p  = handles.data;
%set (gcf,'CurrentAxes',handles.Axes1)
%surf (p);
%    
%set (handles.listbox1,'String',

% --------------------------------------------------------------------
function varargout = plot_data_push_button_Callback(h, eventdata, handles, varargin)
%axes (handles.Axes1);
if strcmp(handles.ploted,'yes')==1
    axes(handles.Axes2)
    p = handles.data;
    surf (p)
else
    Handle = findobj(gcbf,'Tag','Axes2');
    set(Handle,'Visible','on');
    p = handles.data;
    handles.ploted = 'yes';
    %handles.data = p;
    guidata(h,handles);
    set(gcf, 'renderer','zbuffer');
    axes(handles.Axes2)
    %subplot(1,1,1);
    surf (p);
end

% --------------------------------------------------------------------
function varargout = plot_raw_data_Callback(h, eventdata, handles, varargin)
if strcmp(handles.ploted,'yes')==1
    axes(handles.Axes2)
    p = handles.raw_data;
    surf (p)
else
    Handle = findobj(gcbf,'Tag','Axes2');
    set(Handle,'Visible','on');
    handles.ploted = 'yes';
     p =handles.raw_data;
    set(gcf, 'renderer','zbuffer');
    axes(handles.Axes2)
    %subplot(1,1,1);
    surf (p);
end



% --------------------------------------------------------------------
function varargout = normalization_pop_up_button_Callback(h, eventdata, handles, varargin)
val = get (h,'Value');
switch val
case 1
    handles.normalization_type= 'Normalize';
    guidata (h,handles);

case 2
    handles.normalization_type= 'log';
    guidata (h,handles);
case 3
    handles.normalization_type= 'Variance';
    guidata (h,handles);
case 4    
    handles.normalization_type= 'Range';
    guidata (h,handles);
end
guidata (h,handles);

% --------------------------------------------------------------------
function varargout = normalization_ok_push_button_Callback(h, eventdata, handles, varargin)

handles.raw_data = handles.data;
switch handles.normalization_type
case 'Normalize'
    Handle = findobj(gcbf,'Tag','Axes2');
    set(Handle,'Visible','on');
    p  = handles.data;
    handles.ploted = 'yes';
    handles.raw_data = handles.data;
    handles.p = p;
    guidata(h,handles);
    set(gcf, 'renderer','zbuffer');
    axes(handles.Axes2)
    surf (p);
    return
case 'log'
	% Status bar
	h1 = waitbar(.10,'Normalizing Data... please wait..');

    data_normal = som_normalize (handles.data,'log');
    waitbar (.25,h1)
    data_normal = som_normalize (data_normal','log')';
    waitbar (.5,h1)
    handles.data = data_normal;
    guidata (h,handles);
case 'Variance'
	% Status bar
	h1 = waitbar(.10,'Normalizing Data... please wait..');

    data_normal = som_normalize (handles.data,'var');
    waitbar (.25,h1)
    data_normal = som_normalize (data_normal','var')';
    waitbar (.5,h1)
    handles.data = data_normal;
    guidata (h,handles);
case 'Range'
     % Status bar
	h1 = waitbar(.10,'Normalizing Data... please wait..');

    data_normal = som_normalize (handles.data,'range');
    waitbar (.25,h1)
    data_normal = som_normalize (data_normal','range')';
    waitbar (.5,h1)
    handles.data = data_normal;
    guidata (h,handles);
end
%h = waitbar(.25,'Normalizing Data...');
% make an artificial loading sequence
for i = .25:.01:.5
waitbar (i,h1)
end
close(h1)
h1 = waitbar(.75,'Ploting Data...please wait');
axes(handles.Axes2);
p  = handles.data;
surf (p);
for i = .75:.01:1
waitbar (i,h1)
end
%close(h2);
close(h1)
% --------------------------------------------------------------------
function varargout = Specimens_radio_button_Callback(h, eventdata, handles, varargin)
off = [handles.genes_radio_button];
mutual_exclude (off);
% generate the data struct
    if isempty (handles.labels)
    sD = som_data_struct (handles.data, 'comp_names',handles.cnames);
    else 
    sD = som_data_struct (handles.data, 'comp_names',handles.cnames,'labels',handles.labels);
    end
handles.sD = sD;
guidata (h,handles);

Handle = findobj(gcbf,'Tag','partitive_kmeans_push_button');
set(Handle,'Enable','on');
Handle = findobj(gcbf,'Tag','topology_push_button');
set(Handle,'Enable','on');

% --------------------------------------------------------------------
function varargout = genes_radio_button_Callback(h, eventdata, handles, varargin)
off = [handles.Specimens_radio_button];
mutual_exclude (off);
    if isempty (handles.labels)
    sDt = som_data_struct (handles.data','labels',handles.cnames);    
    else
    sDt = som_data_struct (handles.data','labels',handles.cnames,'comp_names',handles.labels);
    end
handles.sD = sDt;
guidata (h,handles);

Handle = findobj(gcbf,'Tag','partitive_kmeans_push_button');
set(Handle,'Enable','on');
Handle = findobj(gcbf,'Tag','topology_push_button');
set(Handle,'Enable','on');

% --------------------------------------------------------------------
function varargout = partitive_kmeans_push_button_Callback(h, eventdata, handles, varargin)
if exist('handles.new_dir') == 0
%make_new_dir 
new = fix(clock);
new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
handles.new_dir = new_dir;
guidata (h,handles);
mkdir (new_dir);
end
% %make_new_dir 
% new = fix(clock);
% new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
% handles.new_dir = new_dir;
% guidata (h,handles);
% mkdir (new_dir);
cd (new_dir);

% if the sDt has been made by the p-kmeans selection buttons
if isfield(handles,'sDt')==1
    sDt = handles.sDt;
else
sDt = som_data_struct (handles.sD.data','labels',handles.sD.comp_names,'comp_names',handles.sD.labels);
end
%axes (handles.axes1);
ptree = p_tree (sDt,'labels','p_tree.txt',1);
handles.sDt = sDt;
handles.ptree = ptree;
guidata (h,handles);
cd ('..');
% check if the SOM has been computed, then activate the BTSVQ
if isfield(handles,'sM')==1
Handle = findobj('Tag','btsvq_push_button');
set(Handle,'Enable','on');
end


% --------------------------------------------------------------------
function varargout = Partitive_kmeans_classify_Callback(h, eventdata, handles, varargin)
if exist('handles.new_dir') == 0
%make_new_dir 
new = fix(clock);
new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
handles.new_dir = new_dir;
guidata (h,handles);
mkdir (new_dir);
end
% %make_new_dir 
% new = fix(clock);
% new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
% handles.new_dir = new_dir;
% guidata (h,handles);
% mkdir (new_dir);
cd (new_dir);

% if the sDt has been made by the p-kmeans selection buttons
if isfield(handles,'sDt')==1
    sDt = handles.sDt;
else
sDt = som_data_struct (handles.sD.data','labels',handles.sD.comp_names,'comp_names',handles.sD.labels);
end
%axes (handles.axes1);
ptree_classify = p_tree_classify(sDt,'labels','p_tree.txt',1);



handles.sDt = sDt;
handles.ptree_classify = ptree_classify;
guidata (h,handles);
cd ('..');
% check if the SOM has been computed, then activate the BTSVQ
if isfield(handles,'sM')==1
Handle = findobj('Tag','btsvq_push_button');
set(Handle,'Enable','on');
end

% --------------------------------------------------------------------
function varargout = Partative_kmeans_Supervised_Callback(h, eventdata, handles, varargin)
if exist('handles.new_dir') == 0
%make_new_dir 
new = fix(clock);
new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
handles.new_dir = new_dir;
guidata (h,handles);
mkdir (new_dir);
end
% %make_new_dir 
% new = fix(clock);
% new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
% handles.new_dir = new_dir;
% guidata (h,handles);
% mkdir (new_dir);
cd (new_dir);

% if the sDt has been made by the p-kmeans selection buttons
if isfield(handles,'sDt')==1
    sDt = handles.sDt;
else
sDt = som_data_struct (handles.sD.data','labels',handles.sD.comp_names,'comp_names',handles.sD.labels);
end

% The whole idea is to cluster a node by k-means, make a mixture model of each child. Find
% posterior probablity of the components for generating the data point from the parent.
% Now we have two scnerios 
% 1) either use all the data points of the parent to a node (which will be in a sense all the genes) or
% 2) only use the codebook of the parent
% 
% There might be some advantages in both-  using one may lead to


%ptree_classify = p_tree_classify_supervised (sDt,'labels','p_tree.txt',1);
[ptree_classify,level0]= p_tree_classify_super_parent (sDt,'labels','p_tree.txt',1);
try
plot_pridictions(ptree_classify, level0)
catch
    % ha ha
end
handles.ptree_classify = ptree_classify;
handles.sDt = sDt;
handles.level0= level0;
guidata (h,handles);
cd ('..');
% check if the SOM has been computed, then activate the BTSVQ
if isfield(handles,'sM')==1
Handle = findobj('Tag','btsvq_push_button');
set(Handle,'Enable','on');
end


% --------------------------------------------------------------------
function varargout = topology_push_button_Callback(h, eventdata, handles, varargin)
	if handles.ext == '.mat'
             ButtonName=questdlg('Do you want to use old topology?', ...
            'Select','Yes','No','No');
        try
        switch ButtonName,
            
        case 'Yes', 
            prompt={'Old map topology:'};
            msize= {num2str(handles.sM.topol.msize)};
            dlgTitle='Want to use old Map topology';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,msize);
            msize = str2num (char(answer(1)));
            handles.map_size = msize;
            
            temp = num2str(msize);
            temp = strcat('[',temp,']');
            set(handles.topology_edit_box,'String',temp);
            handles.map_size = temp;
            guidata (h,handles);
            Handle = findobj('Tag','som_push_button');
            set(Handle,'Enable','on');
            Handle = findobj('Tag','btsvq_push_button');
            set(Handle,'Enable','on');
            return
        case 'No',
			h = waitbar(.25,'SOM Topology...');
			% set the map topology
			sTopol = som_topol_struct('data',handles.sD.data);
			temp = num2str(sTopol.msize);
			temp = strcat('[',temp,']');
			set(handles.topology_edit_box,'String',temp);
			handles.map_size = temp;
			guidata (h,handles);
			% make an artificial loading sequence
			for i = .5:.01:1
			waitbar (i,h)
			end
			close(h);
            Handle = findobj('Tag','som_push_button');
            set(Handle,'Enable','on');
            Handle = findobj('Tag','btsvq_push_button');
            set(Handle,'Enable','on');
            return
        end % switch
    end % try
end
			h = waitbar(.25,'SOM Topology...');
			% set the map topology
			sTopol = som_topol_struct('data',handles.sD.data);
			temp = num2str(sTopol.msize);
			temp = strcat('[',temp,']');
			set(handles.topology_edit_box,'String',temp);
			handles.map_size = temp;
			guidata (h,handles);
			% make an artificial loading sequence
			for i = .5:.01:1
			waitbar (i,h)
			end
			close(h);

%enable the SOM push button
Handle = findobj('Tag','som_push_button');
set(Handle,'Enable','on');

% --------------------------------------------------------------------
function varargout = topology_edit_box_Callback(h, eventdata, handles, varargin)

msize = get (h,'String')
handles.map_size = msize;
guidata (h,handles);

% --------------------------------------------------------------------
function varargout = som_push_button_Callback(h, eventdata, handles, varargin)
if (handles.ext == '.mat') & isfield(handles,'sM')==1 ...
    & length(handles.sD.comp_names) == length(handles.sM.comp_names)
    sM = handles.sM;
    sM_vote= som_autolabel_rank (sM,handles.sD,'add1');
    handles.sM = sM_vote;
    guidata (h,handles);
else
%get the map size from edit box
map_size  =str2num (get(handles.topology_edit_box,'String'));
sM = som_make (handles.sD, 'msize',[map_size]);
% activate this if you want RANKED list
warning off
sM_vote= som_autolabel_rank (sM,handles.sD,'add1');
%sM_vote= som_autolabel (sM,handles.sD,'vote');
handles.sM = sM_vote;
guidata (h,handles);

end
figure (2)
som_show (sM,'umati','all');
som_show_add ('label',sM_vote.labels(:,1),'textsize',6,'textcolor','r');
%title ('Best Matching Label')
figure (3);
hh = som_show (sM,'umati','all');
%som_show_add ('label',sM_vote.labels,'textsize',6,'textcolor','r');
handles.select_label_handle = hh.plane;
guidata (h,handles);

%title ('All Labels')
%make_new_dir 
%new = fix(clock);
%new_dir = strcat (num2str(new(1,1)),'-',num2str(new(1,2)),'_',num2str(new(1,4)),'_',num2str(new(1,5)));
%handles.new_dir = new_dir;
%guidata (h,handles);
%mkdir (new_dir);
%cd (new_dir);

% ........SAVE the UMAT and Level_0

%cd ('..')
figure(4)
% nov 12 2015 make show off
som_show (sM, 'compi','all');
%axes_cplanes = get (gcf, 'Children');



%axes (handles.axes1)
%cla (handles.axes1)
%set(handles.Axes1,'Color',get(0,'defaultUicontrolBackgroundColor'));
%set (axes_cplanes,'Parent',handles.figure1);

%set (gcf,
%cd (handles.new_dir);
%print (gcf, '-dps', '-r200', 'Umat')

%cd ('..');

%set (bb,'Parent',handles.axes1);
%figure (handles.figure1);
%axes (handles.axes1);

% check if the p-kmeans have been done, then activate the BTSVQ
if isfield(handles,'ptree')==1
Handle = findobj('Tag','btsvq_push_button');
set(Handle,'Enable','on');
end
warning on
% --------------------------------------------------------------------
function varargout = btsvq_push_button_Callback(h, eventdata, handles, varargin)

if isfield(handles,'ptree_classify')==1 | isfield(handles,'ptree')==1   
    if isfield(handles,'new_dir') == 0
        %make_new_dir 
        new = fix(clock);
        new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
        handles.new_dir = new_dir;
        guidata (h,handles);
        mkdir(handles.new_dir)    
    end
    
    cd (handles.new_dir);
    %set(gcf, 'renderer','zbuffer');
    %axes(handles.Axes2)
    warning off;
    if isfield(handles,'ptree_classify')==1
        child = p_tree_parent2 (handles.sM,handles.ptree_classify);
    elseif isfield(handles,'ptree')==1
        child = p_tree_parent2 (handles.sM,handles.ptree);
    else
        disp ('p-tree not found');
    end
    
    try
        fid = fopen ('LOG.txt','w');
        fprintf (fid, 'File Name: %s\n Directory: %s\n Normalization: %s\n'... 
            ,handles.filename,pwd,handles.normalization_type);
        fclose(fid);     
    end
    cd ('..');
end

% --------------------------------------------------------------------
function varargout = cancel_push_button_Callback(h, eventdata, handles, varargin)
pos_size = get (handles.figure1,'Position');
dlg_pos = [pos_size(1)+pos_size(3)/5 pos_size(2)+ pos_size(4)/5];
user_response = modaldlg (dlg_pos);

switch user_response
case {'no','cancel'}
return    
case 'yes'
    uiwait
end

% --------------------------------------------------------------------
function varargout = level_edit_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = child_edit_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = load_view_push_button_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = reverse_push_button_Callback(h, eventdata, handles, varargin)

cd ('results');
axes (handles.axes1);
i =1; % initialize the level loop
while 1 % loop till the breaking condition   
%for i = 1 : 8 % 50 levels of Partative tree (may be take it as argument at later stage)     
    try, % exception for the level error
      if isstruct(level(i))
            %fprintf ( 'Level \n',i)
            
            for j = 1:2^i
             try, 
                %if isstruct(level(i).child(j))
                if exist(level(i).child(j).jpg)
                imread (level(i).child(j).jpg)
                    end %if 
            catch, % if the i exceeds the generated level
                fff = 0; % simply take next j value
                %break % break the for loop
            end % try
        end % for j
        i =i+1; 
    else 
        break
    end % if
catch, % if the i exceeds the generated level
    return % break the for loop and return the calling function
end % try
end % for or while

% --------------------------------------------------------------------
function varargout = forward_push_button_Callback(h, eventdata, handles, varargin)


% ------------------------------------------------------------
% Callback for the Prev and Next buttons
% Pass a string argument (str) indicating which object was selected
% ------------------------------------------------------------
function varargout = forward_rev_push_button_Callback(h, eventdata, handles, str)
% Get the index pointer and the addresses
index_level = handles.ptree;

% Depending on whether Prev or Next was clicked change the display
switch str
case '<<'
	% Decrease the index by one
	i = index - 1;	
	% If the index is less then one then set it equal to the index of the 
	% last element in the Addresses array
	if i < 1
		i = length(handles.ptree)-1;
	end
case '>>'
	% Increase the index by one
	i = index + 1;
	
	% If the index is greater than the size of the array then point
	% to the first item in the Addresses array
	if i > length(handles.ptree)-1
		i = 1;
	end	
end
imread = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j)}))

% Get the appropriate data for the index in selected
Current_Name = Addresses(i).Name;
Current_Phone = Addresses(i).Phone;
set(handles.Contact_Name,'string',Current_Name)
set(handles.Contact_Phone,'string',Current_Phone)

% Update the index pointer to reflect the new index
handles.Index = i;
guidata(h,handles)


%-----------------------------------------------------------------------
% FUNCTIONS
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function load_listbox(dir_path,handles)


cd (dir_path) % Change to the specified directory
dir_struct = dir(dir_path); % List contents of directory
[sorted_names,sorted_index] = sortrows({dir_struct.name}'); % Sort names
handles.file_names = sorted_names; % Save the sorted names
handles.is_dir = [dir_struct.isdir]; % Save names of directories
handles.sorted_index = [sorted_index]; % Save sorted index values
guidata(handles.figure1,handles) % Save the handles structure
%set(handles.listbox1,'String',handles.file_names,'Value',1) % Load listbox
set(handles.Text1,'String',pwd) % Display current directory

%-----------------------------------------------------------------------
% to do the mutual function off 
function mutual_exclude (off)
set (off,'Value',0)

% --------------------------------------------------------------------
function varargout = edit_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = edit_copy_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = edit_cut_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = edit_paste_Callback(h, eventdata, handles, varargin)
str =  clipboard('paste');

% --------------------------------------------------------------------
function varargout = edit_select_all_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = file_exit_Callback(h, eventdata, handles, varargin)
clear handles;
close ('all');


% --------------------------------------------------------------------
function varargout = Untitled_1_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = hit_list_Callback(h, eventdata, handles, varargin)

cd (handles.new_dir)
hit_histo_script(handles.sD,handles.sM)
cd ('..')


% % --------------------------------------------------------------------
% function varargout = tw_btsvq_Callback(h, eventdata, handles, varargin)
% 
% cd (handles.new_dir)
% mkdir ('tw_btsvq')
% cd ('tw_btsvq')
% genes = p_tree_btsvq(handles.sD,handles.sM, handles.ptree,'genes_ptree.txt')
% cd('..')
% cd ('..')



% --------------------------------------------------------------------
function varargout = BiTS_CoVQ_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = Hit_List_Callback(h, eventdata, handles, varargin)

cd (handles.new_dir)
hit_histo_script(handles.sD,handles.sM)
cd ('..')

% --------------------------------------------------------------------
function varargout = tw_btsvq_Callback(h, eventdata, handles, varargin)
% try
% newdir = handles.new_dir;
% end
% if exist('newdir','var') ~= 1
% %make_new_dir 
% new = fix(clock);
% new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
% handles.new_dir = new_dir;
% guidata (h,handles);
% mkdir (new_dir);
% end
 
%cd (handles.new_dir)
mkdir ('tw_btsvq2')
cd ('tw_btsvq2')

%gtree = p_tree_btsvq_2(handles.sD,handles.sM, handles.ptree,'genes_ptree.txt')
gtree = p_tree_btsvq_2(handles.sD,handles.sM, handles.ptree,'genes_ptree.txt')


%handles.gtree = gtree;
%guidata(h,handles);
cd ('..')
cd ('..')


% --------------------------------------------------------------------
function varargout = save_mat_file_Callback(h, eventdata, handles, varargin)

sD = handles.sD;
raw_data  = handles.raw_data;

if isfield(handles,'sM')==1 
sM = handles.sM; 
end
save (name ,'raw_data','sD', 'sM')    
if isfield(handles,'ptree')==1 
ptree = handles.ptree;
save (name ,'ptree');%,'-APPEND'); 
elseif isfield(handles,'ptree_classify')==1 
ptree_classify = handles.ptree_classify;
save (name ,'ptree_classify');%,'-APPEND');
end
if isfield(handles,'gtree')==1
    gtree = handles.gtree;
    save (name ,'gtree');%,'-APPEND');
end

% --------------------------------------------------------------------
function varargout = save_as_mat_file_Callback(h, eventdata, handles, varargin)
prompt={'Mat file name? :'};
name= {num2str(handles.name)};
dlgTitle='Save default';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,name);
name = char(answer(1));
sD = handles.sD;
save (name ,'sD');    
if isfield(handles,'raw_data')==1 
raw_data  = handles.raw_data;
save (name ,'raw_data','sD');
end
if isfield(handles,'sM')==1 
sM = handles.sM; 
save (name,'sM','sD');%,'-APPEND')
end


if isfield(handles,'ptree')==1 
ptree = handles.ptree;
save (name ,'ptree');%,'-APPEND'); 
elseif isfield(handles,'ptree_classify')==1 
ptree_classify = handles.ptree_classify;
save (name ,'ptree_classify');%,'-APPEND');
end
if isfield(handles,'gtree')==1
    gtree = handles.gtree;
    save (name ,'gtree');%,'-APPEND');
end

% --------------------------------------------------------------------
function varargout = som_show_Callback(h, eventdata, handles, varargin)
figure(3);
som_show (handles.sM);

% --------------------------------------------------------------------
function varargout = help_menu_Callback(h, eventdata, handles, varargin)
docopt;



% --------------------------------------------------------------------
function varargout = p_kmeans_all_data_Callback(h, eventdata, handles, varargin)
% off = [handles.p_kmeans_codebook];
% mutual_exclude (off);
% generate the data struct
   
    sDt = som_data_struct (handles.data', 'comp_names',handles.cnames,'labels',handles.labels);
    
handles.sDt = sDt;
guidata (h,handles);

Handle = findobj(gcbf,'Tag','partitive_kmeans_push_button');
set(Handle,'Enable','on');
Handle = findobj(gcbf,'Tag','topology_push_button');
set(Handle,'Enable','on');

% --------------------------------------------------------------------
function varargout = p_kmeans_codebook_Callback(h, eventdata, handles, varargin)
% off = [handles.p_kmeans_all_data];
% mutual_exclude (off);

if exist('handles.new_dir') == 0
%make_new_dir 
new = fix(clock);
new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
handles.new_dir = new_dir;
guidata (h,handles);
mkdir (new_dir);
end
% %make_new_dir 
% new = fix(clock);
% new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
% handles.new_dir = new_dir;
% guidata (h,handles);
% mkdir (new_dir);
cd (new_dir);

sM_rank = som_autolabel_rank (handles.sM,handles.sD, 'add1');
sDt = som_data_struct (handles.sM.codebook','labels',handles.sD.comp_names,'comp_names',sM_rank.labels);
handles.sDt = sDt;
guidata (h,handles);

%axes (handles.axes1);
ptree = p_tree (sDt,'labels','p_tree.txt',1);
handles.sDt = sDt;
handles.ptree = ptree;
guidata (h,handles);
cd ('..');
% check if the SOM has been computed, then activate the BTSVQ
if isfield(handles,'sM')==1
Handle = findobj('Tag','btsvq_push_button');
set(Handle,'Enable','on');
end


Handle = findobj(gcbf,'Tag','partitive_kmeans_push_button');
set(Handle,'Enable','on');
Handle = findobj(gcbf,'Tag','topology_push_button');
set(Handle,'Enable','on');

% --------------------------------------------------------------------
function varargout = my_context_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = tw_btsvq_classify_Callback(h, eventdata, handles, varargin)

if isfield(handles,'ptree_classify')==1 | isfield(handles,'ptree')==1   
    if isfield(handles,'new_dir') == 0
        %make_new_dir 
        new = fix(clock);
        new_dir = strcat (date,'_',num2str(new(1,4)),'.',num2str(new(1,5)));
        handles.new_dir = new_dir;
        guidata (h,handles);
        mkdir(handles.new_dir)    
    end
    
    
    cd (handles.new_dir)
    mkdir ('tw_btsvq_classify_rank2')
    cd ('tw_btsvq_classify_rank2')
    
    %gtree = p_tree_btsvq_2(handles.sD,handles.sM, handles.ptree,'genes_ptree.txt')
    %gtree_classify = p_tree_btsvq_classify_rank_oneu(handles.sD,handles.sM, handles.ptree_classify,'genes_ptree_classify.txt')
    %gtree_classify = p_tree_btsvq_classify_rank_one(handles.sD,handles.sM, handles.ptree_classify,'genes_ptree_classify.txt')
    warning off;
    gtree_classify = p_tree_btsvq_classify_rank1up(handles.sD,handles.sM, handles.ptree_classify,'genes_ptree_classify.txt')
    %handles.gtree = gtree;
    warning on;
    guidata(h,handles);
        
    cd ('..')
    cd ('..')
    save probablities gtree_classify
end


% --------------------------------------------------------------------
function varargout = Supervised_btsvq_Callback(h, eventdata, handles, varargin)
sD = handles.sD;
 prompt={'Class 1:';'Class 2:'};
            def={'ADC';'SQCC'};
            dlgTitle='Enter the specimen labels for supervised clustering';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);

            class1 = (char(answer(1)));
            class2 = (char(answer(2)));
            
            
            %$$ following lines are for supervised clustering of lung 11 groups
            %$$ we know the groups in advance, because they have not any fixed words in the labels so we give the indices manually
%             %             % for all 87 samples
%             %             ind1 = [19 20 22 26 37 38 40 53 71 75];
%             %             ind2 = [25 43 48 55 59 67 1 4 5 6 7];
%             %             % for nodes0
%             %             ind1 = [9 10 11 13 17 19 25 28 33 36];
%             %             ind2 = [12 20 23 27 29 30 1 2 3 4 5];
%             % Rec NonRec (VQ1) for 87 samples
%             ind1 = [1 3 4 5 6 7 8 9 13 15 25 41 42 43 48 51 52 55 59 66 67 78 83];
%             ind2 = [19 20 22 26 37 38 40 45 53 58 71 74 75];
%             
%             % good bad (VQ2) for 87 samples
%             %ind1=  [1 4 5 6 7 25 43 48 55 59 67]
%             %ind2 = [19 20 22 26 37 38 40 53 58 71 75]
             
            %$ as ind1 and ind2 are needed fy supervised_btsvq function initialize these by zero, if not comming from above lines
               ind1 = [];
               ind2 = [];


            try
                ButtonName1=questdlg('Do you want to remove house keeping genes','Class1','Yes','No','Yes');
                switch ButtonName1,
                case 'Yes', 
                    prompt={'Labels a:';'Labels b:'};
                    def={'3x SSC';'Arab'};
                    dlgTitle='Enter the specimen labels for supervised clustering';
                    lineNo=1;
                    answer=inputdlg(prompt,dlgTitle,lineNo,def);

                    house1 = strmatch (char(answer(1)),sD.labels);
                    house2 = strmatch (char(answer(2)),sD.labels);
                    sD.labels([house1;house2]) = [];
                    sD.data(([house1;house2]),:) = [];
                    guidata(h,handles);
                    if isempty(ind1) == 1
                        [sD_12, inds_12] = supervised_btsvq (handles.sD, class1, class2)
                    else
                   [sD_12, inds_12] = supervised_btsvq (handles.sD, class1, class2,ind1,ind2)
                   end
               
                case 'No',
                    if isempty(ind1) ==1
                        [sD_12, inds_12] = supervised_btsvq (handles.sD, class1, class2)
                    else
                   [sD_12, inds_12] = supervised_btsvq (handles.sD, class1, class2,ind1,ind2)
                   end

                end % switch
            end % try

% here we need to load the file which contains sD with all genes and the raw_data, as our sD above 
% has been now reduced and also has normalized valuse, so this temp file should have atleast sD             
%cd ('nodes1')
%load temp
%cd ('..')
sD = som_data_struct (sD.data(inds_12,:), 'comp_names',sD.comp_names,'labels',sD.labels(inds_12))

if isstruct (sD)
    handles.data = sD.data;
    handles.cnames = sD.comp_names;
    handles.labels = sD.labels;    
    handles.sD = sD;    
guidata(h,handles);
end


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in sG.
