"""
Data Browser for UV Lidar output

This data browser allows you to open NetCDF files generated as output from the UV lidar, 
plot the data and select areas which you would like to invert. The coordinates of the areas selected can be saved to a 
.lis file for further analysis. 
The script is commented for clarification. 

Built using Python 2.7.11 and Tkinter 8.6

Bugs: 
- The figure only resizes horizontally. 

Author: 
Jessica Brown
MSc Climate Studies, Wageingen University
jessicasbrown92@gmail.com
October 2016

"""
# Changing initial directories
initialdir_new = '/usr/people/brown/Documents/ASCII/UVLidar_data'
initialdir_open = '/usr/people/brown/Documents/ASCII/List_files'
initialdir_save = '/usr/people/brown/Documents/ASCII/List_files'

"""
Imports
"""
import Tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.patches import Rectangle
import tkFileDialog
from netCDF4 import Dataset
import numpy as np
import tkMessageBox
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter, LogFormatterMathtext
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
from decimal import *

number_of_figures = 93
matplotlib.rcParams['figure.max_open_warning'] = number_of_figures # Changes number of figures able to be opened at once

'''
Creating empty lists to store references to objects in 
'''
Filenames = []  # List to store filenames 
frames = []     # List to store frames for each file (parallel)
frames0 = []    # List to store frames for each file (parallel)
frames1 = []
Fig = []        # List to store figures (parallel)
Fig0 = []       # List to store figures (perpendicular)
Fig1 = []
Axes = []       # List to store axes (parallel)
Axes0 = []      # List to store axes (parallel)
Axes1 = []
Im = []         # List to store image (parallel)
Im0 = []        # List to store image (parallel)
Im1 = []
toolbar_frames = []     # List to store toolbar for each file (parallel)
toolbar_frames0 = []    # List to store toolbar for each file (parallel)
toolbar_frames1 = []
List = []                # List to store location of boxes drawn, deleted and remaining for each file
Entry_Frame = []        # List to store frames for changing the colourbar limits
entry_max_index = []    # List to store FileNumberEntry for changing max limit of colourbar
entry_min_index = []    # List to store FileNumberEntry for changing min limit of colourbar
SmallFilenames = []

l_f = LogFormatter(10, labelOnlyBase=False)     # formats colorbar in logspace

class Example(tk.Frame): # Creates root window

    def __init__(self, parent):          
        tk.Frame.__init__(self, parent) # Example class inherits from the Frame container widget
        self.parent = parent            # Save reference to parent widget  
        self.initUI()                   # Delegate creation of user interface to initUI

    '''
    Function to creat main frame which stores all other frames
    '''
    def initUI(self):      

        self.parent.title("Data Browser")       # Title of GUI
        self.grid(sticky=tk.NSEW)               # Make frame visible

        self.var1 = tk.IntVar()                 

        self.var6 = tk.IntVar()
        
        self.var8 = tk.StringVar()          
        self.var8.set('Colourbar lower limit:' )
        
        self.var9 = tk.StringVar()
        self.var9.set('Colourbar upper limit:' )

        self.var10 = tk.StringVar()
        self.var10.set('File name:' )
        
        MainButtons = tk.Frame(self)                                  # Make frame 
        MainButtons.grid(row=0, column=0,  sticky=tk.N+tk.S+tk.E+tk.W)# make frame visible
        self.columnconfigure(0, weight=2)
        self.rowconfigure(0, weight=2)
        
        self.menubar = tk.Menu(self)                # Make reference to menubar
        menu = tk.Menu(self.menubar, tearoff=0)     # Make menubar
        self.menubar.add_cascade(label="File", menu=menu)   # make menubar visible, give name
        menu.add_command(label="New", command=self.AddNewFile, accelerator='Ctrl+N')       # Add command to menubar
        menu.add_command(label="Open..", command=self.OpenSavedFile, accelerator='Ctrl+O')  
        menu.add_command(label="Save", command=self.save, accelerator='Ctrl+S')         
        
        self.bind_all('<Control-n>', self.AddNewFileShortcut)  # Bind shortcut to function
        self.bind_all('<Control-o>', self.OpenSavedFileShortcut) 
        self.bind_all('<Control-s>', self.saveShortcut)

        try:                                            # make menubar visible
            self.master.config(menu=self.menubar)   
        except AttributeError:
            # master is a toplevel window (Python 1.4/Tkinter 1.63)
            self.parent.tk.call(parent, "config", "-menu", self.menubar)
              
        button = tk.Button(MainButtons, text="Show file names", command=self.PrintCases) # Add button to frame
        button.grid(row=0, column=0, sticky=tk.NSEW)                         # Make button visible 

        self.var5 = tk.StringVar()                  # Set variable type (String)
        self.var5.set('Number of files: -' )        # Set variable text 
        label = tk.Label(MainButtons, textvariable = self.var5) # Add label
        label.grid(row=1, column=0, columnspan=1)
        
        delete = tk.Button(MainButtons, text = 'Remove current file', command = self.delete_frame)
        delete.grid(row=0, column=1, sticky=tk.NSEW)
        
        self.var4 = tk.StringVar()
        self.var4.set('Change to file:')
        label = tk.Label(MainButtons, textvariable = self.var4)
        label.grid(row=0, column=2, columnspan=3)#, sticky=tk.NSEW)  
                                   
        minus = tk.Button(MainButtons, text = '-', command = self.Minus)
        minus.grid(row=1, column=2, sticky=tk.NSEW)

        global FileNumberEntry                                  # Making 'FileNumberEntry' a global name
        FileNumberEntry = tk.Entry(MainButtons, width = 10)     # Adding FileNumberEntry to frame
        FileNumberEntry.grid(row=1, column=3, sticky=tk.NSEW)   # Making FileNumberEntry visible
        FileNumberEntry.insert(0, '-')                          # Inserting '-'
        FileNumberEntry.bind('<Return>', self.ChangeFrame)      # Binding <Return> in FileNumberEntry to function 
        
        plus = tk.Button(MainButtons, text = '+', command = self.Plus)
        plus.grid(row=1, column=4, sticky=tk.NSEW)        
        
        label = tk.Label(MainButtons, textvariable = '') # Add blank space
        label.grid(row=1, column=5, columnspan=1) 

        label = tk.Label(MainButtons, textvariable = self.var10) # Add label
        label.grid(row=1, column=6, columnspan=1)
        
        QuitFrame = tk.Frame(self)
        QuitFrame.grid(row=0, column=2, sticky=tk.NSEW)
        
        button = tk.Button(QuitFrame, text="Quit", command=self.quit)
        button.grid(row=0, column=0, sticky=tk.NSEW)        
        
        SelectButtons = tk.Frame(self)
        SelectButtons.grid(row=1, column=2, rowspan=10, sticky=tk.NSEW)

        label = tk.Label(SelectButtons, textvariable = '')
        label.grid(row=0, column=0)
       
        Radio0 = tk.Radiobutton(SelectButtons, text = 'Select areas', variable = self.var1, indicatoron=0, value = 1)
        Radio0.grid(row=1, column=0, sticky=tk.W)
        
        Radio1 = tk.Radiobutton(SelectButtons, text = 'Delete areas', variable = self.var1, indicatoron=0, value = 2)
        Radio1.grid(row=2, column=0, sticky=tk.W)

        label = tk.Label(SelectButtons, textvariable = '')
        label.grid(row=3, column=0)
                
        Radio2 = tk.Radiobutton(SelectButtons, text="Perpendicular", variable = self.var6, indicatoron=0, value = 3)
        Radio2.grid(row=4, column=0, sticky=tk.NSEW)
        Radio2.invoke()
        
        Radio3 = tk.Radiobutton(SelectButtons, text="Parallel", variable = self.var6, indicatoron=0, value = 4)
        Radio3.grid(row=5, column=0, sticky=tk.NSEW)
 
        Radio4 = tk.Radiobutton(SelectButtons, text="Depolarisation", variable = self.var6, indicatoron=0, value = 5)
        Radio4.grid(row=6, column=0, sticky=tk.NSEW)
        
        button = tk.Button(SelectButtons, text="Change view", command=self.view)
        button.grid(row=7, column=0, sticky=tk.W)
    '''
    Function to open previous files based on saved list file
    '''
    def OpenSavedFile(self):
        ftypes = (('List files', '*.lis'), ('All files', '*.*'))          # File types which can be opened  
        dlg = tkFileDialog.Open(filetypes = ftypes, initialdir = initialdir_open) # Tk window for selecting files
        fl = dlg.show()                             # Show Tk window
        if fl is '':
            pass
        else:
			try:
				global oldfilename                       
				oldfilename = fl
				'''
				Setting empty arrays to store variables from list file. 
				Lists are emptied every time a new list file is opened
				'''
				global File_num, File_line, Numbers_line, Numbers_num, Rect_line, Rect_num, Rect_key_list, Diff 
				File_num = [] # Stores index of lines where the file names appear
				File_line = [] # Stores file names
				
				Numbers_line = [] # Stores index of lines where the number of files + boxes appear
				Numbers_num = [] # Stores number of files + boxes
				
				Rect_num = [] # Stores index of lines where boxes appear
				Rect_line = []  # Stores boxes
				Rect_key_list = []
				Diff = [] # Difference between sequential indexes in Rect_num

				with open(oldfilename, 'r') as inFile:  # open list file
					for num, line in enumerate(inFile): # read through file line by line, retrieving index and line
						if '/' in line:                 # Finding lines with filenames
							File_num.append(num)
							File_line.append(line)
						elif len(line) < 10:      # finding lines with number of files/boxes
							Numbers_line.append(line)
							Numbers_num.append(num)
						else:
							Rect_raw = map(float, (map(lambda s: s.strip(','), line.split() )) ) # splitting line into individual floats    
							Rect_line.append(Rect_raw)      # adding to line
							Rect_num.append(num)
					Rect_key_list = [len(Rect_num), 0]      # creating list to store indexes
					for i in range(len(Rect_num)):   
						diff = Rect_num[i]-Rect_num[i-1]    # difference between sequential indexes of boxes
						Diff.append(diff)                   # store the locations
						if Rect_num[i]-Rect_num[i-1] > 1:   
							Rect_key = i                        # finding locations when diff between indexes >1 
							Rect_key_list.insert(1, Rect_key)
							d = {}                               # creating dictionary
							for x in range(len(Rect_key_list)-1):   # In loop, create variable names in dictionary 
								d['Rect{}'.format(x)] = [x]
							for j in range(len(Rect_key_list)-1):
								Rect_key_list_rev = list(reversed(Rect_key_list)) # Reverse list 
								Rect_separate = Rect_line[Rect_key_list_rev[j]:Rect_key_list_rev[j+1]] # Split list into multiple lists containing only the boxes for each file
								d['Rect%s' % j] = Rect_separate # Add separated lists to dictionary
						if max(Diff) <=1:           # If there is only 1 file in lis file
							d = {'Rect0':Rect_line}

					for i in range(len(File_line)):
						full_path = File_line[i]            # Getting full file name path
						partial_path = os.path.split(full_path)[1] # Splitting file name path
						small_path = partial_path[:-1]  # Only getting file name out of path
						
						SmallFilenames.append(small_path)
						Filenames.append(full_path)  # Adding file netcdf file name to list
						if len(Filenames) >= 0:        
							FileNumberEntry.delete(0, tk.END)  # Changing number of file names
							FileNumberEntry.insert(0, '%s' % (len(Filenames)))

						self.var5.set('Number of files: %s' % (len(Filenames)))  # Updating number of open files
						self.var10.set('File name: %s' % small_path)
						self.update_idletasks()

						fig_index = int(FileNumberEntry.get())       # Getting current figure index      
						'''
						Adding sublist to List. 
						One sublist added per file to store boxes.  
						Structure: [[Originally drawn], [Boxes remaining following deletion], [boxes deleted]] 
						Each of these sublists contains a sublist with the coordinates of each box
						'''
						List.append([[], [], []])  

						Toolbar = tk.Frame(self)  # Creating frame to contain matplotlib toolbar
						Toolbar.grid(row=4, column=0,  sticky=tk.W) # Making frame visible
						toolbar_frames.append(Toolbar) # storing location of toolbar frame

						Toolbar1 = tk.Frame(self)
						Toolbar1.grid(row=4, column=0,  sticky=tk.W)
						toolbar_frames1.append(Toolbar1)
						
						Frame=tk.Frame(self)  # Making frame to contain figure (parallel)
						Frame.grid(row=2, column = 0, columnspan = 2, rowspan=2, sticky=tk.NSEW) # Making frame visible
						frames.append(Frame) # Storing location of frame

						Frame1=tk.Frame(self)  #(perpendicular)
						Frame1.grid(row=2, column = 0, columnspan = 2, rowspan=2, sticky=tk.NSEW)
						frames1.append(Frame1)

						Frame0=tk.Frame(self)  #(perpendicular)
						Frame0.grid(row=2, column = 0, columnspan = 2, rowspan=2, sticky=tk.NSEW)
						frames0.append(Frame0)

						Toolbar0 = tk.Frame(self)
						Toolbar0.grid(row=4, column=0,  sticky=tk.W)
						toolbar_frames0.append(Toolbar0)

						nc_file = Dataset(small_path, 'r') # reading net cdf file
						global time_bnds
						global para_beta
						global perp_beta
						global altitude
						global depol
						
						time_bnds = nc_file.variables['time_bnds'][:]  # Getting time data
						altitude = nc_file.variables['altitude'][:]
						para_beta = np.swapaxes((nc_file.variables['para_beta'][:]), 0, 1) # getting para data and swapping axes yo read correctly
						perp_beta = np.swapaxes((nc_file.variables['perp_beta'][:]), 0, 1)
						para_beta[np.isnan(para_beta)]=0 # Changing nan values to 0 
						perp_beta[np.isnan(perp_beta)]=0
						with np.errstate(divide='ignore', invalid='ignore'):
							depol = perp_beta/para_beta
						depol[np.isnan(depol)]=0
						
						fig, ax = plt.subplots()  # Creating figure for parallel figure
						self.canvas = FigureCanvasTkAgg(fig, master=Frame) # Creating tk canvas
						self.canvas.get_tk_widget().grid(row=3, column=0, columnspan=5, rowspan=2, sticky=tk.NSEW) # Making canvas visible
		
						ax.set_xlim([0, 24]) # Setting axes limits
						ax.set_ylim([0, 10])
						self.im  = ax.imshow(para_beta, extent =(time_bnds.min(), time_bnds.max(), altitude.min(), #plotting image
										altitude.max()), aspect = 'auto', cmap = 'jet', interpolation = 'none',
										origin = 'lowest', norm=LogNorm())

						ax.set_xlabel('Time (Decimal hours)') # x axes label 
						ax.set_ylabel('Height (km)')
						ax.grid(True)
						majorLocator    = MultipleLocator(1)
						ax.yaxis.set_major_locator(majorLocator)  
						cax = fig.add_axes([0.89, 0.25, 0.03, 0.5]) # Position of colourbar axes
						fig.colorbar(self.im, label = 'Atten Beta [$km^{-1}\,sr^{-1}$]', cax=cax, format=l_f)#, ticks = lvls) # Making colourbar
						plt.subplots_adjust(left=None, right=0.85, top=None, bottom=None) # adjusting position of figure

						self.im.set_clim(1e-05, 1e-01)

						Rects_get = d['Rect%s' % i]  # Reading boxes from dictionary - for current file only
						Rects_get = np.asarray(Rects_get)  # Changing list to array
						for i in range(len(Rects_get)): # Reading through sublist (now array)
							x0 = Rects_get[i,0] # Retrieving coordinates
							y0 = Rects_get[i,3]
							w  = Rects_get[i,1] - x0
							h  = Rects_get[i,2] - y0
							x1 = Rects_get[i,1]
							y1 = Rects_get[i,2]
							self.rect = Rectangle((x0, y0), w, h) # Setting rectangle coordinates
							ax.add_patch(self.rect) # Adding patch/rectangle
							self.rect.set_picker(5) # Setting picking on 
							self.rect.set_edgecolor('red') 
							self.rect.set_linestyle('solid')
							self.rect.set_facecolor('none')
							self.rect.set_linewidth('5')
							coords =([x0, x1, y1, y0])

							if len(List[fig_index-1][1]) <= 0:  # Appending to List to store boxes
								List[fig_index-1][0].append(coords) # Adding to list which remaining boxes (remaining following deleting boxes)
							else:
								List[fig_index-1][1][0].append(coords) # Else just add to list from which no boxes have been deleted,

						self.canvas.draw()  # Redraw canvas
						Fig.append(fig)     # Store location of figure
						Axes.append(ax)     # Store location of axes
						Im.append(self.im)  # Store location of image

						ax.format_coord = lambda x, y: '' # Stop matplotlib toolbar showing coordinates when mouse is over figure 

						toolbar = NavigationToolbar2TkAgg(self.canvas, Toolbar) # add matplotlib toolbar
						toolbar.update()  
						toolbar.pack(side=tk.TOP) # Make toolbar visible 

						aid = self.canvas.mpl_connect('button_press_event', self.on_press) # Connect function to event on canvas
						bid = self.canvas.mpl_connect('button_release_event', self.on_release)
						cid = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
						fid = self.canvas.mpl_connect('pick_event', self.patch_picker)

						fig1, ax1 = plt.subplots() # Draw depol 
						self.canvas1 = FigureCanvasTkAgg(fig1, master=Frame1)
						self.canvas1.get_tk_widget().grid(row=3, column=0, columnspan=5, rowspan=2, sticky=tk.NSEW)
						Frame1.columnconfigure(0, weight=1)
						Frame1.rowconfigure(3, weight=1)

						ax1.set_xlim([0, 24])
						ax1.set_ylim([0, 10])
						self.im1  = ax1.imshow(depol, extent =(time_bnds.min(), time_bnds.max(), altitude.min(),
										  altitude.max()), aspect = 'auto', cmap = 'jet', interpolation = 'none',
										  origin = 'lowest', norm=LogNorm())

						ax1.set_xlabel('Time (Decimal hours)')
						ax1.set_ylabel('Height (km)')
						ax1.grid(True)
						majorLocator    = MultipleLocator(1)
						ax1.yaxis.set_major_locator(majorLocator)   
						cax1 = fig1.add_axes([0.89, 0.25, 0.03, 0.5])
						fig1.colorbar(self.im1, label = 'Atten Beta [$km^{-1}\,sr^{-1}$]', cax=cax1, format=l_f)
						plt.subplots_adjust(left=None, right=0.85, top=None, bottom=None)
	 
						self.im1.set_clim(1e-05, 1e-01)

						Rects_get = np.asarray(Rects_get)
						for i in range(len(Rects_get)):
							x0 = Rects_get[i,0]
							y0 = Rects_get[i,3]
							w  = Rects_get[i,1] - x0
							h  = Rects_get[i,2] - y0
							self.rect = Rectangle((x0, y0), w, h)
							ax1.add_patch(self.rect)
							self.rect.set_picker(5)
							self.rect.set_edgecolor('red')
							self.rect.set_linestyle('solid')
							self.rect.set_facecolor('none')
							self.rect.set_linewidth('5')

						self.canvas.draw()
						Fig1.append(fig1)
						Axes1.append(ax1)
						Im1.append(self.im1)
						ax1.format_coord = lambda x, y: ''

						aid1 = self.canvas1.mpl_connect('button_press_event', self.on_press)
						bid1 = self.canvas1.mpl_connect('button_release_event', self.on_release)
						cid1 = self.canvas1.mpl_connect('motion_notify_event', self.on_motion)
						fid1 = self.canvas1.mpl_connect('pick_event', self.patch_picker)

						toolbar1 = NavigationToolbar2TkAgg(self.canvas1, Toolbar1)
						toolbar1.update()
						toolbar1.pack(side=tk.BOTTOM)

						fig0, ax0 = plt.subplots() # Draw perpendicular 
						self.canvas0 = FigureCanvasTkAgg(fig0, master=Frame0)
						self.canvas0.get_tk_widget().grid(row=3, column=0, columnspan=5, rowspan=2, sticky=tk.NSEW)
						Frame0.columnconfigure(0, weight=1)
						Frame0.rowconfigure(3, weight=1)

						ax0.set_xlim([0, 24])
						ax0.set_ylim([0, 10])
						self.im0  = ax0.imshow(perp_beta, extent =(time_bnds.min(), time_bnds.max(), altitude.min(),
										  altitude.max()), aspect = 'auto', cmap = 'jet', interpolation = 'none',
										  origin = 'lowest', norm=LogNorm())

						ax0.set_xlabel('Time (Decimal hours)')
						ax0.set_ylabel('Height (km)')
						ax0.grid(True)
						majorLocator    = MultipleLocator(1)
						ax0.yaxis.set_major_locator(majorLocator)                     
						cax0 = fig0.add_axes([0.89, 0.25, 0.03, 0.5])
						fig0.colorbar(self.im0, label = 'Atten Beta [$km^{-1}\,sr^{-1}$]', cax=cax0, format=l_f)
						plt.subplots_adjust(left=None, right=0.85, top=None, bottom=None)

						self.im0.set_clim(1e-05, 1e-01)

						Rects_get = np.asarray(Rects_get)
						for i in range(len(Rects_get)):
							x0 = Rects_get[i,0]
							y0 = Rects_get[i,3]
							w  = Rects_get[i,1] - x0
							h  = Rects_get[i,2] - y0
							self.rect = Rectangle((x0, y0), w, h)
							ax0.add_patch(self.rect)
							self.rect.set_picker(5)
							self.rect.set_edgecolor('red')
							self.rect.set_linestyle('solid')
							self.rect.set_facecolor('none')
							self.rect.set_linewidth('5')

						self.canvas.draw()
						Fig0.append(fig0)
						Axes0.append(ax0)
						Im0.append(self.im0)
						ax0.format_coord = lambda x, y: ''

						aid0 = self.canvas0.mpl_connect('button_press_event', self.on_press)
						bid0 = self.canvas0.mpl_connect('button_release_event', self.on_release)
						cid0 = self.canvas0.mpl_connect('motion_notify_event', self.on_motion)
						fid0 = self.canvas0.mpl_connect('pick_event', self.patch_picker)

						toolbar0 = NavigationToolbar2TkAgg(self.canvas0, Toolbar0)
						toolbar0.update()
						toolbar0.pack(side=tk.BOTTOM)

						clim_min, clim_max = self.im.get_clim()  # Get colourbar limits

						entry_frame = tk.Frame(self) # Create frame to store colourbar editing buttons/entries
						entry_frame.grid(row=4, column=2, sticky=tk.NSEW)
						Entry_Frame.append(entry_frame) # Store location of frame

						label = tk.Label(entry_frame, textvariable = self.var9)
						label.grid(row=0, column=0, columnspan=3)

						button = tk.Button(entry_frame, text="-", command=self.max_colorbar_down)
						button.grid(row=1, column=0, sticky=tk.W)

						global entry_max_entry
						entry_max_entry = tk.Entry(entry_frame, width=6)
						entry_max_entry.grid(row=1, column=1)
						entry_max_entry.insert(0, '%.e' % clim_max)
						entry_max_entry.bind('<Return>', self.change_clim) # Bind <Return> to function
						entry_max_index.append(entry_max_entry)

						button = tk.Button(entry_frame, text="+", command=self.max_colorbar_up)
						button.grid(row=1, column=2, sticky=tk.W)

						label = tk.Label(entry_frame, textvariable = '')
						label.grid(row=2, column=0, columnspan=3)
						
						label = tk.Label(entry_frame, textvariable = self.var8)
						label.grid(row=3, column=0, columnspan=3)

						button = tk.Button(entry_frame, text="-", command=self.min_colorbar_down)
						button.grid(row=4, column=0, sticky=tk.W)
	   
						global entry_min_entry
						entry_min_entry = tk.Entry(entry_frame, width=6)
						entry_min_entry.grid(row=4, column=1)
						entry_min_entry.insert(0, '%.e' % clim_min)
						entry_min_entry.bind('<Return>', self.change_clim)
						entry_min_index.append(entry_min_entry)

						button = tk.Button(entry_frame, text="+", command=self.min_colorbar_up)
						button.grid(row=4, column=2, sticky=tk.W)
					
			except RuntimeError:
				pass
				tkMessageBox.showinfo("Error", "Invalid pathname in list file") # Option to have error window    
    '''
    Fuction to link shortcut to OpenSavedFile function
    '''
    def OpenSavedFileShortcut(self, event):
        return self.OpenSavedFile()
    
    '''
    Function to add/plot new file. For comments see OpenSavedFile
    '''
    def AddNewFile(self): 
        ftypes = (('Netcdf files', '*.nc'), ('All files', '*.*'))
        dlg = tkFileDialog.askopenfilenames(filetypes = ftypes, initialdir = initialdir_new)
        list_filenames =  self.tk.splitlist(dlg)
        oldlen_filenames = len(Filenames)
        if len(list_filenames)+len(Filenames) > 31:
            tkMessageBox.showinfo("Memory overload", "Please select %s files or less." %(31-oldlen_filenames))
        elif len(list_filenames) > 31:
            tkMessageBox.showinfo("Memory overload", "Please select 31 files or less.")
        else: 
            for i in range(len(list_filenames)):
                global filename
                filename = list_filenames[i]
                Filenames.append(filename)
                partial_path = os.path.split(filename)[1] # Splitting file name path
                small_path = partial_path[:-1]  # Only getting file name out of path
                SmallFilenames.append(small_path)

                if len(Filenames) >= 0:
                    FileNumberEntry.delete(0, tk.END)
                    FileNumberEntry.insert(0, '%s' % (len(Filenames)))
                    
                self.var10.set('File name: %s' % small_path)
                self.var5.set('Number of files: %s' % (len(Filenames)))
                self.update_idletasks()

                List.append([[], [], []])

                Toolbar = tk.Frame(self)
                Toolbar.grid(row=4, column=0,  sticky=tk.W)
                toolbar_frames.append(Toolbar)

                Toolbar1 = tk.Frame(self)
                Toolbar1.grid(row=4, column=0,  sticky=tk.W)
                toolbar_frames1.append(Toolbar1)

                Toolbar0 = tk.Frame(self)
                Toolbar0.grid(row=4, column=0,  sticky=tk.W)
                toolbar_frames0.append(Toolbar0)

                Frame=tk.Frame(self)
                Frame.grid(row=2, column = 0, columnspan = 2, rowspan=2, sticky=tk.NSEW)
                frames.append(Frame)

                Frame1=tk.Frame(self)  #(perpendicular)
                Frame1.grid(row=2, column = 0, columnspan = 2, rowspan=2, sticky=tk.NSEW)
                frames1.append(Frame1)

                Frame0=tk.Frame(self)
                Frame0.grid(row=2, column = 0, columnspan = 2, rowspan=2, sticky=tk.NSEW)
                frames0.append(Frame0)

                nc_file = Dataset(filename, 'r')
                global time_bnds
                global para_beta
                global perp_beta
                global altitude
                global depol
                time_bnds = nc_file.variables['time_bnds'][:]
                altitude = nc_file.variables['altitude'][:]
                para_beta = np.swapaxes((nc_file.variables['para_beta'][:]), 0, 1)
                perp_beta = np.swapaxes((nc_file.variables['perp_beta'][:]), 0, 1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    depol = perp_beta/para_beta
                para_beta[np.isnan(para_beta)]=0
                perp_beta[np.isnan(perp_beta)]=0
                depol[np.isnan(depol)]=0

                fig, ax = plt.subplots()
                self.canvas = FigureCanvasTkAgg(fig, master=Frame)
                self.canvas.get_tk_widget().grid(row=3, column=0, columnspan=5, rowspan=2, sticky=tk.NSEW)
                Frame.columnconfigure(0, weight=1)
                Frame.rowconfigure(3, weight=1)

                ax.set_xlim([0, 24])
                ax.set_ylim([0, 10])
                self.im  = ax.imshow(para_beta, extent =(time_bnds.min(), time_bnds.max(), altitude.min(),
                                altitude.max()), aspect = 'auto', cmap = 'jet', interpolation = 'none',
                                origin = 'lowest', norm=LogNorm())
                ax.grid(True)
                majorLocator    = MultipleLocator(1)
                ax.yaxis.set_major_locator(majorLocator)
                
                ax.set_xlabel('Time (Decimal hours)')
                ax.set_ylabel('Height (km)')
                cax = fig.add_axes([0.89, 0.25, 0.03, 0.5])
                fig.colorbar(self.im, label = 'Atten Beta [$km^{-1}\,sr^{-1}$]', cax=cax, format=l_f)#, ticks = lvls)
                plt.subplots_adjust(left=None, right=0.85, top=None, bottom=None)

                self.im.set_clim(1e-05, 1e-01)

                self.canvas.draw()
                Fig.append(fig)
                Axes.append(ax)
                Im.append(self.im)

                ax.format_coord = lambda x, y: ''

                toolbar = NavigationToolbar2TkAgg(self.canvas, Toolbar)
                toolbar.update()
                toolbar.pack(side=tk.TOP)

                aid = self.canvas.mpl_connect('button_press_event', self.on_press)
                bid = self.canvas.mpl_connect('button_release_event', self.on_release)
                cid = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
                fid = self.canvas.mpl_connect('pick_event', self.patch_picker)

                fig1, ax1 = plt.subplots() # Draw depol 
                self.canvas1 = FigureCanvasTkAgg(fig1, master=Frame1)
                self.canvas1.get_tk_widget().grid(row=3, column=0, columnspan=5, rowspan=2, sticky=tk.NSEW)
                Frame1.columnconfigure(0, weight=1)
                Frame1.rowconfigure(3, weight=1)

                ax1.set_xlim([0, 24])
                ax1.set_ylim([0, 10])
                self.im1  = ax1.imshow(depol, extent =(time_bnds.min(), time_bnds.max(), altitude.min(),
                                  altitude.max()), aspect = 'auto', cmap = 'jet', interpolation = 'none',
                                  origin = 'lowest', norm=LogNorm())
                ax1.grid(True)
                majorLocator    = MultipleLocator(1)
                ax1.yaxis.set_major_locator(majorLocator)

                ax1.set_xlabel('Time (Decimal hours)')
                ax1.set_ylabel('Height (km)')
                cax1 = fig1.add_axes([0.89, 0.25, 0.03, 0.5])
                fig1.colorbar(self.im1, label = 'Atten Beta [$km^{-1}\,sr^{-1}$]', cax=cax1, format=l_f)
                plt.subplots_adjust(left=None, right=0.85, top=None, bottom=None)

                self.im1.set_clim(1e-05, 1e-01)

                self.canvas.draw()
                Fig1.append(fig1)
                Axes1.append(ax1)
                Im1.append(self.im1)
                ax1.format_coord = lambda x, y: ''

                aid1 = self.canvas1.mpl_connect('button_press_event', self.on_press)
                bid1 = self.canvas1.mpl_connect('button_release_event', self.on_release)
                cid1 = self.canvas1.mpl_connect('motion_notify_event', self.on_motion)
                fid1 = self.canvas1.mpl_connect('pick_event', self.patch_picker)

                toolbar1 = NavigationToolbar2TkAgg(self.canvas1, Toolbar1)
                toolbar1.update()
                toolbar1.pack(side=tk.BOTTOM)

                fig0, ax0 = plt.subplots()
                self.canvas0 = FigureCanvasTkAgg(fig0, master=Frame0)
                self.canvas0.get_tk_widget().grid(row=3, column=0, columnspan=5, rowspan=2, sticky=tk.NSEW)
                Frame0.columnconfigure(0, weight=1)
                Frame0.rowconfigure(3, weight=1)

                ax0.set_xlim([0, 24])
                ax0.set_ylim([0, 10])
                self.im0  = ax0.imshow(perp_beta, extent =(time_bnds.min(), time_bnds.max(), altitude.min(),
                                  altitude.max()), aspect = 'auto', cmap = 'jet', interpolation = 'none',
                                  origin = 'lowest', norm=LogNorm())

                ax0.set_xlabel('Time (Decimal hours)')
                ax0.set_ylabel('Height (km)')
                ax0.grid(True)
                majorLocator    = MultipleLocator(1)
                ax0.yaxis.set_major_locator(majorLocator)
                cax0 = fig0.add_axes([0.89, 0.25, 0.03, 0.5])
                fig0.colorbar(self.im0, label = 'Atten Beta [$km^{-1}\,sr^{-1}$]', cax=cax0, format=l_f)#, ticks = lvls)
                plt.subplots_adjust(left=None, right=0.85, top=None, bottom=None)

                self.im0.set_clim(1e-05, 1e-01)

                self.canvas0.draw()
                Fig0.append(fig0)
                Axes0.append(ax0)
                Im0.append(self.im0)
                ax0.format_coord = lambda x, y: ''

                aid0 = self.canvas0.mpl_connect('button_press_event', self.on_press)
                bid0 = self.canvas0.mpl_connect('button_release_event', self.on_release)
                cid0 = self.canvas0.mpl_connect('motion_notify_event', self.on_motion)
                fid0 = self.canvas0.mpl_connect('pick_event', self.patch_picker)
    
                toolbar0 = NavigationToolbar2TkAgg(self.canvas0, Toolbar0)
                toolbar0.update()
                toolbar0.pack(side=tk.BOTTOM)

                clim_min, clim_max = self.im.get_clim()

                entry_frame = tk.Frame(self)
                entry_frame.grid(row=4, column=2, sticky=tk.NSEW)
                Entry_Frame.append(entry_frame)
            
                label = tk.Label(entry_frame, textvariable = self.var9)
                label.grid(row=0, column=0, columnspan=3)
            
                button = tk.Button(entry_frame, text="-", command=self.max_colorbar_down)
                button.grid(row=1, column=0, sticky=tk.W)
            
                global entry_max_entry
                entry_max_entry = tk.Entry(entry_frame, width=6)
                entry_max_entry.grid(row=1, column=1)
                entry_max_entry.insert(0, '%.e' % clim_max)
                entry_max_entry.bind('<Return>', self.change_clim)
                entry_max_index.append(entry_max_entry)
            
                button = tk.Button(entry_frame, text="+", command=self.max_colorbar_up)
                button.grid(row=1, column=2, sticky=tk.W)

                label = tk.Label(entry_frame, textvariable = '')
                label.grid(row=2, column=0, columnspan=3)
            
                label = tk.Label(entry_frame, textvariable = self.var8)
                label.grid(row=3, column=0, columnspan=3)
            
                button = tk.Button(entry_frame, text="-", command=self.min_colorbar_down)
                button.grid(row=4, column=0, sticky=tk.W)
                       
                global entry_min_entry
                entry_min_entry = tk.Entry(entry_frame, width=6)
                entry_min_entry.grid(row=4, column=1)
                entry_min_entry.insert(0, '%.e' % clim_min)
                entry_min_entry.bind('<Return>', self.change_clim)
                entry_min_index.append(entry_min_entry)

                button = tk.Button(entry_frame, text="+", command=self.min_colorbar_up)
                button.grid(row=4, column=2, sticky=tk.W)

    '''
    Fuction to link shortcut to OpenSavedFile function
    '''
    def AddNewFileShortcut(self, event):
        return self.AddNewFile()

    '''
    Fuction to redraw figure when limits of the colourbar have been changed
    (Minimum colorbar limit has been increased by button)
    '''        
    def min_colorbar_up(self):
        frame_index   = (int(FileNumberEntry.get()) - 1) # Get index of current figure
        if self.var6.get() ==4: # Para
           toolbar = toolbar_frames[frame_index]
           frame = frames[frame_index]
           ax = Axes[frame_index]
           fig = Fig[frame_index]
           im = Im[frame_index]

        if self.var6.get() == 3: # Perp
           toolbar = toolbar_frames0[frame_index]
           frame = frames0[frame_index]
           ax = Axes0[frame_index]
           fig = Fig0[frame_index]
           im = Im0[frame_index]

        if self.var6.get() == 5: # Depol
           toolbar = toolbar_frames1[frame_index]
           frame = frames1[frame_index]
           ax = Axes1[frame_index]
           fig = Fig1[frame_index]
           im = Im1[frame_index]
        
        clim_min, clim_max = im.get_clim() # Get current colourbar limits
        new_clim_min = clim_min*10 # Get new colourbar limits
        im.set_clim(new_clim_min, clim_max) # Draw figure with new colourbar limits
        fig.canvas.draw()
        entry_min = entry_min_index[frame_index] # Change relevant entry in colourbar limits frame
        entry_min.delete(0, tk.END)
        entry_min.insert(0, '%.e' % new_clim_min)

    '''
    Fuction to redraw figure when limits of the colourbar have been changed
    (Minimum colorbar limit has been decreased by button)
    For comments see min_colorbar_up
    '''  
    def min_colorbar_down(self):
        frame_index   = (int(FileNumberEntry.get()) - 1)
        if self.var6.get() ==4: # Para
           toolbar = toolbar_frames[frame_index]
           frame = frames[frame_index]
           ax = Axes[frame_index]
           fig = Fig[frame_index]
           im = Im[frame_index]

        if self.var6.get() == 3: # Perp
           toolbar = toolbar_frames0[frame_index]
           frame = frames0[frame_index]
           ax = Axes0[frame_index]
           fig = Fig0[frame_index]
           im = Im0[frame_index]

        if self.var6.get() == 5: # Depol
           toolbar = toolbar_frames1[frame_index]
           frame = frames1[frame_index]
           ax = Axes1[frame_index]
           fig = Fig1[frame_index]
           im = Im1[frame_index]

        clim_min, clim_max = im.get_clim()
        new_clim_min = clim_min/10
        im.set_clim(new_clim_min, clim_max)
        fig.canvas.draw()
        entry_min = entry_min_index[frame_index]
        entry_min.delete(0, tk.END)
        entry_min.insert(0, '%.e' % new_clim_min)

    '''
    Fuction to redraw figure when limits of the colourbar have been changed
    (Maximum colorbar limit has been increased by button)
    For comments see min_colorbar_up
    '''  
    def max_colorbar_up(self):
        frame_index   = (int(FileNumberEntry.get()) - 1)
        if self.var6.get() ==4: # Para
           toolbar = toolbar_frames[frame_index]
           frame = frames[frame_index]
           ax = Axes[frame_index]
           fig = Fig[frame_index]
           im = Im[frame_index]

        if self.var6.get() == 3: # Perp
           toolbar = toolbar_frames0[frame_index]
           frame = frames0[frame_index]
           ax = Axes0[frame_index]
           fig = Fig0[frame_index]
           im = Im0[frame_index]

        if self.var6.get() == 5: # Depol
           toolbar = toolbar_frames1[frame_index]
           frame = frames1[frame_index]
           ax = Axes1[frame_index]
           fig = Fig1[frame_index]
           im = Im1[frame_index]
        
        clim_min, clim_max = im.get_clim()
        new_clim_max = clim_max*10
        im.set_clim(clim_min, new_clim_max)
        fig.canvas.draw()
        entry_max = entry_max_index[frame_index]
        entry_max.delete(0, tk.END)
        entry_max.insert(0, '%.e' % new_clim_max)

    '''
    Fuction to redraw figure when limits of the colourbar have been changed
    (Maximum colorbar limit has been decreased by button)
    For comments see min_colorbar_up
    '''  
    def max_colorbar_down(self):
        frame_index   = (int(FileNumberEntry.get()) - 1)
        if self.var6.get() ==4: # Para
           toolbar = toolbar_frames[frame_index]
           frame = frames[frame_index]
           ax = Axes[frame_index]
           fig = Fig[frame_index]
           im = Im[frame_index]

        if self.var6.get() == 3: # Perp
           toolbar = toolbar_frames0[frame_index]
           frame = frames0[frame_index]
           ax = Axes0[frame_index]
           fig = Fig0[frame_index]
           im = Im0[frame_index]

        if self.var6.get() == 5: # Depol
           toolbar = toolbar_frames1[frame_index]
           frame = frames1[frame_index]
           ax = Axes1[frame_index]
           fig = Fig1[frame_index]
           im = Im1[frame_index]
        
        clim_min, clim_max = im.get_clim()
        new_clim_max = clim_max/10
        im.set_clim(clim_min, new_clim_max)
        fig.canvas.draw()
        entry_max = entry_max_index[frame_index]
        entry_max.delete(0, tk.END)
        entry_max.insert(0, '%.e' % new_clim_max)

    '''
    Fuction to redraw figure when limits of the colourbar have been changed
    (Maximum or mimumum colorbar limit has altered by input into entry box followed by <Return>)
    '''  
    def change_clim(self, event):
        frame_index   = (int(FileNumberEntry.get()) - 1)
        if self.var6.get() ==4: # Para
           toolbar = toolbar_frames[frame_index]
           frame = frames[frame_index]
           ax = Axes[frame_index]
           fig = Fig[frame_index]
           im = Im[frame_index]

        if self.var6.get() == 3: # Perp
           toolbar = toolbar_frames0[frame_index]
           frame = frames0[frame_index]
           ax = Axes0[frame_index]
           fig = Fig0[frame_index]
           im = Im0[frame_index]

        if self.var6.get() == 5: # Depol
           toolbar = toolbar_frames1[frame_index]
           frame = frames1[frame_index]
           ax = Axes1[frame_index]
           fig = Fig1[frame_index]
           im = Im1[frame_index]
        
        entry_max = entry_max_index[frame_index]
        entry_min = entry_min_index[frame_index]
        
        clim_min = float(entry_min.get()) # Retrieve values in entry boxes
        clim_max = float(entry_max.get())
        im.set_clim(clim_min, clim_max) # Set colourbar limits to these values
        fig.canvas.draw()  # Redraw canvas

    '''
    Fuction to initialise a rectange patch when 'select areas' radiobutton is selected and canvas is clicked on
    '''              
    def on_press(self, event):
        self.is_pressed = True  # If canvas is clicked on
        if event.xdata is not None and event.ydata is not None:  # If canvas is clicked on
            fig_index = int(FileNumberEntry.get()) # Refer to relevant frame
            if self.var1.get() == 1:  # If 'Select Areas' is on
                if self.var6.get() ==4:   # Parallel
                    ax = Axes[fig_index-1] # Refer to relevant frame
                    fig = Fig[fig_index-1] # Refer to relevant frame
                if self.var6.get() == 3: # Perpendicular
                    ax = Axes0[fig_index-1]
                    fig = Fig0[fig_index-1]
                if self.var6.get() == 5: # Perpendicular
                    ax = Axes1[fig_index-1]
                    fig = Fig1[fig_index-1]
                self.rect = Rectangle((0,0), 0, 0) # Set rectangle to be invisible
                ax.add_patch(self.rect)  # Add patch(rectangle)
                self.rect.set_picker(5) # Turn on picking
                self.x0, self.y0 = event.xdata, event.ydata  # Set inital corner of rectangle to be positioned where you clicked
                self.rect.set_width(0) 
                self.rect.set_height(0)
                self.rect.set_xy((self.x0, self.y0))
                self.rect.set_edgecolor('red')
                self.rect.set_linestyle('solid')
                self.rect.set_facecolor('none')
                self.rect.set_linewidth('5')
                fig.canvas.draw()  # Draw canvas

    '''
    Fuction to add a moving (updating) rectange patch when 'select areas' radiobutton is selected and there is motion of the canvas. 
    For comments see on_press
    '''  
    def on_motion(self, event):
        try:
            if self.is_pressed:
                if event.xdata is not None and event.ydata is not None:
                    fig_index = int(FileNumberEntry.get())
                    if self.var1.get() == 1:
                        if self.var6.get() ==4:
                            ax = Axes[fig_index-1]
                            fig = Fig[fig_index-1]
                        if self.var6.get() == 3:
                            ax = Axes0[fig_index-1]
                            fig = Fig0[fig_index-1]
                        if self.var6.get() == 5: # Perpendicular
                            ax = Axes1[fig_index-1]
                            fig = Fig1[fig_index-1]
                        self.x1, self.y1 = event.xdata, event.ydata # Set other corner of rectangle to be where the mouse is now moved to
                        self.rect.set_width(self.x1 - self.x0) 
                        self.rect.set_height(self.y1 - self.y0)
                        self.rect.set_xy((self.x0, self.y0))
                        fig.canvas.draw()
        except AttributeError:
            pass

    '''
    Fuction to draw rectange patch mouse button is released. Coordinates of rectangle is added to List. 
    For comments see on_press 
    '''  
    def on_release(self, event):
        self.is_pressed = False # If button is released
        fig_index = int(FileNumberEntry.get())
        if self.var1.get() == 1:
            if self.var6.get() ==4:
                ax = Axes[fig_index-1]
                fig = Fig[fig_index-1]
            if self.var6.get() == 3:
                ax = Axes0[fig_index-1]
                fig = Fig0[fig_index-1]
            if self.var6.get() == 5: # Perpendicular
                ax = Axes1[fig_index-1]
                fig = Fig1[fig_index-1]
            self.x1, self.y1 = event.xdata, event.ydata
            self.rect.set_xy((self.x0, self.y0))
            coords =([self.x0, self.x1, self.y1, self.y0])
            fig.canvas.draw()
            if len(List[fig_index-1][1]) <= 0: # For comments see OpenSavedFile
                List[fig_index-1][0].append(coords) 
            else:
                List[fig_index-1][1][0].append(coords)  

    '''
    Fuction to be able to delete rectangles on click
    '''
    def patch_picker(self, event):
        self.is_pressed = True # If clicked
        if self.var1.get() == 2:  # If 'delete areas' radiobutton is selected
            fig_index = int(FileNumberEntry.get())
            if self.var6.get() == 4: 
                fig = Fig[fig_index-1]
            if self.var6.get() == 3:
                fig = Fig0[fig_index-1]
            if self.var6.get() == 5:
                fig = Fig1[fig_index-1]
            patch = event.artist            # Select patch
            height = patch.get_height()     # Get height of patch
            width = patch.get_width()
            x = patch.get_x()
            y = patch.get_y()
            box = ([x, x+width, y+height, y])   # Define box as patch
            List[fig_index-1][2].append(box) # Add coordinates of box to List 
            patch.remove()      # Remove box from canvas
            Remaining = [x for x in List[fig_index-1][0] if x not in (List[fig_index-1][2])]
            if len(List[fig_index-1][1]) > 0:  # If the Remaining sublist of List is not empty 
                List[fig_index-1][1] = []   # Make it empty 
            List[fig_index-1][1].append(Remaining) # Store location of remaining  rectangles 
            fig.canvas.show() # Redraw canvas

    '''
    Fuction to switch between parallel/perpendicular/depol figures
    '''
    def view(self):
        try:
            frame_index = int(FileNumberEntry.get())  # Get relevant file
            Remain = List[frame_index-1][1]  # Get coordinates of remaining boxes (remaining following deletion)
            Drawn = List[frame_index-1][0] # Get coordinates of boxes (only up to date if no boxes have been deleted)

            if len(Remain) <= 0:    
                Rerect = np.asarray(Drawn)
            else:
                Rerect = np.asarray(Remain)

            if Rerect.ndim >2:
                Rerect_flat = np.asarray([item for sublist in Rerect for item in sublist]) # Flatten sublist and convert to array if needed
            else:
                Rerect_flat = Rerect

            if self.var6.get() ==4:
               toolbar = toolbar_frames[frame_index-1]
               frame = frames[frame_index-1]
               ax = Axes[frame_index-1]
               fig = Fig[frame_index-1]
               im = Im[frame_index-1]
               clim_max = 1e-01
               clim_min = 1e-05
          
            if self.var6.get() == 3: # Perp
               toolbar = toolbar_frames0[frame_index-1]
               frame = frames0[frame_index-1]
               ax = Axes0[frame_index-1]
               fig = Fig0[frame_index-1]
               im = Im0[frame_index-1]
               clim_max = 1e-01
               clim_min = 1e-05

            if self.var6.get() == 5: # Depol
               toolbar = toolbar_frames1[frame_index-1]
               frame = frames1[frame_index-1]
               ax = Axes1[frame_index-1]
               fig = Fig1[frame_index-1]
               im = Im1[frame_index-1]
               clim_max = 1e+01
               clim_min = 1e-01
               
            toolbar.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
            frame.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'

            im.set_clim(clim_min, clim_max)
            entry_min = entry_min_index[frame_index-1] # Change relevant entry in colourbar limits frame
            entry_min.delete(0, tk.END)
            entry_min.insert(0, '%.e' % clim_min)
            entry_max = entry_max_index[frame_index-1] # Change relevant entry in colourbar limits frame
            entry_max.delete(0, tk.END)
            entry_max.insert(0, '%.e' % clim_max)
            fig.canvas.draw()
            
            for i in range(len(Rerect_flat)): # Redraw rectangles on other figure, based on stored coordinates
                x0 = Rerect_flat[i,0]
                y0 = Rerect_flat[i,3]
                w  = Rerect_flat[i,1] - x0
                h  = Rerect_flat[i,2] - y0
                self.rect = Rectangle((x0, y0), w, h)
                ax.add_patch(self.rect)
                self.rect.set_picker(5)
                self.rect.set_edgecolor('red')
                self.rect.set_linestyle('solid')
                self.rect.set_facecolor('none')
                self.rect.set_linewidth('5')
                fig.canvas.draw()
        except ValueError:
            pass

    '''
    Fuction to switch between files based on Plus button
    '''        
    def Plus(self):
        try:
            value = int(FileNumberEntry.get()) # Get index of current file
            if (value+1) > len(Filenames): # Set limit 
                new_value = value
            else:
                new_value = value+1

                if self.var6.get() ==4:
                    ax = Axes[new_value-1]
                    fig = Fig[new_value-1]
                    frame = frames[new_value-1]
                    toolbar = toolbar_frames[new_value-1]
                    im = Im[new_value-1]
                    clim_max = 1e-01
                    clim_min = 1e-05                    
                if self.var6.get() == 3:
                    ax = Axes0[new_value-1]
                    fig = Fig0[new_value-1]
                    frame = frames0[new_value-1]
                    toolbar = toolbar_frames0[new_value-1]
                    im = Im0[new_value-1]
                    clim_max = 1e-01
                    clim_min = 1e-05                     
                if self.var6.get() == 5:
                    ax = Axes1[new_value-1]
                    fig = Fig1[new_value-1]
                    frame = frames1[new_value-1]
                    toolbar = toolbar_frames1[new_value-1]
                    im = Im1[new_value-1]
                    clim_max = 1e+01
                    clim_min = 1e-01
                                         
                clim_entry = Entry_Frame[new_value-1]
                clim_entry.tkraise()                
                im.set_clim(clim_min, clim_max)
                entry_min = entry_min_index[new_value-1] # Change relevant entry in colourbar limits frame
                entry_min.delete(0, tk.END)
                entry_min.insert(0, '%.e' % clim_min)
                entry_max = entry_max_index[new_value-1] # Change relevant entry in colourbar limits frame
                entry_max.delete(0, tk.END)
                entry_max.insert(0, '%.e' % clim_max)
                fig.canvas.draw()
                toolbar.tkraise()
                toolbar.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
                frame.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
                #self.var5.set('Number of files: %d' % new_value)

                Remain = List[new_value-1][1]  # Get coordinates of remaining boxes (remaining following deletion)
                Drawn = List[new_value-1][0] # Get coordinates of boxes (only up to date if no boxes have been deleted)
                if len(Remain) <= 0:    
                    Rerect = np.asarray(Drawn)
                else:
                    Rerect = np.asarray(Remain)
    
                if Rerect.ndim >2:
                    Rerect_flat = np.asarray([item for sublist in Rerect for item in sublist]) # Flatten sublist and convert to array if needed
                else:
                    Rerect_flat = Rerect
                for i in range(len(Rerect_flat)): # Redraw rectangles on other figure, based on stored coordinates
                    x0 = Rerect_flat[i,0]
                    y0 = Rerect_flat[i,3]
                    w  = Rerect_flat[i,1] - x0
                    h  = Rerect_flat[i,2] - y0
                    self.rect = Rectangle((x0, y0), w, h)
                    ax.add_patch(self.rect)
                    self.rect.set_picker(5)
                    self.rect.set_edgecolor('red')
                    self.rect.set_linestyle('solid')
                    self.rect.set_facecolor('none')
                    self.rect.set_linewidth('5')
                    fig.canvas.draw()
                filename= SmallFilenames[new_value-1]
                self.var10.set('File Name: %s' % filename)
                self.update_idletasks()
            FileNumberEntry.delete(0, tk.END)
            FileNumberEntry.insert(0, new_value)
        except ValueError:
            pass

    '''
    Fuction to switch between files based on changing value in Entry box
    '''    
    def ChangeFrame(self, event):
        value = int(FileNumberEntry.get())
        if value>(len(Filenames)):
            tkMessageBox.showinfo("Incorrect case number", "There are only %s cases" %(len(Filenames)))
        if len(Filenames) >= 1: 
            new_value = value
            if self.var6.get() ==4:
                ax = Axes[new_value-1]
                fig = Fig[new_value-1]
                frame = frames[new_value-1]
                toolbar = toolbar_frames[new_value-1]
                im = Im[new_value-1]
                clim_max = 1e-01
                clim_min = 1e-05
            if self.var6.get() == 3:
                ax = Axes0[new_value-1]
                fig = Fig0[new_value-1]
                frame = frames0[new_value-1]
                toolbar = toolbar_frames0[new_value-1]
                im = Im0[new_value-1]
                clim_max = 1e-01
                clim_min = 1e-05
            if self.var6.get() == 5:
                ax = Axes1[new_value-1]
                fig = Fig1[new_value-1]
                frame = frames1[new_value-1]
                toolbar = toolbar_frames1[new_value-1]
                im = Im1[new_value-1]
                clim_max = 1e+01
                clim_min = 1e-01
                
            clim_entry = Entry_Frame[new_value-1]
            clim_entry.tkraise()
            im.set_clim(clim_min, clim_max)
            entry_min = entry_min_index[new_value-1] # Change relevant entry in colourbar limits frame
            entry_min.delete(0, tk.END)
            entry_min.insert(0, '%.e' % clim_min)
            entry_max = entry_max_index[new_value-1] # Change relevant entry in colourbar limits frame
            entry_max.delete(0, tk.END)
            entry_max.insert(0, '%.e' % clim_max)
            fig.canvas.draw()
            toolbar.tkraise()
            toolbar.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
            frame.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
            #self.var5.set('Number of files: %d' % new_value)

            Remain = List[new_value-1][1]  # Get coordinates of remaining boxes (remaining following deletion)
            Drawn = List[new_value-1][0] # Get coordinates of boxes (only up to date if no boxes have been deleted)
            if len(Remain) <= 0:    
                Rerect = np.asarray(Drawn)
            else:
                Rerect = np.asarray(Remain)

            if Rerect.ndim >2:
                Rerect_flat = np.asarray([item for sublist in Rerect for item in sublist]) # Flatten sublist and convert to array if needed
            else:
                Rerect_flat = Rerect
            for i in range(len(Rerect_flat)): # Redraw rectangles on other figure, based on stored coordinates
                x0 = Rerect_flat[i,0]
                y0 = Rerect_flat[i,3]
                w  = Rerect_flat[i,1] - x0
                h  = Rerect_flat[i,2] - y0
                self.rect = Rectangle((x0, y0), w, h)
                ax.add_patch(self.rect)
                self.rect.set_picker(5)
                self.rect.set_edgecolor('red')
                self.rect.set_linestyle('solid')
                self.rect.set_facecolor('none')
                self.rect.set_linewidth('5')
                fig.canvas.draw()
            filename= SmallFilenames[new_value-1]
            self.var10.set('File Name: %s' % filename)
            self.update_idletasks()
            
    '''
    Fuction to switch between files based on Minus button
    '''            
    def Minus(self):
        try:
            value = int(FileNumberEntry.get())
            if len(Filenames) == 0:
                new_value = 0
            elif value-1 < 1:
                new_value = 1
            else:
                new_value = value-1
                if self.var6.get() ==4:
                    ax = Axes[new_value-1]
                    fig = Fig[new_value-1]
                    frame = frames[new_value-1]
                    toolbar = toolbar_frames[new_value-1]
                    im = Im[new_value-1]
                    clim_max = 1e-01
                    clim_min = 1e-05
                                        
                if self.var6.get() == 3:
                    ax = Axes0[new_value-1]
                    fig = Fig0[new_value-1]
                    frame = frames0[new_value-1]
                    toolbar = toolbar_frames0[new_value-1]
                    im = Im0[new_value-1]
                    clim_max = 1e-01
                    clim_min = 1e-05
                                        
                if self.var6.get() == 5:
                    ax = Axes1[new_value-1]
                    fig = Fig1[new_value-1]
                    frame = frames1[new_value-1]
                    toolbar = toolbar_frames1[new_value-1]
                    im = Im1[new_value-1]
                    clim_max = 1e+01
                    clim_min = 1e-01
                    
                clim_entry = Entry_Frame[new_value-1]
                clim_entry.tkraise()
                im.set_clim(clim_min, clim_max)
                entry_min = entry_min_index[new_value-1] # Change relevant entry in colourbar limits frame
                entry_min.delete(0, tk.END)
                entry_min.insert(0, '%.e' % clim_min)
                entry_max = entry_max_index[new_value-1] # Change relevant entry in colourbar limits frame
                entry_max.delete(0, tk.END)
                entry_max.insert(0, '%.e' % clim_max)                
                toolbar.tkraise()
                toolbar.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
                frame.tkraise() # Raise appropriate frame depending on selction of Radiobutton 'Parallel' or 'Perpendicular'
                #self.var5.set('Number of files: %d' % new_value)

                Remain = List[new_value-1][1]  # Get coordinates of remaining boxes (remaining following deletion)
                Drawn = List[new_value-1][0] # Get coordinates of boxes (only up to date if no boxes have been deleted)
                if len(Remain) <= 0:    
                    Rerect = np.asarray(Drawn)
                else:
                    Rerect = np.asarray(Remain)
    
                if Rerect.ndim >2:
                    Rerect_flat = np.asarray([item for sublist in Rerect for item in sublist]) # Flatten sublist and convert to array if needed
                else:
                    Rerect_flat = Rerect
                for i in range(len(Rerect_flat)): # Redraw rectangles on other figure, based on stored coordinates
                    x0 = Rerect_flat[i,0]
                    y0 = Rerect_flat[i,3]
                    w  = Rerect_flat[i,1] - x0
                    h  = Rerect_flat[i,2] - y0
                    self.rect = Rectangle((x0, y0), w, h)
                    ax.add_patch(self.rect)
                    self.rect.set_picker(5)
                    self.rect.set_edgecolor('red')
                    self.rect.set_linestyle('solid')
                    self.rect.set_facecolor('none')
                    self.rect.set_linewidth('5')
                    fig.canvas.draw()
                filename= SmallFilenames[new_value-1]
                self.var10.set('File Name: %s' % filename)
                self.update_idletasks()
            FileNumberEntry.delete(0, tk.END)
            FileNumberEntry.insert(0, new_value)
        except ValueError:
            pass

    '''
    Fuction to save rectangle coordinates as .lis file
    '''            
    def save(self):
        file = tkFileDialog.asksaveasfilename(initialdir = initialdir_save) # Open tk window for saving file
        try:
            if file is '':
                pass
            else:
                outfile = open('%s' % file, 'a')  # Open file saved as, opened in append mode 
                outfile.write('   %s\n' % len(Filenames)) # Write number of files saved in list file 
                for i in range(len(Filenames)):
                    fig_index = i+1
                    if file is None:
                        return
                    if len(List[fig_index-1][1]) <= 0:  # If no rectangles have been deleted
                        Array = List[fig_index-1][0]
                        outfile.write('%s\n' % Filenames[fig_index-1]) # Write file name
                        outfile.write('   %s\n' % len(Array)) # Write number of rectangles
                        for i in range(len(Array)):
                            Array1 = [float(Decimal('%.4f' % e)) for e in Array[i] ] # Convert string to float 
                            list_flat  = ",  ".join( repr(e) for e in Array1) # Join floats with comma
                            outfile.write('  %s\n' % list_flat) # Write coordinates to .lis file
                    else: # If rectangles have been deleted
                        Array = ((List[fig_index-1][1]))
                        array_flat = [item for sublist in Array for item in sublist]
                        outfile.write('%s\n' % Filenames[fig_index-1])
                        outfile.write('   %s\n' % len(array_flat))
                        for i in range(len(array_flat)):
                            Array1 = [float(Decimal('%.4f' % e)) for e in array_flat[i]]
                            list_flat  = ",  ".join( repr(e) for e in Array1)
                            outfile.write('  %s\n' % list_flat)
                outfile.close() # Close file
        except TypeError:
            pass
            

    '''
    Fuction for shortcut to save function
    '''    
    def saveShortcut(self, event):  
        return self.save()

    '''
    Fuction to delete current frame
    '''    
    def delete_frame(self):
        try:
            value = int(FileNumberEntry.get()) # Get index of current frame
            index = value-1
            frame = frames[index]
            frame0 = frames0[index]
            frame1 = frames1[index]

            toolbar = toolbar_frames[index]
            toolbar0 = toolbar_frames0[index]
            toolbar1 = toolbar_frames1[index]
            entryframe = Entry_Frame[index]
            fig = Fig[index]
            fig0 = Fig0[index]
            fig1 = Fig1[index]
            
            frame.destroy() # Destroy frames
            frame0.destroy()
            frame1.destroy()
            toolbar.destroy()
            toolbar0.destroy()
            toolbar1.destroy()
            entryframe.destroy()
            plt.close(fig) # Remove plots from memory
            plt.close(fig0)
            plt.close(fig1)
            
            del Filenames[index] # Delete stored locations
            del frames[index]
            del frames0[index]
            del frames1[index]
            del Fig[index]
            del Fig0[index]
            del Fig1[index]
            del Axes[index]
            del Axes0[index]
            del Axes1[index]
            del toolbar_frames[index]
            del toolbar_frames0[index]
            del toolbar_frames1[index]
            del List[index]
            del Im[index]
            del Im0[index]
            del Im1[index]
            del Entry_Frame[index]
            del entry_max_index[index] 
            del entry_min_index[index] 
            
            FileNumberEntry.delete(0, tk.END) # update entry
            if index < 1:
                FileNumberEntry.insert(0, 1)            
            else:
                FileNumberEntry.insert(0, '%d' % index)              
            self.var5.set('Number of files: %d' % len(Filenames)) # Update label
            self.var10.set('File name: %s' % SmallFilenames[index-1])
            self.update_idletasks()
            del SmallFilenames[index]

            new_value = int(FileNumberEntry.get())            
            if self.var6.get() ==4:
                ax = Axes[new_value-1]
                fig = Fig[new_value-1]
                frame = frames[new_value-1]
                toolbar = toolbar_frames[new_value-1]
            if self.var6.get() == 3:
                ax = Axes0[new_value-1]
                fig = Fig0[new_value-1]
                frame = frames0[new_value-1]
                toolbar = toolbar_frames0[new_value-1]
            if self.var6.get() == 5:
                ax = Axes1[new_value-1]
                fig = Fig1[new_value-1]
                frame = frames1[new_value-1]
                toolbar = toolbar_frames1[new_value-1]


            clim_entry = Entry_Frame[new_value-1]
            toolbar.tkraise()
            clim_entry.tkraise()
            toolbar.tkraise() 
            frame.tkraise() 
            if len(List)>1:
                Remain = List[new_value-1][1]  # Get coordinates of remaining boxes (remaining following deletion)
                Drawn = List[new_value-1][0] # Get coordinates of boxes (only up to date if no boxes have been deleted)
                if len(Remain) <= 0:    
                    Rerect = np.asarray(Drawn)
                else:
                    Rerect = np.asarray(Remain)
        
                if Rerect.ndim >2:
                    Rerect_flat = np.asarray([item for sublist in Rerect for item in sublist]) # Flatten sublist and convert to array if needed
                else:
                    Rerect_flat = Rerect
                
                Rerect_flat_list = [item for sublist in Rerect_flat for item in sublist]
                if not Rerect_flat_list:
                    pass
                else:
                    for i in range(len(Rerect_flat)): # Redraw rectangles on other figure, based on stored coordinates
                        x0 = Rerect_flat[i,0]
                        y0 = Rerect_flat[i,3]
                        w  = Rerect_flat[i,1] - x0
                        h  = Rerect_flat[i,2] - y0
                        self.rect = Rectangle((x0, y0), w, h)
                        ax.add_patch(self.rect)
                        self.rect.set_picker(5)
                        self.rect.set_edgecolor('red')
                        self.rect.set_linestyle('solid')
                        self.rect.set_facecolor('none')
                        self.rect.set_linewidth('5')
                        fig.canvas.draw()            
                
        except ValueError:
            pass
            #tkMessageBox.showinfo("Error", "No files to remove") # Option to have error window

    '''
    Fuction to show popup window with names of files
    '''    
    def PrintCases(self):
        toplevel = tk.Toplevel() # Create toplevel window
        Files_print = '\n'.join(Filenames) # Get filenames
        label1 = tk.Label(toplevel, text='%s' %Files_print, height = 0, width=100) # Add filenames to label
        label1.pack()

'''
Fuction to quit GUI
'''    
def quit():
    root.destroy()

'''
Main function to create GUI
'''    
def main():
    root = tk.Tk()
    root.resizable(True, True)
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight = 1)
    root.rowconfigure(2, weight = 1)
    app = Example(root)
    root.mainloop()

if __name__ == '__main__':
    main()

'''
01
02
05
06
08
13
14
'''
