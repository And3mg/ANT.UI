#pyinstaller --windowed --icon="Icon.ico" --paths="." ANT.UI.py
import numpy as np
import matplotlib.pyplot as plt
from os import listdir, mkdir
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from customtkinter import CTkButton as Button
from customtkinter import CTkCanvas as Canvas
from customtkinter import CTkFrame as Frame
from customtkinter import CTkComboBox as Combobox
from customtkinter import CTkLabel as Label
from customtkinter import CTkSlider as Slider
from customtkinter import CTkCheckBox as Checkbutton
from customtkinter import CTkToplevel as Toplevel
from customtkinter import CTkOptionMenu as Menu
from customtkinter import CTkTextbox as Textbox
from customtkinter import CTkEntry as Text


import customtkinter
customtkinter.set_appearance_mode("System")  # Modes: system (default), light, dark
customtkinter.set_default_color_theme("blue")  # Themes: blue (default), dark-blue, green

__version__ = '1.0.0'

class Scale:

    def __init__(self, frame, name, from_n, to_n, resol):
       
        self.label = Label(frame,  text = f'{name} = 0', font = ('Arial', 12))
        self.slider = Slider(frame, from_= from_n , to= to_n, number_of_steps=int((to_n-from_n)/resol),
                            width = 150, orientation='horizontal', command=self.update)
        self.slider_plus = Button(frame,  text = '+', width=5, command = self.Move_slider_plus)
        self.slider_minus = Button(frame,  text = '-', width=5, command = self.Move_slider_minus)  

        self.label.pack(side = 'top', expand = True)
        self.slider_minus.pack(side = 'left')
        self.slider.pack(side = 'left', expand = True)
        self.slider_plus.pack(side = 'left')

        self.resol = resol
        self.name = name
        if resol == 1:
            self.decimals = 0
        elif resol == .1:
            self.decimals = 1

    def Move_slider_plus(self):
        self.slider.set(self.slider.get()+self.resol)
        self.update(self.slider.get()+self.resol)

    def Move_slider_minus(self):
        self.slider.set(self.slider.get()-self.resol)
        self.update(self.slider.get()-self.resol)

    def get(self):
        return float(self.slider.get())
    def update(self, val):
        self.label.configure(text= f'{self.name} = {np.round(val, self.decimals)}')
    def set(self,x):
        self.slider.set(x)
        self.update(x)




#%% GLOBAl MEMORY
root = None
fig = None
ax = None
#%% DICTIONARY INDEX

from _DICS import *

#%% CODIGO PLOTEAR


def sphere(axes, poss, r, c):
    u, v = np.mgrid[0:2*np.pi:10j, 0:np.pi:6j]
    x = r*np.cos(u)*np.sin(v) +poss[0]
    y = r*np.sin(u)*np.sin(v) +poss[1]
    z = r*np.cos(v) +poss[2]
    #axes.plot_wireframe(x, y, z, color= c, linewidth = 2)
    axes.plot_surface(x, y, z, color= c, linewidth = 2)
    return

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.25*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    ax.set_xlabel(r'X $[\AA]$')
    ax.set_ylabel(r'Y $[\AA]$')
    ax.set_zlabel(r'Z $[\AA]$')

    ax.set_facecolor(bg_color)

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False  
    ax.xaxis.pane.set_edgecolor(bg_color)
    ax.yaxis.pane.set_edgecolor(bg_color)
    ax.zaxis.pane.set_edgecolor(bg_color)
    #ax.axis('off')
    ax.grid(False)
    return

#%% ATOM MANAGMENT

def List_from_dir(Dir, extensions = ''):
    archives = []
    for x in listdir(f'./{Dir}'):
        temp = x.split('.')
        if len(temp) > 1:
            if extensions:
                if temp[-1] == extensions:
                    archives.append(temp[0])
            else:
                archives.append(x)
    return archives

def find_center(atoms):
    if len(atoms) >= 0:
        position = []
        for atom in atoms:
            position.append(np.array(atom[-3:]))
        return np.mean(np.array(position), axis = 0)
    else: return None

class xyz_object():
    
    def __init__(self, route):
        self.name = route.split('/')[-1].split('.')[0]
        self.mol_flag = False
        self.U = {}
        if route.split('/')[0] == 'Molecules':
            self.mol_flag = True
        with open(route) as text: #Leemos la molecula
            self.N = int(text.readline())
            temp = text.readline()
            if self.mol_flag and temp != '\n' and not ('PureMetalic' in self.name):
                [self.U.update(y) for y in [{x.split()[0]:x.split()[1]} for x in temp.split('|')]]
            self.poss = np.zeros((self.N, 3))
            self.elements = np.zeros((self.N), dtype = '<U2')
            line = text.readline()
            i = 0
            while line != '':
                self.elements[i] = line.split()[0]
                self.poss[i] = [float(x) for x in line.split()[1:] ]
                line = text.readline()
                i+=1
        self.CM = find_center(self.poss)
        self.poss = self.poss-self.CM
        self.CM = find_center(self.poss)
        
        if self.mol_flag:
            self.Molecular_atoms = np.unique(self.elements)
            
        self.semi_infinite_layer = np.zeros(self.N, dtype = 'bool')    
        self.extremes = np.zeros(self.N, dtype = 'bool')                                                     
                                                 
        if not self.mol_flag:
            self.poss = self.poss-self.CM
            self.CM = find_center(self.poss)
            if self.name.split('_')[0] == 'Block':
                self.semi_infinite_layer = np.logical_xor(self.poss[:,2].min() == self.poss[:,2], True)
                self.extremes = np.logical_xor(self.poss[:,2].min() == self.poss[:,2], True)
            else:
                self.semi_infinite_layer = (self.poss[:,2].max() == self.poss[:,2]) + \
                    (np.unique(self.poss[:,2])[np.unique(self.poss[:,2]).argsort()][-2] == self.poss[:,2])
                self.extremes = (self.poss[:,2].max() == self.poss[:,2]) 


        self.spatial_coords = np.array([0,0,0])
        self.euler_angles = np.array([0,0,0])
        self.angles_mem = np.array([0,0,0])
        self.find_neighbors()
        
    def chiral_mol(self): 
        temp = np.zeros((self.N, 4))
        for i,x in enumerate(self.poss):
            temp[i] = np.array([x[0],x[1],x[2],1])
        specular_matrix = np.array([[1,0,0,0],
                                   [0,-1,0,0],
                                   [0,0,1,0],
                                   [0,0,0,1]])
        temp = temp@specular_matrix
        for i,x in enumerate(temp):
            self.poss[i] = np.array([x[0],x[1],x[2]])

    def reflection_z(self):
        refZ = np.array([[1, 0, 0],[0, 1, 0], [0, 0, -1]])
        self.poss = self.poss@refZ
       
    def true_coordinates(self):

        phi = (self.euler_angles[0]- self.angles_mem[0])*2*np.pi/360
        theta = (self.euler_angles[1]- self.angles_mem[1])*2*np.pi/360
        chi = (self.euler_angles[2]- self.angles_mem[2])*2*np.pi/360
        
        rotationXi = np.array([[1,0,0],
                             [0, np.cos(phi), np.sin(phi)],
                             [0, -np.sin(phi), np.cos(phi)]])
        
        rotationYi = np.array([[np.cos(theta), 0 , np.sin(theta)],
                             [0,1,0],
                             [-np.sin(theta), 0, np.cos(theta)]])
    
        rotationZi = np.array([[np.cos(chi), np.sin(chi), 0],
                             [-np.sin(chi), np.cos(chi), 0],
                             [0, 0, 1]])
        
        for i,x in enumerate(self.poss):
            self.poss[i] = x@rotationXi@rotationYi@rotationZi
        
        self.angles_mem = self.euler_angles.copy()
        
        return self.poss+self.spatial_coords
  
    def change_elem(self, nw_elem):
        rate = radius[nw_elem]/radius[self.elements[0]]
        
        for i in range(self.N):
            self.elements[i] = nw_elem
            
        self.poss = self.poss*rate
        
    def find_neighbors(self):
        self.neigh_count = np.zeros(self.N, dtype = 'int')
        for i,r0 in enumerate(self.poss):
            temp = np.linalg.norm(self.poss - r0, axis = 1)
            if not self.extremes[i]:
                self.neigh_count[i] = np.sum(temp<2*radius[f'{self.elements[i]}'])-1
            else:
                self.neigh_count[i] = 10        
        
    def plot_object(self, ax):
        
        temp = self.true_coordinates()
        '''
        for i in range(self.N):
            sphere(ax, temp[i],
                    radius[self.elements[i]], colores[self.elements[i]])
        '''
        if not self.mol_flag:
            for i in range(self.N):
                if self.neigh_count[i]<10:
                    sphere(ax, temp[i],
                           radius[self.elements[i]], colores[self.elements[i]])
        elif self.mol_flag:
            for i in range(self.N):
                sphere(ax, temp[i],
                       radius[self.elements[i]], colores[self.elements[i]])
        


#%% Output generation


from _DFT import *

SOC_FLAG = True
try:
    from _SOC import *
except:
    SOC_FLAG = False

from _OUT import *


#%% TKINTER WINDOW 

class object_controls():
    global root, ax
    
    def __init__(self, xyz_object, frame):
        self.obj = xyz_object
                
        #Frame 1 for controling the electrode in cartesian axes
        self.frame_1 = Frame(frame)
        self.frame_1.configure(width=720,height=50)
        self.frame_1.pack(pady = 5, fill = 'both')
        
        self.frame_x = Frame(self.frame_1, width = 240)
        self.slider_x = Scale(self.frame_x, 'x', -20, 20, 0.1)
        self.frame_x.pack(side = 'left', expand = True)

        self.frame_y = Frame(self.frame_1, width = 240)
        self.slider_y = Scale(self.frame_y, 'y', -20, 20, 0.1)
        self.frame_y.pack(side = 'left', expand = True)

        self.frame_z = Frame(self.frame_1, width = 240)
        self.slider_z = Scale(self.frame_z, 'z', -20, UserInfo['ZMAX'], 0.1)
        self.frame_z.pack(side = 'left', expand = True)
        
        #Frame 2 for controling the electrode in euler angles
        self.frame_2 = Frame(frame)
        self.frame_2.configure(width=720,height=50)
        self.frame_2.pack(pady = 5, fill = 'both')
        
        self.frame_phi = Frame(self.frame_2, width = 240)
        self.slider_phi = Scale(self.frame_phi, 'phi', -180, 180, 1)
        self.frame_phi.pack(side = 'left', expand = True)

        self.frame_theta = Frame(self.frame_2, width = 240)
        self.slider_theta = Scale(self.frame_theta, 'theta', -180, 180, 1)
        self.frame_theta.pack(side = 'left', expand = True)

        self.frame_psi = Frame(self.frame_2, width = 240)
        self.slider_psi = Scale(self.frame_psi, 'psi', -180, 180, 1)
        self.frame_psi.pack(side = 'left', expand = True)
        
        if self.obj.mol_flag:
            self.frame_3 = Frame(frame)
            self.frame_3.configure(width=720,height=25)
            self.frame_3.pack(pady = 5, fill = 'both')
            
            self.chiral_btn = Button(self.frame_3,  text = 'Specular Reflection', command=self.chiral)
            self.chiral_btn.pack(expand = True)
            self.chiral_flag = False
            
    def chiral(self):
        if not self.chiral_flag:
            self.chiral_btn.configure(text_color = 'black')
            self.chiral_flag = not self.chiral_flag
        else:
            self.chiral_btn.configure(text_color = 'white')
            self.chiral_flag = not self.chiral_flag
        self.obj.chiral_mol()     

    
    def destroy(self):
        self.frame_1.destroy()
        self.frame_2.destroy()
        if self.obj.mol_flag: self.frame_3.destroy()
    
    def update(self):
        self.obj.spatial_coords = np.array([self.slider_x.get(),self.slider_y.get(),self.slider_z.get()])
        self.obj.euler_angles = np.array([self.slider_phi.get(),self.slider_theta.get(),self.slider_psi.get()])
        self.obj.plot_object(ax)


class toplevel_opt():

    def __init__(self, opt_data, molecule, returner):
        self.window = Toplevel()
        self.window.geometry("500x360")
        self.window.title("Optimization options")
        
        #Frames for electrode data
        self.frame_1 = Frame(self.window)
        self.frame_1.configure(width=500,height=50)
        self.frame_1.pack(pady = 0, fill = 'both')

        self.frame_2 = Frame(self.window)
        self.frame_2.configure(width=500,height=50)
        self.frame_2.pack(pady = 10, fill = 'both')

        self.frame_3 = Frame(self.window)
        self.frame_3.configure(width=500,height=50)
        self.frame_3.pack(pady = 0, fill = 'both')

        self.frame_4 = Frame(self.window)
        self.frame_4.configure(width=500,height=50)
        self.frame_4.pack(pady = 10, fill = 'both')


        #Top electrode
        self.label_top = Label(self.frame_1,  text = 'Top electrode: Bethe Layer', font = ('Arial', 18))
        self.menu_bethe_top = Menu(self.frame_2, values=['fix','z'])
        self.label_top1 = Label(self.frame_3,  text = 'Top electrode', font = ('Arial', 18))
        self.menu_top = Menu(self.frame_4, values=['fix','all','z'])

        self.label_top.pack(side= 'left', expand=True)
        self.menu_bethe_top.pack(side = 'left', expand  =True)
        self.label_top1.pack(side='left', expand = True)
        self.menu_top.pack(side='left', expand = True)
        
        #Bottom electrode


        self.label_bottom = Label(self.frame_1,  text = 'Bottom tip: Bethe Layer', font = ('Arial', 18))
        self.menu_bethe_bottom = Menu(self.frame_2, values=['fix','z'])
        self.label_bottom1 = Label(self.frame_3,  text = 'Bottom electrode', font = ('Arial', 18))
        self.menu_bottom = Menu(self.frame_4, values=['fix','all','z'])

        self.label_bottom.pack(side= 'right', expand=True)
        self.menu_bethe_bottom.pack(side = 'right', expand  =True)
        self.label_bottom1.pack(side= 'right', expand=True)
        self.menu_bottom.pack(side='right', expand = True)

        #Frame for molecule
        self.frame_mol = Frame(self.window)
        self.frame_mol.configure(width=500,height=100)
        self.frame_mol.pack(pady = 25, fill = 'both')

        self.label_mol = Label(self.frame_mol,  text = 'Molecule: ', font = ('Arial', 18))
        self.menu_mol = Menu(self.frame_mol, values=['fix','all','z','all, some only z'], command= self.mol_choose_if)

        self.label_mol.pack(side = 'left', expand = True)
        self.menu_mol.pack(side = 'right', expand = True)

        #Frame for chains
        self.frame_chains = Frame(self.window)
        self.frame_chains.configure(width=500,height=100)
        self.frame_chains.pack(pady = 10, fill = 'both')

        self.label_chains = Label(self.frame_chains,  text = 'Chains: ', font = ('Arial', 18))
        self.menu_chains = Menu(self.frame_chains, 
                                values=['none', 'first step', 'every step'])

        self.label_chains.pack(side = 'left', expand = True)
        self.menu_chains.pack(side = 'right', expand = True)
        
        #Frame for Nopt
        self.frame_Nopt = Frame(self.window)
        self.frame_Nopt.configure(width=500,height=100)
        self.frame_Nopt.pack(pady = 10, fill = 'both')

        self.label_Nopt = Label(self.frame_Nopt,  text = 'Number of optimization steps: ', font = ('Arial', 14))
        self.text_Nopt = Text(self.frame_Nopt, width = 60, height = 2, font=('consolas', 14))

        self.label_Nopt.pack(side = 'left', expand=True)
        self.text_Nopt.pack(side = 'left', expand=True)


        self.opt_in_ANT = BooleanVar()
        self.check_opt_in_ANT = Checkbutton(self.frame_Nopt, text = 'Opt. in ANT', variable = self.opt_in_ANT, 
                                   onvalue = True, offvalue = False)
        self.check_opt_in_ANT.pack(side = 'left', expand=True)

        #Frame for Ok button
        self.frame_btn = Frame(self.window)
        self.frame_btn.configure(width=500,height=50)
        self.frame_btn.pack(side = 'bottom', fill = 'both')

        self.ok_btn = Button(self.frame_btn,  text = 'Ok', command = self.ok )
        self.ok_btn.pack(side = 'bottom', expand = True)

        #Setting the current configuration
        self.menu_bethe_top.set(opt_data['Top_Bethe'])
        self.menu_bethe_bottom.set(opt_data['Bottom_Bethe'])
        self.menu_top.set(opt_data['Top'])
        self.menu_bottom.set(opt_data['Bottom'])
        self.menu_mol.set(opt_data['Molecule'])
        self.atom_list = opt_data['Molecule_z_atoms']
        self.menu_chains.set(opt_data['Chains'])
        self.text_Nopt.insert(END, str(opt_data['Opt_steps']))


        self.returner = returner
        self.molecule = molecule

    def mol_choose_if(self, option):
        if option == 'all, some only z':
            subtoplevel_opt(self.molecule, self.aux_mol_choose)
            pass

    def aux_mol_choose(self, zlist):
        self.atom_list = zlist
        
    
    def ok(self):

        self.returner({'Top_Bethe': self.menu_bethe_top.get(), 'Top': self.menu_top.get(), 
                       'Bottom_Bethe': self.menu_bethe_bottom.get(), 'Bottom': self.menu_bottom.get(),
                         'Molecule': self.menu_mol.get(), 'Molecule_z_atoms': self.atom_list, 
                         'Chains': self.menu_chains.get(), 'Opt_steps': int(self.text_Nopt.get()),
                         'Opt_in_ANT': self.check_opt_in_ANT.get()})
        self.window.destroy()

        
class subtoplevel_opt():
    
    def __init__(self, molecule, returner):

        self.window = Toplevel()
        self.window.geometry("500x360")
        self.window.title("Choose z atoms")

        self.data_returner = returner

        #Display and write frames
        self.frame_text = Frame(self.window)
        self.frame_text.configure(width=500,height=300)
        self.frame_text.pack(pady = 0, fill = 'both')

        self.textbox_mol = Textbox(self.frame_text, width=480, height= 260, activate_scrollbars=True, 
                                   state = 'normal', font=('Arial', 14))

        #Write the text
        self.textbox_mol.insert('0.0', 'List of atoms in molecule, write index to select: \n \n')
        self.textbox_mol.insert('end', 'Index   Elem   x   y   z  \n')
        i = 0
        for elem, poss in zip(molecule.elements, molecule.poss):
            self.textbox_mol.insert('end', f'{i}   {elem}   {np.round(poss[0],3)}   {np.round(poss[1],3)}   {np.round(poss[2],3)}  \n')
            i+=1
        self.textbox_mol.state = 'disable'
        self.textbox_mol.pack(side='top', expand = True, fill = 'both')

        self.entry_id = Text(self.frame_text, width = 200, height = 2, font=('consolas', 14))
        self.entry_id.pack(side='bottom', expand = True, fill = 'both')
        
        #OK button
        self.frame_btn = Frame(self.window)
        self.frame_btn.configure(width=500,height=50)
        self.frame_btn.pack(side = 'bottom', fill = 'both')

        self.ok_btn = Button(self.frame_btn,  text = 'Ok', command = self.ok )
        self.ok_btn.pack(side = 'bottom', expand = True)


     
    def ok(self):
        self.data_returner([int(x) for x in self.entry_id.get().split()])
        self.window.destroy()
        return
        
class toplevel_grid():
    
    def __init__(self, xmin, xmax, ymin, ymax, z, bottom_electrode, molecule, data_returner):
        self.window = Toplevel()
        self.window.geometry("500x500")
        self.window.title("Grid DFT input asistant")
        
        #Frame for xdata
        self.frame_x = Frame(self.window)
        self.frame_x.configure(width=500,height=20)
        self.frame_x.pack(side = 'top', fill = 'x')
        
        h = 2
        self.label_xfrom = Label(self.frame_x,  text = 'From xmin =', font = ('Arial', 14))
        self.text_xmin = Text(self.frame_x, width = 60, height = h, font=('consolas', 14))
        self.label_xto = Label(self.frame_x,  text = ' to xmax =', font = ('Arial', 14))
        self.text_xmax = Text(self.frame_x, width = 60, height = h, font=('consolas', 14))
        self.label_nx = Label(self.frame_x,  text = ' nx =', font = ('Arial', 14))
        self.text_nx = Text(self.frame_x, width = 24, height = h, font=('consolas', 14))

        self.label_xfrom.pack(side = 'left', expand=True)
        self.text_xmin.pack(side = 'left', expand=True)
        self.label_xto.pack(side = 'left', expand=True)
        self.text_xmax.pack(side = 'left', expand=True)
        self.label_nx.pack(side = 'left', expand=True)
        self.text_nx.pack(side = 'left', expand=True)

        self.text_xmax.insert(END, str(round(xmax, 2)))
        self.text_xmin.insert(END, str(round(xmin, 2)))
        self.text_nx.insert(END, str(2))

        #Frame for ydata
        self.frame_y = Frame(self.window)
        self.frame_y.configure(width=500,height=50)
        self.frame_y.pack(side = 'top', fill = 'x')
        
        self.label_yfrom = Label(self.frame_y,  text = 'From ymin =', font = ('Arial', 14))
        self.text_ymin = Text(self.frame_y, width = 60, height = h, font=('consolas', 14))
        self.label_yto = Label(self.frame_y,  text = ' to ymax =', font = ('Arial', 14))
        self.text_ymax = Text(self.frame_y, width = 60, height = h, font=('consolas', 14))
        self.label_ny = Label(self.frame_y,  text = ' ny =', font = ('Arial', 14))
        self.text_ny = Text(self.frame_y, width = 24, height = h, font=('consolas', 14))

        self.label_yfrom.pack(side = 'left', expand=True)
        self.text_ymin.pack(side = 'left', expand=True)
        self.label_yto.pack(side = 'left', expand=True)
        self.text_ymax.pack(side = 'left', expand=True)
        self.label_ny.pack(side = 'left', expand=True)
        self.text_ny.pack(side = 'left', expand=True)

        self.text_ymax.insert(END, str(round(ymax, 2)))
        self.text_ymin.insert(END, str(round(ymin, 2)))
        self.text_ny.insert(END, str(2))

        #Frame for ok btn
        self.frame_btn = Frame(self.window)
        self.frame_btn.configure(width=500,height=50)
        self.frame_btn.pack(side = 'top', fill = 'both')

        self.ok_btn = Button(self.frame_btn,  text = 'Ok', command = self.ok )
        self.ok_btn.pack(side = 'left', expand = True)


        #Frame for canvas
        self.frame_canvas = Frame(self.window)
        self.frame_canvas.configure(width=500,height=350)
        self.frame_canvas.pack(side = 'top', fill = 'both')
        
        Fig = plt.Figure(figsize=(5, 5), dpi=100)
        self.ax = Fig.add_subplot(111)
        Fig.set_facecolor('#000000')
        self.ax.set_facecolor('#000000')
        self.ax.tick_params(axis='x', colors='white')    
        self.ax.tick_params(axis='y', colors='white')  
        self.ax.spines['left'].set_color('white')
        self.ax.spines['top'].set_color('white') 
        self.ax.spines['right'].set_color('white')
        self.ax.spines['bottom'].set_color('white') 

        self.bottom_electrode = bottom_electrode
        self.molecule = molecule
        
        
        self.canvas = FigureCanvasTkAgg(Fig, master = self.frame_canvas)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        
        
        self.after_process = None
        
        self.xmin, self.xmax, self.nx = 0,0,2
        self.ymin, self.ymax, self.ny = 0,0,2
        self.z = z

        self.data_returner = data_returner
        
    def ok(self):
        self.data_returner([[self.xmin, self.xmax, self.nx], [self.ymin, self.ymax, self.ny], self.z])
        self.window.after_cancel(self.after_process)
        self.window.destroy()
        return
        
    def update(self):
        self.ax.clear()
        self.ax.plot(self.bottom_electrode[:,0],self.bottom_electrode[:,1], 'or', markersize = 10)
        self.ax.plot(self.molecule[:,0],self.molecule[:,1], 'ob', markersize = 5)
        
        try:
            self.xmin, self.xmax, self.nx = float(self.text_xmin.get()), float(self.text_xmax.get()), int(self.text_nx.get())
            self.ymin, self.ymax, self.ny = float(self.text_ymin.get()), float(self.text_ymax.get()), int(self.text_ny.get())
        except: pass
        
        x = np.linspace(self.xmin, self.xmax, self.nx)
        y = np.linspace(self.ymin, self.ymax, self.ny)
        X,Y = np.meshgrid(x,y)

        self.ax.plot(X,Y,'g')
        self.ax.plot(X.T, Y.T,'g')
        self.canvas.draw()
        self.canvas.blit()
        self.canvas.flush_events()        
        
        self.after_process = self.window.after(500, self.update)

class toplevel_pull():
    
    def __init__(self, data_returner):
        self.window = Toplevel()
        self.window.geometry("500x200")
        
        self.window.title("Pull electrodes DFT asistant")
        
        #Frame for zdata
        self.frame_z = Frame(self.window)
        self.frame_z.configure(width=500,height=20)
        self.frame_z.pack(pady = 50, fill = 'x')
        
        h = 2
        self.label_dz = Label(self.frame_z,  text = 'dz =', font = ('Arial', 14))
        self.text_dz = Text(self.frame_z, width = 60, height = h, font=('consolas', 14))
        self.label_nz = Label(self.frame_z,  text = '  nz =', font = ('Arial', 14))
        self.text_nz = Text(self.frame_z, width = 24, height = h, font=('consolas', 14))


        self.label_dz.pack(side = 'left', expand=True)
        self.text_dz.pack(side = 'left', expand=True)
        self.label_nz.pack(side = 'left', expand=True)
        self.text_nz.pack(side = 'left', expand=True)

        self.text_dz.insert(END, str(0.5))
        self.text_nz.insert(END, str(5))

        #Frame for checkbutton
        self.frame_chk = Frame(self.window)
        self.frame_chk.configure(width=500,height=50)
        self.frame_chk.pack(side = 'top', fill = 'x')
        
        self.label_both = Label(self.frame_chk,  text = 'Both electrodes', font = ('Arial', 14))
        self.flag_both = BooleanVar()
        self.check_both = Checkbutton(self.frame_chk, text='', variable = self.flag_both, onvalue = True, offvalue = False)

        self.label_both.pack(side = 'left', expand=True)
        self.check_both.pack(side = 'left', expand=True)

        #Frame for ok btn
        self.frame_btn = Frame(self.window)
        self.frame_btn.configure(width=500,height=50)
        self.frame_btn.pack(side = 'top', fill = 'both')

        self.ok_btn = Button(self.frame_btn,  text = 'Ok', command = self.ok )
        self.ok_btn.pack(side = 'left', expand = True)
        
        self.data_returner = data_returner
        
    def ok(self):
        self.data_returner(float(self.text_dz.get()), int(self.text_nz.get()), self.flag_both.get())
        self.window.destroy()
        return
        
class toplevel_rotate():
    
    def __init__(self, data_returner):
        self.window = Toplevel()
        self.window.geometry("500x500")
        self.window.title("Rotatation DFT asistant")
        
        #Frame for zdata
        self.frame_rot = Frame(self.window)
        self.frame_rot.configure(width=500,height=20)
        self.frame_rot.pack(pady = 50, fill = 'x')
        
        h = 2
        self.label_from = Label(self.frame_rot,  text = 'From ', font = ('Arial', 14))
        self.text_from = Text(self.frame_rot, width = 60, height = h, font=('consolas', 14))
        self.label_to = Label(self.frame_rot,  text = '  to ', font = ('Arial', 14))
        self.text_to = Text(self.frame_rot, width = 24, height = h, font=('consolas', 14))
        self.label_by = Label(self.frame_rot,  text = '  by', font = ('Arial', 14))
        self.text_by = Text(self.frame_rot, width = 24, height = h, font=('consolas', 14))


        self.label_from.pack(side = 'left', expand=True)
        self.text_from.pack(side = 'left', expand=True)
        self.label_to.pack(side = 'left', expand=True)
        self.text_to.pack(side = 'left', expand=True)
        self.label_by.pack(side = 'left', expand=True)
        self.text_by.pack(side = 'left', expand=True)

        self.text_from.insert(END, str(0))
        self.text_to.insert(END, str(90))
        self.text_by.insert(END, str(15))

        #Frame for checkbutton
        self.frame_chk_top = Frame(self.window)
        self.frame_chk_top.configure(width=500,height=50)
        self.frame_chk_top.pack(side = 'top', fill = 'x')
        
        self.flag_top = BooleanVar()
        self.check_top = Checkbutton(self.frame_chk_top, text = 'Top electrode', variable = self.flag_top, onvalue = True, offvalue = False)
        self.check_top.pack(side = 'left', expand=True)

        self.frame_chk_mol = Frame(self.window)
        self.frame_chk_mol.configure(width=500,height=50)
        self.frame_chk_mol.pack(side = 'top', fill = 'x')
        
        self.flag_mol = BooleanVar()
        self.check_mol = Checkbutton(self.frame_chk_mol, text = 'Molecule', variable = self.flag_mol, onvalue = True, offvalue = False)
        self.check_mol.pack(side = 'left', expand=True)

        self.frame_chk_bottom = Frame(self.window)
        self.frame_chk_bottom.configure(width=500,height=50)
        self.frame_chk_bottom.pack(side = 'top', fill = 'x')
        
        self.flag_bottom = BooleanVar()
        self.check_bottom = Checkbutton(self.frame_chk_bottom, text = 'Bottom electrode', variable = self.flag_bottom, onvalue = True, offvalue = False)
        self.check_bottom.pack(side = 'left', expand=True)

        self.frame_chk_SOC = Frame(self.window)
        self.frame_chk_SOC.configure(width=500,height=50)
        self.frame_chk_SOC.pack(side = 'top', fill = 'x')
        
        self.flag_SOC = BooleanVar()
        self.check_SOC = Checkbutton(self.frame_chk_SOC, text = 'SOC', variable = self.flag_SOC, onvalue = True, offvalue = False)
        self.check_SOC.pack(side = 'left', expand=True)

        #Frame for ok btn
        self.frame_btn = Frame(self.window)
        self.frame_btn.configure(width=500,height=50)
        self.frame_btn.pack(side = 'top', fill = 'both')

        self.ok_btn = Button(self.frame_btn,  text = 'Ok', command = self.ok )
        self.ok_btn.pack(side = 'left', expand = True)
        
        self.data_returner = data_returner
        
    def ok(self):
        self.data_returner(float(self.text_from.get()), float(self.text_to.get()), float(self.text_by.get()),
                           self.flag_top.get(), self.flag_mol.get(), self.flag_bottom.get(), self.flag_SOC.get())
        self.window.destroy()
        return


class tkinter_window():
    global root, ax, fig
     
    def __init__(self):
        self.after_process = None
        self.grid = None
        
        #Frame for choosing electrodes
        self.frame_label_combo_elec = Frame(root)
        self.frame_label_combo_elec.configure(width=720,height=100)
        self.frame_label_combo_elec.pack(pady = 0, fill = 'both')
        
        self.frame_combo_elec = Frame(root)
        self.frame_combo_elec.configure(width=720,height=100)
        self.frame_combo_elec.pack(pady = 10, fill = 'both')
        
        
        #Top
        self.label_combo_top = Label(self.frame_label_combo_elec,  text = 'Top Electrode', font = ('Arial', 16))
        self.combo_top_elem = Combobox(self.frame_combo_elec, state="readonly", values = Electrode_elem, justify='center', text_color='black', width=70)
        self.combo_top_type = Combobox(self.frame_combo_elec, state="readonly", values = List_from_dir('ElectrodesDFT','xyz'), justify='center', text_color='black')
        
        self.label_combo_top.pack(side = 'left', expand = True)
        self.combo_top_type.pack(side = 'left', expand = True)
        self.combo_top_elem.pack(side = 'left', expand = True)
        
        self.combo_top_elem.set('Au')
        self.combo_top_type.set('Point_FCC_[111]')
        
        self.up_electrode = xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')

        
        #Bottom
        self.label_combo_bottom = Label(self.frame_label_combo_elec,  text = 'Bottom Electrode', font = ('Arial', 16))
        self.combo_bottom_elem = Combobox(self.frame_combo_elec, state="readonly", values = Electrode_elem, justify='center', text_color='black', width=70)
        self.combo_bottom_type = Combobox(self.frame_combo_elec, state="readonly", values = List_from_dir('ElectrodesDFT','xyz'), justify='center', text_color='black')
        
        self.label_combo_bottom.pack(side = 'right', expand = True)
        self.combo_bottom_type.pack(side = 'left', expand = True)
        self.combo_bottom_elem.pack(side = 'left', expand = True)
        
        self.combo_bottom_elem.set('Au')
        self.combo_bottom_type.set('Point_FCC_[111]')
        
        self.down_electrode = xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')
        self.down_electrode.euler_angles = np.array([0,180, 60])
        self.down_electrode.spatial_coords = np.array([0,0,-11])


        #Frame 0 for title of Electrode controls
        self.frame_0 = Frame(root)
        self.frame_0.configure(width=720,height=100)
        self.frame_0.pack(pady = 30, fill = 'both')

        self.label_electrode = Label(self.frame_0,  text = 'Top electrode controls', font = ('Arial', 16))
        self.label_electrode.pack(fill = 'x')
        
        #Frames 1 and 25
        self.top_controls = object_controls(self.up_electrode, self.frame_0)
        self.top_controls.slider_z.set(20) #For a cleaner initial view
        
        
        #Frame for choosing molecule
        self.frame_label_combo = Frame(root)
        self.frame_label_combo.configure(width=720,height=100)
        self.frame_label_combo.pack(side = 'top', fill = 'both')
        
        self.frame_combo = Frame(root)
        self.frame_combo.configure(width=720,height=100)
        self.frame_combo.pack(side = 'top', fill = 'both')

        self.label_combo = Label(self.frame_label_combo,  text = 'Molecule', font = ('Arial', 16))
        self.combo = Combobox(self.frame_combo, state="readonly", values=List_from_dir('Molecules','xyz'), justify='center' , text_color='black', width=400)
        self.label_combo.pack(side = 'left', expand = True)
        self.combo.pack(side = 'left', expand = True)
        init_mol = List_from_dir('Molecules','xyz')[0]
        self.combo.set(init_mol)

        self.molecule = xyz_object(f'Molecules/{init_mol}.xyz')

        #Optimization options
        self.opt_data = {'Top_Bethe': 'fix', 'Top': 'fix', 'Bottom_Bethe': 'fix', 'Bottom': 'fix',
                         'Molecule': 'fix', 'Molecule_z_atoms': [], 'Chains': 'none', 'Opt_steps': 10,
                         'Opt_in_ANT': False}

        self.btn_opt = Button(self.frame_combo,  text = 'Optimization Options', command = self.opt_toplevel)
        self.btn_opt.pack(side = 'right', expand = True)


        #self.flag_opt = BooleanVar()
        #self.check_opt = Checkbutton(self.frame_combo, variable = self.flag_opt, onvalue = True, offvalue = False, text='Optimize molecule')
        #self.check_opt.pack(side = 'right', expand=True)

        #Frame 3 for title of Molecule controls       
        self.frame_3 = Frame(root)
        self.frame_3.configure(width=720,height=100)
        self.frame_3.pack(side = 'top', fill = 'both')

        self.label_electrode = Label(self.frame_3,  text = 'Molecule controls', font = ('Arial', 16))
        self.label_electrode.pack(fill = 'x')
        
        #Frames 4 and 5
        self.mol_controls = object_controls(self.molecule, self.frame_3)
        
        #Frame for DFT options
        self.frame_DFT_opts = Frame(root)
        self.frame_DFT_opts.configure(width=720,height=100)
        self.frame_DFT_opts.pack(side = 'top', fill = 'both')
        
        self.DFT_opts = {'Functional': 'HSE06', 'Shell': 'Close Shell', 'DFTU': False, 'EABS': False}
        
        self.combo_funct = Combobox(self.frame_DFT_opts, state="readonly", command=self.update_DFT,
                                    values = [x for x in Functionals.keys()], justify='center', text_color='black', width=100)
        self.combo_shell = Combobox(self.frame_DFT_opts, state="readonly", command=self.update_DFT,
                                    values = [x for x in Shell.keys()], justify='center', text_color='black', width=140)

        self.combo_shell.pack(side = 'left', expand = True)
        self.combo_shell.set('Close Shell')
        self.combo_funct.pack(side = 'left', expand = True)
        self.combo_funct.set('HSE06')
        
        self.DFT_U = BooleanVar()
        self.check_U = Checkbutton(self.frame_DFT_opts, text = 'DFT+U', variable = self.DFT_U, 
                                   onvalue = True, offvalue = False, command=self.upt_dft_chk,)
        self.check_U.pack(side = 'left', expand=True)

        self.EABS = BooleanVar()
        self.check_EABS = Checkbutton(self.frame_DFT_opts, text = 'Abs. Energy', variable = self.EABS, 
                                   onvalue = True, offvalue = False, command=self.upt_dft_chk,)
        self.check_EABS.pack(side = 'left', expand=True)

        #Frame with the buttons for assintant tools
        self.frame_assis = Frame(root)
        self.frame_assis.configure(width=720,height=200)
        self.frame_assis.pack(pady = 25, fill = 'both')

        self.btn_grid = Button(self.frame_assis,  text = 'Grid Assistant', command = self.grid_toplevel)
        self.btn_rotate = Button(self.frame_assis,  text = 'Rotation Assistant', command = self.rotate_toplevel)
        self.btn_pull = Button(self.frame_assis,  text = 'Pull Assistant', command = self.pull_toplevel)

        self.btn_grid.pack(side = 'left', expand = True)
        self.btn_rotate.pack(side = 'left', expand = True)
        self.btn_pull.pack(side = 'left', expand = True)

        #Frame with the buttons for output creation
        self.frame_btn = Frame(root)
        self.frame_btn.configure(width=720,height=200)
        self.frame_btn.pack(side = 'bottom', fill = 'both')
        
        self.xyz_btn = Button(self.frame_btn,  text = 'Create xyz', command = self.create_xyz)
        self.xyz_btn.pack(side = 'left', expand = True)
        
        self.ant_btn = Button(self.frame_btn,  text = 'Create ANT', command = self.create_ant)
        self.ant_btn.pack(side = 'left', expand = True)
        
        #self.reax_btn = Button(self.frame_btn,  text = 'REAX Input', command = self.create_reax)
        #self.reax_btn.pack(side = 'right', expand = True)
        
        if SOC_FLAG:
            self.soc_btn = Button(self.frame_btn,  text = 'SOC ANT', command = self.create_soc)
            self.soc_btn.pack(side = 'left', expand = True)

        
        
    def create_xyz(self):
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)
        
        Create_xyz_output([self.down_electrode, self.molecule, self.up_electrode], 'Outputs/'+name, name)

        return
    
    def create_ant(self):
        #Create Dir tree
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)
        
        Create_output([self.down_electrode, self.molecule, self.up_electrode], 
                      'Outputs/'+name, name, self.opt_data, self.DFT_opts)
        return
        
    def create_reax(self):
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)
        
        Create_reax_output([self.down_electrode, self.molecule, self.up_electrode], 'Outputs/'+name, name)
        return

    def create_soc(self):
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)
        
        Create_output_soc([self.down_electrode, self.molecule, self.up_electrode], 'Outputs/'+name,
                          name, self.opt_data, self.DFT_opts)
        return
    
    def update_opt_options(self, opt_data):
        self.opt_data = opt_data

    def opt_toplevel(self):
        toplevel_opt(self.opt_data, self.molecule, self.update_opt_options)

    def grid_toplevel(self):
        temp = find_center(self.up_electrode.true_coordinates())
        
        tpg = toplevel_grid(-abs(temp[0]), abs(temp[0]), -abs(temp[1]), abs(temp[1]), temp[2],
                      self.down_electrode.true_coordinates(), self.molecule.true_coordinates(), self.create_grid_dft)
        tpg.update()

    
    def create_grid_dft(self, grid):
        #Create Dir tree
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)
        
        Create_grid_output([self.down_electrode, self.molecule, self.up_electrode], 'Outputs/'+name, name, 
                           grid, self.opt_data, self.DFT_opts)
    
        
    def pull_toplevel(self):
        toplevel_pull(self.create_pull_dft)

    
    def create_pull_dft(self, dz, nz, both):
        #Create Dir tree
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)
        
        Create_pull_output([self.down_electrode, self.molecule, self.up_electrode], 'Outputs/'+name, 
                           name, dz, nz, both, self.opt_data, self.DFT_opts)
    
    def rotate_toplevel(self):
        toplevel_rotate(self.create_rot_dft)

    def create_rot_dft(self, frm, to, by, top, mol, bottom, SOC):
        angles = np.arange(frm, to+by, by)
        boolean_data = (bottom, mol, top)

        #Create Dir tree
        exisisting_outputs = listdir('Outputs')
        i = 1
        for x in exisisting_outputs:
            if x.split('-')[0] == self.molecule.name.replace(' ', '_'):
                i+=1
        i = f'{i:04}'
            
        name = self.molecule.name.replace(' ', '_') + '-' + i
        mkdir('Outputs/'+name)

        Create_rot_output([self.down_electrode, self.molecule, self.up_electrode], 
                          'Outputs/'+name, name, angles, boolean_data, SOC, self.opt_data, self.DFT_opts)
        return
    
    def upt_dft_chk(self): self.update_DFT(None)
    def update_DFT(self, _):
        self.DFT_opts['Functional'] = self.combo_funct.get()
        self.DFT_opts['Shell'] = self.combo_shell.get()
        self.DFT_opts['EABS'] = self.EABS.get()

        self.DFT_opts['DFTU'] = self.DFT_U.get()
        if self.DFT_opts['Functional'] not in self.molecule.U.keys():
            self.check_U.deselect()
            self.check_U['state'] = DISABLED
            self.DFT_opts['DFTU'] = False
        if self.DFT_opts['Functional'] in self.molecule.U.keys(): 
            self.check_U['state'] = NORMAL
        return

    def Update(self):
        ax.clear()
        
        change_flag = False
        
        if self.molecule.name != self.combo.get():
            self.mol_controls.destroy()
            self.molecule = xyz_object('Molecules/'+self.combo.get()+'.xyz')
            self.mol_controls = object_controls(self.molecule, self.frame_3)
            
            change_flag = True
                    
        if self.up_electrode.elements[0] != self.combo_top_elem.get():
            self.up_electrode.change_elem(self.combo_top_elem.get())
            
            change_flag = True
        
        if self.up_electrode.name != self.combo_top_type.get():
            self.top_controls.destroy()
            self.up_electrode = xyz_object('ElectrodesDFT/'+self.combo_top_type.get()+'.xyz')
            self.top_controls = object_controls(self.up_electrode, self.frame_0)

        if self.down_electrode.elements[0] != self.combo_bottom_elem.get():
            self.down_electrode.change_elem(self.combo_bottom_elem.get())
            
            change_flag = True

        if self.down_electrode.name != self.combo_bottom_type.get():
            self.down_electrode = xyz_object('ElectrodesDFT/'+self.combo_bottom_type.get()+'.xyz')
            #self.down_electrode.euler_angles = np.array([0,180, 0])
            self.down_electrode.reflection_z()
            self.down_electrode.spatial_coords = np.array([0,0,-11])

                    
        if change_flag and SOC_FLAG:
            if not self.up_electrode.elements[0] in SOC_PARAMS.keys() or not self.down_electrode.elements[0] in SOC_PARAMS.keys():
                self.soc_btn['state'] = DISABLED

            elif self.up_electrode.elements[0] in SOC_PARAMS.keys() and self.down_electrode.elements[0] in SOC_PARAMS.keys():
                self.soc_btn['state'] = NORMAL

        self.update_DFT(None)
        self.mol_controls.update()
        self.top_controls.update()
        
        self.down_electrode.plot_object(ax)
        
        set_axes_equal(ax)
        fig.canvas.draw()
        fig.canvas.blit()
        fig.canvas.flush_events()
        self.after_process = root.after(500, self.Update) 

#%% MAIN

def main():
    global root, ax, fig
    plt.ion()
    fig = plt.figure(figsize = (8,8))
    fig.patch.set_facecolor(bg_color)

    ax = plt.subplot(projection='3d')

    root = customtkinter.CTk()
    root.title(f'ANT.UI v{__version__}')
    root.geometry('720x720')

    a = tkinter_window()
    a.Update()
    root.mainloop()
    return

def call_back_dft(rt, axes, figure):
    global root, ax, fig
    
    root = rt
    ax = axes
    fig = figure
    
    a = tkinter_window()
    a.Update()
    return a

    

if __name__ == '__main__': main()



