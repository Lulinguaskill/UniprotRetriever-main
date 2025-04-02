"""  File for Uniprot UI"""
import tkinter

import customtkinter as ctk  # Wrapper for tkinter
from PIL import Image  # Image manipulation module
import os  # Who knows ?
from pyperclip import copy  # Clipboard manager
import threading
from functools import partial
from tktooltip import ToolTip
from PIL import Image, ImageTk
from tkinter import filedialog

from RequestTools import *  # Our request functions !
from BasicClasses import *  # Our objects
from utilities import *  # Bazar
from CTkGif import *
from blastUtils import get_blast_result, blast
from inflexion_points import *

ctk.set_appearance_mode("light")  # select light theme
themepath = os.path.join('ressources', "custom-theme.json")  # path to json theme file
ctk.set_default_color_theme(themepath)  # theme application
active_protein = None

root = None  # application

font = ('AvenirNextLTPro-Bold', 20)  # Font for text

colors = dict(blue='#006F9D', hblue='#1E93C3', dark_blue='#03045E', white='#FBFBFE', nwhite='#DEDCFF', red='#DD534A',
              green='#6EBB3F', light_blue='#9AE1FF', gray='#D9D9D9', hdb='#080A9D')

test_protein = 'P53601'

options = None

parameters = [('Database access path :', 'db_path'),
              ('CPU cores for blast : ', 'cpu')]


class MainApp(ctk.CTk):
    """ root for UniprotRetriever"""

    def __init__(self, **kwargs):
        global root
        super().__init__(**kwargs)  # Calling parent
        self.title('Uniprot Retriever')
        self.geometry(
            center_window(972, 648, int(self.winfo_screenwidth()), int(self.winfo_screenheight())))  # centering
        self.minsize(432, 520)
        self.iconbitmap('ressources/icon.ico')  # Logo of the window
        self.mainFrame = MainFrame(self)  # Plotting our first frame
        self.prev_frame = None
        root = self  # Storing self in the root
        self.display()  # Display
        self.mainloop()  # Running

    def display(self):
        self.grid_rowconfigure(0, weight=1)  # Centering
        self.grid_columnconfigure(0, weight=1)
        self.mainFrame.grid(row=0, column=0, padx=0, pady=0, sticky='nsew')

    def go_menu(self, target_frame, option=None, asynchronous_load=False):
        """ remove current active frame and place the with """
        if asynchronous_load:
            threading.Thread(target=self.go_menu(target_frame, option, asynchronous_load=False)).start()
        else:
            self.prev_frame = self.mainFrame
            if option is not None:
                self.mainFrame = target_frame(self, option)
            else:
                self.mainFrame = target_frame(self)
            self.mainFrame.grid_remove()
            self.mainFrame.grid(row=0, column=0, padx=0, pady=0, sticky='nsew')

    def home(self):
        self.go_menu(MainFrame)

    def prev(self):
        """ Remove current frame and goes to the previous one (if changed with go_menu)"""
        self.mainFrame.grid_remove()
        self.mainFrame = self.prev_frame
        self.mainFrame.grid(row=0, column=0, padx=0, pady=0, sticky='nsew')


class TaxonomyFrame(ctk.CTkFrame):
    """ Frame for displaying the whole taxonomy of a specie """

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        active_protein.build_taxonomy()
        self.taxo = active_protein.taxonomy.get_taxonomy()
        print(self.taxo)
        self.infoFrame = TaxonomyRankFrame(self, self.taxo[-1], fg_color=colors['nwhite'], corner_radius=10)
        self.treeFrame = TaxonomyTreeFrame(self, self.taxo, fg_color=colors['nwhite'], corner_radius=10)
        self.homeButton = IconButton(self, 'return', command=root.prev)
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=0)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure((0, 1), weight=1)
        self.homeButton.grid(row=0, column=0, padx=(10, 0), pady=(10, 0), sticky='w')
        self.treeFrame.grid(row=1, column=0, padx=15, pady=15, sticky='nsew')
        self.infoFrame.grid(row=1, column=1, padx=15, pady=15, sticky='nsew')


class TaxonomyTreeFrame(ctk.CTkFrame):
    """ Frame for displaying the taxonomy tree of a specie """

    def __init__(self, master: any, taxo: list[str], **kwargs):
        super().__init__(master, **kwargs)
        self.taxonomy = taxo  # List of species
        self.nb_specie = len(self.taxonomy)  # Number of species in taxonomy
        self.buttons = {}
        self.gradient = get_color_gradient('#03045E', '#1E93C3', self.nb_specie)  # Buttons colors
        self.build_buttons()

    def build_buttons(self):
        offset = 4
        button_height = 25
        additionnal_padding = 10
        self.grid_columnconfigure((0, 1), weight=0)
        for index, name in enumerate(self.taxonomy):
            if index == 0:
                ymorepad = 10
            else:
                ymorepad = 0
            self.buttons[index] = ctk.CTkButton(self, height=button_height, width=250, text=name, font=font,
                                                border_width=2, fg_color=self.gradient[index], border_color='#000000',
                                                text_color='#FFFFFF', hover_color=get_hover_color(self.gradient[index]),
                                                command=partial(self.master.infoFrame.update_labels_callback, name))
            self.buttons[index].grid(row=index, column=0, padx=((index + 1) * offset + additionnal_padding, 0),
                                     pady=(3 + ymorepad, 3), sticky='w')


class TaxonomyRankFrame(ctk.CTkFrame):
    """ Frame for displaying data about a specific rank in a specie taxonomy """

    def __init__(self, master: any, rank: str, **kwargs):
        super().__init__(master, **kwargs)
        self.rank = rank  # Taxonomy rank
        self.descriptions = None
        self.title = ctk.CTkLabel(self, text=self.rank, font=('Avenir LT 65 Medium Bold', 35), justify='left',
                                  text_color='#000000')  # Name of the rank for display
        self.descrTextBox = ctk.CTkTextbox(self, font=font, corner_radius=5, text_color='#000000',
                                           border_width=2, border_color='#000000',
                                           wrap='word')  # Container for the descritpion
        self.update_labels_callback(self.rank)  # Callback in subprocess for faster first load
        self.bar = ctk.CTkLabel(self, text='',
                                image=ctk.CTkImage(Image.open('ressources/bar.png'), size=(500, 5)))  # Sep bar
        self.display()
        self.build_descriptions()

    def display(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure((0, 1), weight=0)
        self.grid_rowconfigure(2, weight=1)
        self.title.grid(row=0, column=0, pady=(15, 0), padx=10, sticky='s')
        self.bar.grid(row=1, column=0, pady=0, padx=10, sticky='ew')
        self.descrTextBox.grid(row=2, column=0, pady=(0, 15), padx=15, sticky='nsew')

    def update_labels(self, new_rank, mh=None):
        """ This is called if a new button is clicked before this class ended to download all the descriptions"""
        self.rank = new_rank
        self.title.configure(text=self.rank)
        self.descrTextBox.delete("0.0", "end")
        self.descrTextBox.insert("0.0", get_rank_data(new_rank))
        self.lock_entry()

    def update_labels_callback(self, new_rank):
        self.descrTextBox.configure(state='normal')
        self.descrTextBox.delete("0.0", "end")
        if new_rank is None:
            self.descrTextBox.insert("0.0", "No description available")
        else:
            self.rank = new_rank
            self.title.configure(text=self.rank)
            if self.descriptions is not None and self.rank in self.descriptions.keys():
                self.descrTextBox.insert("0.0", self.descriptions[self.rank])
                self.lock_entry()
            else:
                threading.Thread(target=partial(self.update_labels, self.rank)).start()

    def lock_entry(self):
        self.descrTextBox.configure(state='disabled')

    def build_descriptions(self):
        threading.Thread(target=partial(self.build_descriptions_callback, self.master.taxo)).start()

    def build_descriptions_callback(self, taxo):
        self.descriptions = get_taxonomy_data(taxo)
        active_protein.taxonomy.save_descr(self.descriptions)


class MainFrame(ctk.CTkFrame):
    """ Welcome page"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        ratio = 0.7  # Logo size ratio
        self.logo = ctk.CTkLabel(self, text='', image=ctk.CTkImage(Image.open(os.path.join('ressources', 'logo.png')),
                                                                   size=(
                                                                       496 * ratio,
                                                                       227 * ratio)))  # Empty label with image
        self.select = SelectFrame(self, fg_color='transparent')
        self.search = SearchFrame(self, fg_color='transparent')
        self.load = LoadingAnimationFrame(self, fg_color='transparent')
        self.orLabel = ctk.CTkLabel(self, text='Or', text_color=colors['dark_blue'], font=font)
        self.paramButton = IconButton(self, 'param', command=self.disp_param)
        self.dbButton = IconButton(self, 'db')
        self.display()

    def display(self):
        """ Configure the layout and place widgets on the window"""
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((0, 2), weight=0)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure((1, 2, 3, 4), weight=1)
        self.paramButton.grid(row=0, column=0, padx=8, pady=8, sticky='nw')
        self.logo.grid(row=0, column=1, padx=0, pady=8, sticky='n')
        self.dbButton.grid(row=0, column=2, padx=8, pady=8, sticky='ne')
        self.search.grid(row=1, column=1, padx=0, pady=0, sticky='ew')
        self.orLabel.grid(row=2, column=1, padx=0, pady=0, sticky='')
        self.select.grid(row=3, column=1, padx=0, pady=0, sticky='ew')
        self.load.grid(row=4, column=1, padx=0, pady=0, sticky='ew')

    @staticmethod
    def disp_param():
        root.go_menu(OptionsFrame)


class SearchFrame(ctk.CTkFrame):
    """ frame to enter and search an ID on Uniprot"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.label = ctk.CTkLabel(self, text='Enter an Uniprot ID : ', text_color=colors['dark_blue'], font=font)
        self.entry = ctk.CTkEntry(self, text_color=colors['dark_blue'], font=font, border_width=0, height=36,
                                  fg_color=colors['gray'], width=305, justify='center')
        self.entry.insert(0, 'Q9NZT1')
        self.button = ctk.CTkButton(self, text_color=colors['gray'], text='Search',
                                    command=threading.Thread(target=self.button_callback).start,
                                    corner_radius=5, fg_color=colors['dark_blue'], width=100, font=font, height=36,
                                    hover_color=colors['hdb'])
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure((0, 1, 2), weight=1)
        self.label.grid(row=0, column=0, padx=0, pady=0, sticky='w')
        self.entry.grid(row=0, column=1, padx=0, pady=0, sticky='w')
        self.button.grid(row=0, column=2, padx=(30, 0), pady=0, sticky='e')

    def button_callback(self, event=None):
        global active_protein
        entry = self.entry.get()  # Uniprot ID
        active_protein = Protein(entry, self.master.load)  # Create protein class and retrieve necessary information
        bind_options()  # Bind the options of the next frame to new active protein
        root.go_menu(ProteinDashBoardFrame, option=None, asynchronous_load=True)  # Change UI
        # active_protein.get_2D_render()


class SelectFrame(ctk.CTkFrame):
    """ frame to select an ID from local DB"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.available_locals = get_local_proteins()
        self.label = ctk.CTkLabel(self, text='Select from local : ', text_color=colors['dark_blue'], font=font)
        self.select = ctk.CTkOptionMenu(self, text_color=colors['dark_blue'], font=font, corner_radius=5,
                                        fg_color=colors['gray'], values=self.available_locals,
                                        button_color=colors['dark_blue'], button_hover_color=colors['hdb'],
                                        width=305, height=36)
        self.select._canvas.itemconfig("dropdown_arrow", fill=colors['gray'])
        self.select.set('')
        self.button = ctk.CTkButton(self, text_color=colors['gray'], text='Load', command=self.button_callback,
                                    corner_radius=5, fg_color=colors['dark_blue'], width=100, font=font, height=36,
                                    hover_color=colors['hdb'])
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure((0, 1, 2), weight=1)
        self.label.grid(row=0, column=0, padx=0, pady=0, sticky='w')
        self.select.grid(row=0, column=1, padx=0, pady=0, sticky='w')
        self.button.grid(row=0, column=2, padx=0, pady=0, sticky='e')

    def button_callback(self, event=None):
        global active_protein
        entry = self.select.get()  # Uniprot ID
        active_protein = Protein(entry, ui=self.master.load)  # Create protein class and retrieve necessary information
        bind_options()  # Bind the options of the next frame to new active protein
        root.go_menu(ProteinDashBoardFrame, option=None, asynchronous_load=True)  # Change UI
        # active_protein.get_2D_render()


class LoadingAnimationFrame(ctk.CTkFrame):
    """ frame to display loading steps """

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.steps = ['Creating local repository ...', 'Retrieving fasta file from Uniprot ...',
                      'Extracting data from fasta file ...']
        self.loading_bar = ctk.CTkProgressBar(self, progress_color=colors['hdb'], fg_color=colors['gray'],
                                              border_width=3, height=22,
                                              border_color=colors['dark_blue'], width=762)
        self.label = ctk.CTkLabel(self, text='Waiting for user input ...', text_color=colors['dark_blue'],
                                  font=font)
        self.loading_bar.set(0)
        self.index = 0
        self.increm = 0
        self.distance = 1 / len(self.steps)
        self.display()

    def display(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure((0, 1), weight=1)
        self.loading_bar.grid(row=0, column=0, padx=0, pady=25, sticky='')
        self.label.grid(row=1, column=0, padx=0, pady=0, sticky='')

    def progress(self):
        try:
            self.label.configure(text=self.steps[self.index])
            self.index += 1
            self.after(50, self.step)
        except IndexError:
            pass

    def step(self):
        if self.increm < 20:
            if self.loading_bar.get() < self.distance * self.index:
                self.loading_bar.set(self.distance * (self.index - 1) + self.increm * self.distance / 20)
                self.increm += 1
                self.after(50, self.step)
        else:
            self.loading_bar.set(self.distance * self.index)
            self.increm = 0


class ProteinDashBoardFrame(ctk.CTkFrame):
    """ frame to display all the information about a protein and to select actions to perform"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.infoFrame = ProteinInfoFrame(self, fg_color=colors['gray'])
        self.optionFrame = ProteinOptionsFrame(self, fg_color=colors['gray'])
        self.renderFrame = None
        self.homeButton = IconButton(self, 'home', root.home)
        self.renderFrame = ImageFrame(self, [220, 220], task=Task(active_protein.build_render), fg_color=colors['gray'])
        self.photoFrame = ImageFrame(self, [220, 220], task=Task(active_protein.download_photo),
                                     fg_color=colors['gray'])
        self.title = ctk.CTkLabel(self, text=active_protein.prot_id, text_color=colors['dark_blue'], font=(font[0], 60))
        self.title.bind('<Button-1>', active_protein.open_uniprot_webpage)
        self.bar = ctk.CTkLabel(self, text='',
                                image=ctk.CTkImage(Image.open('ressources/large_bar.png'), size=(776, 10)))  # Sep bar
        self.display()

    def display(self):
        self.grid_rowconfigure((0, 1), weight=0)
        self.grid_rowconfigure((2, 3), weight=1)
        self.grid_columnconfigure((0, 1, 2), weight=1)
        self.homeButton.grid(row=0, column=0, padx=(20, 0), pady=(10, 0), sticky='nw')
        self.title.grid(row=0, column=0, padx=0, pady=0, sticky='', columnspan=3)
        self.bar.grid(row=1, column=0, padx=0, pady=0, sticky='', columnspan=3)
        self.infoFrame.grid(row=2, column=0, padx=(20, 0), pady=20, sticky='nsew', columnspan=2)
        self.optionFrame.grid(row=2, column=2, padx=20, pady=20, sticky='nsew', rowspan=2)
        self.renderFrame.grid(row=3, column=0, padx=20, pady=(0, 20), sticky='nsew')
        self.photoFrame.grid(row=3, column=1, padx=0, pady=(0, 20), sticky='nsew')


class ProteinInfoFrame(ctk.CTkFrame):
    """ frame to display basic info about the protein"""

    def __init__(self, master: any, cjustify='left', **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.cjustify = cjustify  # custom justification
        self.data = active_protein.get_info()
        self.names = ['Protein', 'Organism', 'Amino acids']
        self.labels = concatenate_infos(self.names, self.data)
        self.labels_widget = []
        self.build_widgets()

    def build_widgets(self):
        self.grid_columnconfigure(0, weight=1)
        for index, label in enumerate(self.labels):
            self.labels_widget.append(
                ctk.CTkLabel(self, text=str(label), text_color=colors['dark_blue'], font=font, justify=self.cjustify))
            self.grid_rowconfigure(index, weight=1)
            self.labels_widget[index].grid(row=index, column=0, padx=(20, 0), pady=0, sticky='w')


class ProteinOptionsFrame(ctk.CTkFrame):
    """ frame to display the possible actions to perform on the active protein"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.iconFrame = OptionsIconFrame(self, fg_color='transparent')
        self.optionLabel = ctk.CTkLabel(self, text=self.iconFrame.active_option, text_color=colors['dark_blue'],
                                        font=(font[0], 30))
        self.describeLabel = ctk.CTkLabel(self, text=options[self.iconFrame.active_option][0],
                                          text_color=colors['dark_blue'], font=font,
                                          justify='center', width=250)
        self.button = ctk.CTkButton(self, text_color=colors['gray'], text='Select', command=self.button_callback,
                                    corner_radius=5, fg_color=colors['dark_blue'], width=100, font=font, height=36,
                                    hover_color=colors['hdb'])
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=0)
        self.grid_rowconfigure((1, 2, 3), weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.iconFrame.grid(row=0, column=0, padx=0, pady=0, sticky='')
        self.optionLabel.grid(row=1, column=0, padx=0, pady=0, sticky='')
        self.describeLabel.grid(row=2, column=0, padx=0, pady=0, sticky='')
        self.button.grid(row=3, column=0, padx=0, pady=0, sticky='')

    def button_callback(self):
        func = options[self.iconFrame.active_option][2]
        if func is not None:
            func()

    def set_active_option(self, option_name):
        self.optionLabel.configure(text=option_name)
        self.describeLabel.configure(text=options[option_name][0])

    @staticmethod
    def disp_2D_render():
        root.go_menu(MatplotlibRenderFrame, option=active_protein.get_2D_render, asynchronous_load=True)

    @staticmethod
    def disp_taxo():
        root.go_menu(TaxonomyFrame)

    @staticmethod
    def disp_blast():
        root.go_menu(BlastFrame)


class OptionsIconFrame(ctk.CTkFrame):
    """ icons for options"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.buttons = {}
        self.active_option = 'Save'  # key of the option dictionnary of current selected option
        self.build_buttons()

    def build_buttons(self):
        self.grid_rowconfigure((0, 1), weight=1)
        self.grid_columnconfigure((0, 1, 2), weight=1)
        i = 0
        for name, params in options.items():
            self.buttons[name] = IconButton(self, file=params[1], width=60)
            self.buttons[name].bind('<Button-1>', self.button_callback)
            self.buttons[name].grid(row=1 if i > 2 else 0, column=i % 3, padx=25, pady=25, sticky='')
            i += 1
        self.toggle(self.active_option)

    def button_callback(self, event: tkinter.Event):
        name_pressed = find_key_by_value(self.buttons, event.widget.master)  # finding the clicked button
        if name_pressed == self.active_option:
            self.master.button_callback()
        else:
            self.toggle(self.active_option)
            self.toggle(name_pressed)
            self.active_option = name_pressed
            self.master.set_active_option(self.active_option)

    def toggle(self, name):
        button = self.buttons[name]
        if button.cget('fg_color') == colors['hblue']:  # on
            button.configure(fg_color=colors['white'])
        else:
            button.configure(fg_color=colors['hblue'])


class ImageFrame(ctk.CTkFrame):
    """ class for displaying images inside a frame. Adapts to whatever image size"""

    def __init__(self, master: any, constraints, task: Task, autostart=True, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.constraints = constraints
        self.gif = CTkGif(self, 'ressources/loading.gif', loop=True, acceleration=3, size=(110, 110))
        self.task = task
        self.image_path = None
        self.auto_start = autostart
        self.img = ctk.CTkLabel(self, text='')
        self.img.bind('<Button-1>', self.open_img)
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.gif.grid(row=0, column=0, pady=0, padx=0, sticky='')
        if self.auto_start:
            self.gif.start()
            self.start_image_render()

    def display_img(self, image_path):
        self.image_path = image_path
        new_size = get_constrained_image_dimensions(self.constraints, image_path)
        self.img.configure(image=ctk.CTkImage(Image.open(image_path), size=new_size))
        self.gif.stop()
        self.gif.grid_remove()
        self.img.grid(row=0, column=0, pady=0, padx=0, sticky='')
        print('image displayed')

    def open_img(self, event=None):
        os.startfile(self.image_path)

    def start_image_render(self):
        if not self.auto_start:
            self.gif.start()
        print(f'tache : {self.task.target_func}')
        self.task.start(self)

    def end(self, output):
        print('end')
        self.display_img(output)


class IconButton(ctk.CTkButton):
    """ Custom button for small icon buttons"""

    def __init__(self, master: any, file: str, command=None, width=40, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.configure(width=width, height=width, text='', command=command, border_width=2.5,
                       border_color=colors['dark_blue'], fg_color=colors['white'], corner_radius=5,
                       image=ctk.CTkImage(Image.open(f'ressources/{file}.png'), size=(30, 30)))


# plt.savefig("plot.svg", format="svg")


class MatplotlibRenderFrame(ctk.CTkFrame):
    """ Frame for displaying a render done with mpl. Offers multiple possbilities and recall proteins info"""

    def __init__(self, master: any, renderfunc, size=(640, 480), autostart=True, info=True, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.info = info
        self.imageFrame = ImageFrame(self, size, Task(renderfunc), autostart, fg_color='transparent')
        self.optionsFrame = RenderIconFrame(self, fg_color='transparent',
                                            copymethod=active_protein.copy_2D_render_to_clipboard)
        self.topInfoFrame = ProteinInfoFrame(self, fg_color=colors['gray'], cjustify='center')
        self.display()

    def display(self):
        self.grid_rowconfigure((0, 1), weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=0)
        self.imageFrame.grid(row=1, column=0, padx=20, pady=(0, 20), sticky='nsew')
        self.optionsFrame.grid(row=1, column=1, padx=(0, 20), pady=(0, 20), sticky='nsew')
        if self.info:
            self.topInfoFrame.grid(row=0, column=0, padx=20, pady=20, sticky='nsew', columnspan=2)


class RenderIconFrame(ctk.CTkFrame):
    """ icons for options in rendering mpl frame"""

    def __init__(self, master: any, copymethod, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.options = [('copy', copymethod), ('save', self.save_img), ('return', root.prev)]
        self.buttons = {}
        self.build_buttons()

    def build_buttons(self):
        self.grid_rowconfigure([i for i in range(len(options))], weight=1)
        self.grid_columnconfigure(0, weight=1)
        for index, (icon, func) in enumerate(self.options):
            self.buttons[index] = IconButton(self, file=icon, width=60, command=func)
            self.buttons[index].grid(row=index, column=0, padx=25, pady=25, sticky='')

    @staticmethod
    def save_img():
        output = filedialog.asksaveasfilename(defaultextension='png')
        if output is not None:
            active_protein.save_as(output)


class BlastFrame(ctk.CTkFrame):
    """ Frame containing all the display tools for blast, including sequence and the graph"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.sequenceFrame = SequenceFrame(self, fg_color=colors['gray'])
        self.blastOptionsFrame = BlastOptionsFrame(self, fg_color=colors['gray'])
        self.resultFrame = BlastResultFrame(self, fg_color=colors['gray'])
        self.display()

    def display(self):
        self.grid_rowconfigure((0, 1), weight=1)
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.sequenceFrame.grid(row=0, column=0, padx=20, pady=20, sticky='ew', columnspan=2)
        self.blastOptionsFrame.grid(row=1, column=0, padx=20, pady=(0, 20), sticky='nsw')
        self.resultFrame.grid(row=1, column=1, padx=(0, 20), pady=(0, 20), sticky='nsew')


class SequenceFrame(ctk.CTkFrame):
    """ frame to display the raw sequence and copy it"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.label = ctk.CTkLabel(self, text='Sequence', text_color=colors['dark_blue'], font=font)
        self.entry = ctk.CTkEntry(self, text_color=colors['dark_blue'], font=font, border_width=0, height=40,
                                  fg_color=colors['white'], width=305, justify='left')
        self.entry.insert(0, str(active_protein.sequence))
        self.copyButton = IconButton(self, 'copy', self.copy_sequence)
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=0)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=0)
        self.label.grid(row=0, column=0, pady=20, padx=25, sticky='w')
        self.entry.grid(row=1, column=0, pady=(0, 20), padx=20, sticky='ew')
        self.copyButton.grid(row=1, column=1, pady=(0, 20), padx=(0, 20), sticky='e')

    def copy_sequence(self):
        copy(self.entry.get())


class BlastOptionsFrame(ctk.CTkFrame):
    """ frame to display the raw sequence and copy it"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.pil_img_graph = None  # Storing the produced graph of first blast
        self.dbLabel = ctk.CTkLabel(self, text='Database', text_color=colors['dark_blue'], font=font)
        self.boundariesLabel = ctk.CTkLabel(self, text='Graph boundaries (power of 10)', text_color=colors['dark_blue'],
                                            font=font)
        self.toLabel = ctk.CTkLabel(self, text='to', text_color=colors['dark_blue'], font=font)
        self.databaseChoice = ctk.CTkOptionMenu(self, text_color=colors['dark_blue'], font=font, corner_radius=5,
                                                fg_color=colors['white'], values=['Swissprot', 'NR'],
                                                button_color=colors['dark_blue'], button_hover_color=colors['hdb'],
                                                width=320, height=36)
        self.downEntry = ctk.CTkEntry(self, text_color=colors['dark_blue'], font=font, border_width=0, height=40,
                                      fg_color=colors['white'], width=80, justify='center')
        self.upEntry = ctk.CTkEntry(self, text_color=colors['dark_blue'], font=font, border_width=0, height=40,
                                    fg_color=colors['white'], width=80, justify='center')
        self.downEntry.bind('<Return>', self.set_focus)
        self.databaseChoice._canvas.itemconfig("dropdown_arrow", fill=colors['gray'])
        self.goButton = ctk.CTkButton(self, text_color=colors['gray'], text='Generate graph',
                                      command=self.button_callback,
                                      corner_radius=5, fg_color=colors['dark_blue'], width=200, font=font, height=36,
                                      hover_color=colors['hdb'])
        self.display()

    def display(self):
        self.grid_rowconfigure((1, 3, 4), weight=1)
        self.grid_rowconfigure((0, 2), weight=0)
        self.grid_columnconfigure((0, 2), weight=1)
        self.grid_columnconfigure(1, weight=0)
        self.dbLabel.grid(row=0, column=0, pady=20, padx=20, sticky='sw')
        self.databaseChoice.grid(row=1, column=0, pady=20, padx=20, sticky='nw', columnspan=3)
        self.boundariesLabel.grid(row=2, column=0, pady=20, padx=20, sticky='w', columnspan=3)
        self.downEntry.grid(row=3, column=0, pady=20, padx=20, sticky='w')
        self.toLabel.grid(row=3, column=1, pady=20, padx=20, sticky='ew')
        self.upEntry.grid(row=3, column=2, pady=20, padx=20, sticky='e')
        self.goButton.grid(row=4, column=0, pady=20, padx=20, sticky='', columnspan=3)

    def button_callback(self):
        # up, down = int(self.downEntry.get()), int(self.upEntry.get())
        self.master.resultFrame.imageFrame.start_image_render()

    def set_focus(self, event=None):
        self.upEntry.focus_set()


class BlastResultFrame(ctk.CTkFrame):
    """ frame to display the info of a family"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.imageFrame = ImageFrame(self, constraints=(380, 380), task=Task(active_protein.blast), autostart=False,
                                     fg_color='transparent')
        self.optionsFrame = RenderIconFrame(self, fg_color='transparent',
                                            copymethod=active_protein.copy_blastrender_to_clipboard)
        self.recursiveButton = ctk.CTkButton(self, text_color=colors['gray'], text='Recursive blast on first family',
                                             command=self.go_recursive_blast,
                                             corner_radius=5, fg_color=colors['dark_blue'], width=200, font=font,
                                             height=36,
                                             hover_color=colors['hdb'])
        self.display()

    def display(self):
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=0)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=0)
        self.imageFrame.grid(row=0, column=0, padx=20, pady=20, sticky='nsew')
        self.optionsFrame.grid(row=0, column=1, padx=0, pady=20, sticky='nsew')
        self.recursiveButton.grid(row=1, column=0, padx=0, pady=(0, 20), sticky='', columnspan=2)

    @staticmethod
    def go_recursive_blast():
        root.go_menu(RecursiveBlastFrame, option=None)


class FamilyDashboardFrame(ctk.CTkFrame):
    """ frame to display the info of a family"""

    def __init__(self, master: any, family, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.family = family
        self.label = ctk.CTkLabel(self, text=f'Family size : {len(family.sequences)}, '
                                             f'Mean sequence lenght : {mean(family.get_aa_dispersion)}',
                                  text_color=colors['dark_blue'], font=font)


class RecursiveBlastFrame(ctk.CTkFrame):
    """ frame to display the info of a family"""

    def __init__(self, master: any, family,  **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.family = family
        self.title = ctk.CTkLabel(self, text=f'Recursive Blast on {active_protein.prot_id} first family',
                                  font=('Avenir LT 65 Medium Bold', 35), justify='left',
                                  text_color='#000000')
        self.bar = ctk.CTkLabel(self, text='',
                                image=ctk.CTkImage(Image.open('ressources/bar.png'), size=(500, 5)))  # Sep bar
        self.familyInfoFrame = FamilyDashboardFrame(self, self.family, fg_color=colors['gray'])
        self.rblastInfoFrame = ctk.CTkFrame(self, fg_color=colors['gray'])
        self.rblastGraphFrame = ctk.CTkFrame(self, fg_color=colors['gray'])
        self.waffleFrame = ImageFrame(self, (300, 300), Task(family.build_waffle), True, fg_color='transparent')
        self.display()

    def display(self):
        self.grid_columnconfigure((0, 1), weight=1)
        self.grid_rowconfigure((0, 1), weight=0)
        self.grid_rowconfigure((2, 3), weight=1)
        self.title.grid(row=0, column=0, padx=0, pady=0, sticky='ew', columnpsan=2)
        self.bar.grid(row=1, column=0, padx=0, pady=0, sticky='ew', columnpsan=2)
        self.familyInfoFrame.grid(row=2, column=0, padx=0, pady=0, sticky='nsew')
        self.rblastInfoFrame.grid(row=3, column=0, padx=0, pady=0, sticky='nsew')
        self.rblastGraphFrame.grid(row=2, column=1, padx=0, pady=0, sticky='nsew')
        self.waffleFrame.grid(row=3, column=1, padx=0, pady=0, sticky='nsew')


class OptionsFrame(ctk.CTkFrame):
    """ parameters of the software"""

    def __init__(self, master: any, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.title = ctk.CTkLabel(self, text='Parameters', text_color=colors['dark_blue'], font=(font[0], 60))
        self.bar = ctk.CTkLabel(self, text='',
                                image=ctk.CTkImage(Image.open('ressources/large_bar.png'), size=(600, 15)))  # Sep bar
        self.homeButton = IconButton(self, 'home', command=root.home)
        self.lines = []
        self.display()

    def display(self):
        self.grid_rowconfigure((0, 1, 2), weight=0)
        self.grid_columnconfigure(0, weight=1)
        self.homeButton.grid(row=0, column=0, padx=(10, 0), pady=(10, 0), sticky='sw')
        self.title.grid(row=1, column=0, padx=0, pady=20, sticky='new')
        self.bar.grid(row=2, column=0, padx=0, pady=(0, 20), sticky='new')
        i = 0
        for label, key in parameters:
            self.grid_rowconfigure(i+3, weight=1)
            self.lines.append(ParameterLine(self, label, key, fg_color='transparent'))
            self.lines[-1].grid(row=i+3, column=0, padx=20, pady=(0, 20), sticky='nsew')
            i += 1



class ParameterLine(ctk.CTkFrame):
    """ line to change a parameter"""

    def __init__(self, master: any, label, config_key, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.label = label
        self.config = get_config('config.yml')
        self.config_key = config_key
        self.label = ctk.CTkLabel(self, text=label, text_color=colors['dark_blue'], font=font)
        self.entry = ctk.CTkEntry(self, text_color=colors['dark_blue'], font=font, border_width=0, height=36,
                                  fg_color=colors['gray'], width=305, justify='center')
        self.entry.insert(0, str(self.config[self.config_key]))
        self.entry.bind('<Return>', self.save)
        self.display()

    def display(self):
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.label.grid(row=0, column=0, padx=20, pady=20, sticky='')
        self.entry.grid(row=0, column=1, padx=(0, 20), pady=20, sticky='ew')

    def save(self, event=None):
        self.config[self.config_key] = self.entry.get()
        save_config(self.config, 'config.yml')
        self.green_border()

    def green_border(self):
        self.entry.configure(border_width=2, border_color=colors['green'])
        self.after(3000, self.clear_border)

    def clear_border(self):
        self.entry.configure(border_width=0, border_color=colors['green'])


def bind_options():
    """ bind the functions of IconButtons after an active_protein change"""
    global options
    options = dict(Save=('Save the current proteins \nand loaded informations\n in a local database', 'db', None),
                   Photo=('Display the photo of specie\n in bigger screen', 'photo', None),
                   Taxonomy=('Get the whole specie \ntaxonomy with treewiew\n and descriptions', 'treeview',
                             ProteinOptionsFrame.disp_taxo),
                   Structure=('Informations about 3D \nstructure of the molecule', 'aa', None),
                   Sequence=(
                       'Sequence tools such as \nBLAST and alignement', 'sequence', ProteinOptionsFrame.disp_blast),
                   Render=('Build 2D render : \nDij graph with \nsequence alignment', 'copy',
                           ProteinOptionsFrame.disp_2D_render))
