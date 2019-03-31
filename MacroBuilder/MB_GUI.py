from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from MacroB import build_macrocomplex, unify_ids, read_pdbs, get_interaction_dict, update_interactions_dict, main_loop, save_results
from textwrap import dedent
import sys
import os
from subprocess import Popen, PIPE
from queue import Queue, Empty # Python 3
import threading
import subprocess

class StdRedirector():
    def __init__(self, text_widget):
        self.text_space = text_widget

    def write(self, string):
        self.text_space.config(state=NORMAL)
        self.text_space.insert("end", string)
        self.text_space.see("end")
        self.text_space.config(state=DISABLED)

class MB(Frame):
    def show_help(self):
        messagebox.showinfo(title="Help", message="Here I should show the document with the help of the program")

    def show_about(self):
        messagebox.showinfo(title="About", message="This program has been done in the PYT subject")
    def quit(self):
        if messagebox.askyesno("Quit", "Are you sure you want to exit?"):
            Frame.quit(self)

    def get_dir(self):
        directory.set(filedialog.askdirectory())

    def get_template(self):
        template_path.set(filedialog.askopenfilename(title="Select a template structure file", filetypes=[
            ("PDB files", "*.pdb"), ("mmCIF files", "*.cif")]))

    def create_menu(self):
        self.menubar = Menu(self)

        # CREATE THE FILEMENU
        filemenu = Menu(self.menubar)
        filemenu.add_command(label="Select Directory", command=self.get_dir)
        filemenu.add_separator()
        filemenu.add_command(label="QUIT", command=self.quit)

        # CREATE THE HELP MENU
        helpmenu = Menu(self.menubar)
        helpmenu.add_command(label="Help", command=self.show_help)
        helpmenu.add_command(label="About", command=self.show_about)

        self.menubar.add_cascade(label="File", menu=filemenu)
        self.menubar.add_cascade(label="Help", menu=helpmenu)

        self.master.config(menu=self.menubar)

    def create_options(self):
        self.options = LabelFrame(self, text="Options", padx=5, pady=5)
        self.create_options_frame()
        self.options.grid(row=0, column=0)

    def create_options_frame(self):
        frame = Frame(self.options)
        label_output = Label(frame, text="Output name")
        label_num_chains = Label(frame, text="Max number of chains")
        label_num_models = Label(frame, text="Number of models")
        label_current_dir = Label(frame, text="Selected directory:")
        label_current_dir_path = Label(frame, textvariable=directory)
        self.label_dirty = Checkbutton(frame, text="dirty mode (generates tmp files)", onvalue=True, offvalue=False, variable=dirty)
        self.label_verbose = Checkbutton(frame, text="Verbose (Show all steps)", onvalue=True, offvalue=False, variable=verbose)
        label_template = Label(frame, text="Select template (Optional)")
        entry_template = Button(frame, text="Browse", command=self.get_template)
        label_template_path = Label(frame, textvariable=template_path)
        entry_output = Entry(frame, textvariable=output_name)
        self.entry_max_chains = Spinbox(frame, from_=100, to=1000)
        self.entry_num_models = Spinbox(frame, from_=1, to=100)
        self.entry_run = Button(frame, text="RUN", command=self.thread_MB)
        label_current_dir.grid(row=0, column=0, sticky="w")
        label_current_dir_path.grid(row=0, column=1, sticky="w", columnspan=6)
        label_output.grid(row=1, column=0, sticky="w")
        entry_output.grid(row=1, column=1, sticky="w")
        label_num_chains.grid(row=2, column=0, sticky="w")
        self.entry_max_chains.grid(row=2, column=1, sticky="w")
        label_num_models.grid(row=3, column=0, sticky="w")
        self.entry_num_models.grid(row=3, column=1, sticky="w")
        self.label_dirty.grid(row=1, column=2, sticky="w", columnspan=2)
        self.label_verbose.grid(row=2, column=2, sticky="w", columnspan=2)
        label_template.grid(row=4, column=0, sticky="w")
        entry_template.grid(row=4, column=1, sticky="w")
        label_template_path.grid(row=4, column=2, columnspan=8)
        self.entry_run.grid(row=2, column=5, sticky="w")
        frame.grid(row=0)

    def create_sequence_dict(self):
        self.seq_dict_frame = LabelFrame(self, text="Sequences", padx=5, pady=5)
        self.create_sequence_dict_frame()
        self.seq_dict_frame.grid(row=1, column=0)

    def create_sequence_dict_frame(self):
        seq_frame = Frame(self.seq_dict_frame)
        scrollbar_v = Scrollbar(seq_frame, orient=VERTICAL)
        scrollbar_h = Scrollbar(seq_frame, orient=HORIZONTAL)
        self.seq_listbox = Text(seq_frame, yscrollcommand=scrollbar_v.set, xscrollcommand=scrollbar_h.set, width=70,
                                height=6, state=DISABLED, wrap="none", background="black", foreground="white")
        scrollbar_v.config(command=self.seq_listbox.yview)
        scrollbar_h.config(command=self.seq_listbox.xview)
        scrollbar_v.pack(side=RIGHT, fill=Y)
        scrollbar_h.pack(side=BOTTOM, fill=X)
        self.seq_listbox.pack(side=LEFT, expand=True, fill=BOTH)
        seq_frame.grid(row=0)

    def create_estequiometry(self):
        self.estequiometry_frame = LabelFrame(self, text="Estequiometry", padx=5, pady=5)
        self.create_estequiometry_frame()
        self.estequiometry_frame.grid(row=2, column=0)

    def create_estequiometry_frame(self):
        frame = Frame(self.estequiometry_frame)
        scrollbar = Scrollbar(frame, orient=HORIZONTAL)
        self.estequiometry_listbox = Text(frame, xscrollcommand=scrollbar.set, width=70, height=2,
                                          state=DISABLED, background="black", foreground="white")
        # scrollbar.config(command=self.seq_listbox.yview)
        scrollbar.pack(side=BOTTOM, fill=X)
        self.estequiometry_listbox.pack(side=LEFT, expand=True, fill=BOTH)
        frame.grid(row=0)

    def create_console(self):
        self.console_frame = LabelFrame(self, text="Console", padx=5, pady=5)
        self.create_console_frame()
        self.console_frame.grid(row=0, column=1, rowspan=5)

    def create_console_frame(self):
        frame = Frame(self.console_frame)
        scrollbar = Scrollbar(frame, orient=VERTICAL)
        self.console = Text(frame, yscrollcommand=scrollbar.set, width=60, height=50, state=DISABLED, wrap='word',
                            background="black", foreground="white")
        scrollbar.pack(side=RIGHT, fill=Y)
        scrollbar.config(command=self.console.yview)
        self.console.pack(side=LEFT, expand=True, fill=BOTH)
        sys.stdout = StdRedirector(self.console)
        sys.stderr = StdRedirector(self.console)
        frame.grid(row=0)

    def update_image(self, best_model):
        best_model_name = output_name.get()+"_1.cif"
        os.system("pymol %s -c -d 'hide all;show ribbon;util.cbc' -g /tmp/tmp.png " % best_model_name)
        #subprocess.run(["pymol", best_model_name, "-c", "-g", "tmp.png"])

        image = PhotoImage(file="/tmp/tmp.png")
        #self.model_image.config(image=image)
        #self.model_image.image = image

        self.model_image.create_image(300, 250, anchor=CENTER, image=image)
        self.model_image.image = image


        #self.model_image.configure(image=image, width=80, height=30)


    def create_image(self):
        self.image_frame = LabelFrame(self, text="Structure", padx=5, pady=5)
        self.create_image_frame()
        self.image_frame.grid(row=3, column=0)

    def create_image_frame(self):
        frame = Frame(self.image_frame)
        self.model_image = Canvas(frame, width=600, height=500, background="black")
        self.model_image.pack()
        frame.grid(row=0, column=1)


    def update_estequiometry(self, best_model):
        stq_dict = {}
        for chain in best_model:
            stq_dict.setdefault(chain.id, 0)
            stq_dict[chain.id] += 1
        self.estequiometry_listbox.config(state=NORMAL)
        self.estequiometry_listbox.delete(1.0, END)
        for key in sorted(stq_dict.keys()):
            self.estequiometry_listbox.insert("end", "%s:%s " % (key, stq_dict[key]))
        self.estequiometry_listbox.config(state=DISABLED)


    def run_MB(self):
        self.in_pdbmodels = read_pdbs(directory.get()+"/", verbose.get())
        self.seq_dict = unify_ids(self.in_pdbmodels)
        self.interaction_dict = get_interaction_dict(self.in_pdbmodels, verbose.get())
        self.update_seq_dict()
        update_interactions_dict(self.interaction_dict, verbose)
        out_pdbmodels = main_loop(int(self.entry_num_models.get()), output_name.get(), self.seq_dict, self.interaction_dict,
                  verbose=verbose.get(), max_chains=int(self.entry_max_chains.get()), dirty=dirty.get(), template=template_path.get())
        self.update_estequiometry(out_pdbmodels[0])
        save_results(out_pdbmodels, output_name.get(), self.seq_dict, template=template_path.get(), verbose=verbose.get())
        self.update_image(out_pdbmodels[0])
        self.entry_run.config(state='normal')


    def thread_MB(self):
        self.console.config(state=NORMAL)
        self.console.delete(1.0, END)
        self.console.config(state=DISABLED)
        self.estequiometry_listbox.config(state=NORMAL)
        self.estequiometry_listbox.delete(1.0, END)
        self.estequiometry_listbox.config(state=DISABLED)
        self.seq_listbox.config(state=NORMAL)
        self.seq_listbox.delete(1.0, END)
        self.seq_listbox.config(state=DISABLED)
        self.model_image.delete("all")
        self.entry_run.config(state='disabled')
        #q = Queue()  # limit output buffering (may stall subprocess)
        t = threading.Thread(target=self.run_MB)
        t.daemon = True # close pipe if GUI process exits
        t.start()



    def createWidgets(self):
        self.create_menu()
        self.create_console()
        self.create_options()
        self.create_sequence_dict()
        self.create_estequiometry()

        self.create_image()
        self.grid(row=0)

    def update_seq_dict(self):
        self.seq_listbox.config(state=NORMAL)
        self.seq_listbox.delete(1.0, END)
        for key in self.seq_dict:
            self.seq_listbox.insert("end", "%s: %s\n" % (self.seq_dict[key], key))
        self.seq_listbox.config(state=DISABLED)



    def __init__(self, master=None, **kwargs):
        Frame.__init__(self, master, **kwargs)
        self.master.wm_title("Macrocomplex Builder")
        self.master.resizable(width=False, height=False)
        self.createWidgets()


if not os.path.exists('tmp'):
    os.mkdir('tmp')
if not os.path.exists('models'):
    os.mkdir('models')
root = Tk()
directory = StringVar()
template_path = StringVar()
output_name = StringVar()
max_number_chains = IntVar()
verbose = BooleanVar()
dirty = BooleanVar()
directory.set('None')
app = MB(master=root, padx=10, pady=10)
root.mainloop()