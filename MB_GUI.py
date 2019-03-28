import tkinter
import tkinter.messagebox
import tkinter.filedialog

entry_var = tkinter.StringVar()

class Application(tkinter.Frame):
    def open_directory(self):
        fasta_filename_path = tkinter.filedialog.askdirectory()

    def show_help(self):
        tkinter.messagebox.showinfo(title="Help",
                                    message="Here I should show the document with the help of the program")

    def show_about(self):
        tkinter.messagebox.showinfo(title="About", message="This program has been done in the PYT subject")

    def create_menu(self):
        self.menubar = tkinter.Menu(self)

        # CREATE THE FILEMENU
        filemenu = tkinter.Menu(self.menubar)
        filemenu.add_command(label="Open Directory", command=self.open_directory)
        filemenu.add_separator()
        filemenu.add_command(label="QUIT", command=self.quit)

        # CREATE THE HELP MENU
        helpmenu = tkinter.Menu(self.menubar)
        helpmenu.add_command(label="Help", command=self.show_help)
        helpmenu.add_command(label="About", command=self.show_about)

        self.menubar.add_cascade(label="File", menu=filemenu)
        self.menubar.add_cascade(label="Help", menu=helpmenu)

        self.master.config(menu=self.menubar)

    def create_options(self):
        frame = tkinter.Frame(self.header)

        #scrollbar = tkinter.Scrollbar(frame, orient=tkinter.VERTICAL)

        label = tkinter.Label(frame, text="Enter UniprotAccession code:")


        #self.seq_listbox = tkinter.Listbox(frame, selectmode=tkinter.SINGLE,
        #                                   height=5,width=100,
        #                                   yscrollcommand=scrollbar.set)

        #scrollbar.config(command=self.seq_listbox.yview)
        #scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)
        label.pack(side=tkinter.LEFT)

        self.seq_listbox.pack(side=tkinter.LEFT, expand=True, fill=tkinter.BOTH)
        entry = tkinter.Entry(root, bd=3, width=6, textvariable=entry_var)
        entry.pack(side=tkinter.LEFT)


        #self.seq_listbox.bind('<<ListboxSelect>>', self.select_seq_from_listbox)

        frame.pack(fill=tkinter.BOTH)


    def create_header(self):
        self.header = tkinter.LabelFrame(self, text="Options", padx=5, pady=5)
        self.create_options()
        self.header.grid(row=0, column=0)



    def createWidgets(self):
        self.create_menu()
        self.create_header()

        #self.create_left_frame()
        #self.create_barplot_options_frame()
        #self.create_sequence_textbox()
        #self.create_graphic_canvas()

        self.grid(row=0)

    def __init__(self, master=None, **kwargs):
        tkinter.Frame.__init__(self, master, **kwargs)

        self.master.wm_title("Program to manage FASTA files with BioPython.")
        self.master.resizable(width=False, height=False)

        self.config(width=600)
        self.config(height=600)

            # DEFINE ATTRIBUTES
        self.seq_list = []
        self.global_frequencies = {}
        self.seq_listbox = None
        self.menubar = None
        self.current_SeqRecord = None
        self.database_freq_checkbox = None
        self.single_seq_freq_checkbox = None
        self.sequence_text = None
        self.regular_expression_var = None
        self.createWidgets()



root = tkinter.Tk()
app = Application(master=root,padx=10, pady=10)

app.mainloop()