"""
Bioinformatics App for Analyzing UniProtKB Data

This module defines a Tkinter-based application for analyzing UniProtKB
data. It allows users to browse and load UniProtKB files,
perform various analyses, and visualize results through plots.

The application includes functionality to extract UniProtKB entries from
text files, find consensus sequences using regular expressions,
and analyze data based on user-selected options.

Author: Daniel Smit
Date: 18-01-2024
"""
import tkinter as tk
from tkinter import filedialog
from tkinter import scrolledtext, messagebox
import regex as re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class BioAPP:
    """
    Bioinformatics App class for analyzing UniProtKB data.
    """

    def __init__(self, master):
        """
        Initialize the BioAPP.

        Args:
            master (tk.Tk): The root Tkinter window.
        """

        self.master = master
        self.master.title("Bioinformatics App")

        self.label = tk.Label(master, text="Select Input file:")
        self.label.pack(pady=10)

        self.entry = tk.Entry(master, width=40)
        self.entry.pack(pady=10)

        self.browse_button = tk.Button(
            master, text="Browse", command=self.browse_file
        )
        self.browse_button.pack(pady=10)

        self.load_button = tk.Button(
            master, text="Load File", command=self.load_uniprotkb_file
        )
        self.load_button.pack(pady=20)

        self.option_var = tk.StringVar()

        self.radio_a = tk.Radiobutton(
            master,
            text="Per organism: top 10 protein sequences",
            variable=self.option_var,
            value="a",
        )
        self.radio_a.pack()
        self.radio_a["state"] = "disabled"  # Disabled

        self.radio_b = tk.Radiobutton(
            master,
            text="Per organism: found consensus sequences",
            variable=self.option_var,
            value="b",
        )
        self.radio_b.pack()
        self.radio_b["state"] = "disabled"  # Disabled

        self.radio_c = tk.Radiobutton(
            master,
            text="Per accession code: corresponding sequence length",
            variable=self.option_var,
            value="c",
        )
        self.radio_c.pack()
        self.radio_c["state"] = "disabled"  # Disabled

        self.radio_a.select()

        self.analyze_button = tk.Button(
            self.master, text="Analyze Data", command=self.analyze_data
        )

        self.analyze_button.pack()
        self.analyze_button["state"] = "disabled"  # Disabled

        self.text_output = scrolledtext.ScrolledText(
            master, wrap="none", height=20, width=80
        )
        self.text_output.pack(side=tk.LEFT, fill=tk.BOTH)

        self.plot_frame = tk.Frame(master)
        self.plot_frame.pack(pady=20, side=tk.RIGHT, fill=tk.BOTH)

        # Attributes:
        self.filepath = ""
        self.ids = []
        self.ac_numbers = []
        self.dates = []
        self.rec_names = []
        self.gn_names = []
        self.os_names = []
        self.sequence_headers = []
        self.sequences = []
        self.consensus_sequences = ()
        self.entries = []

    def browse_file(self):
        """
        Open a file dialog to browse and select a file.

        Returns:
            self.entry: filled with the filepath
        """
        filepath = filedialog.askopenfilename(
            filetypes=[("Text Files", "*.txt"),
                       ("All Files", "*.*")]
        )
        if filepath:
            self.filepath = filepath
            self.entry.delete(0, tk.END)
            self.entry.insert(0, filepath)

    def process_entries_from_file(self):
        """
        Extracts UniProtKB entries from a file.

        Args:
            self.filepath (str): The path to the UniProtKB text file.

        Returns:
            list: List of dictionaries representing UniProtKB entries.
        """
        entries = []
        sequence_header = ""
        sequence = ""
        in_sequence_block = False

        with open(self.filepath, "r", encoding="UTF-8") as file:
            entry = {}
            for line in file:
                line = line.strip()

                if line.startswith("//"):
                    # Add sequence header and sequence to the entry
                    entry["SQ Header"] = sequence_header.strip()
                    entry["SQ"] = sequence.strip()
                    entries.append(entry)

                    # Reset variables for the next entry
                    entry = {}
                    sequence_header = ""
                    sequence = ""
                    in_sequence_block = False
                elif line.startswith("SQ   "):
                    # Start of the sequence block,
                    # set in_sequence_block flag
                    in_sequence_block = True
                    sequence_header = line[5:].strip()
                elif in_sequence_block:
                    # If in the sequence block,
                    # append line to the sequence
                    sequence += line.strip().replace(" ", "")
                else:
                    # If not in the sequence block,
                    # process the line as a regular entry field
                    key, value = line[:5].strip(), line[5:].strip()
                    if key not in entry:
                        entry[key] = []
                    entry[key].append(value)

            self.entries = entries
            return entries

    def extract_entries(self):
        """
        Extract Entries from filedata.

        """
        entries = self.process_entries_from_file()

        # Separate information into lists
        for entry in entries:
            self.ids.append(entry.get("ID", [""])[0])
            self.ac_numbers.append(entry.get("AC", [""])[0])
            self.dates.append(entry.get("DT", []))
            self.rec_names.append(entry.get("DE", [""])[0])
            self.gn_names.append(entry.get("GN", [""])[0])
            self.os_names.append(entry.get("OS", [""])[0])
            self.sequence_headers.append(entry.get("SQ Header", ""))
            self.sequences.append(entry.get("SQ", ""))

    def load_uniprotkb_file(self):
        """
        Load File, and enable buttons when file loaded successful.

        """
        try:
            self.extract_entries()
            self.consensus_sequences = self.find_consensus_sequences()
            if self.filepath:
                print(f"File Path: {self.filepath}")

                # The analyze_button + the select radio buttons
                # will now be enabled
                self.analyze_button["state"] = "normal"
                self.radio_a["state"] = "normal"
                self.radio_b["state"] = "normal"
                self.radio_c["state"] = "normal"

            else:
                # Keep the buttons disabled
                self.analyze_button["state"] = "disabled"
                self.radio_a["state"] = "disabled"
                self.radio_b["state"] = "disabled"
                self.radio_c["state"] = "disabled"

        # Exceptions
        except FileNotFoundError:
            messagebox.showerror("Error",
                                 f"File '{self.filepath}' "
                                 f"not found!")
            print(f"Error: File '{self.filepath}' not Found!.")
        except IsADirectoryError:
            messagebox.showerror(
                "Error", f"'{self.filepath}' "
                         f"is a directory, not a file."
            )
            print(f"Error: '{self.filepath}' is a directory, "
                  f"not a file.")
        except PermissionError:
            messagebox.showerror("Error",
                                 f"No permission "
                                 f"to access '{self.filepath}'.")
            print(f"Error: No permission to access '{self.filepath}'.")
        # except Exception as e:
        #     messagebox.showerror("Error", f"Error while "
        #                                   f"reading file: {e}")
        #     print(f"Error while reading file: {e}")

    def find_consensus_sequences(self):
        """
        Find consensus sequences in UniProtKB entries using regex.

        Returns:
            list: List of tuples containing organism
            and matching sequences.
        """
        try:
            # Searching for consensus sequences with regex
            prosite_pattern = (r"[GSTNP][A-Z]{6}[FYVHR]"
                               r"[IVN][KEP][A-Z]G[STIVKRQ]"
                               r"Y[DNQKRMV][EP][A-Z]{3}"
                               r"[LIMVA]")

            matches = []
            for i in self.entries:
                match_sequences = set(re.findall(prosite_pattern,
                                                 i["SQ"]))
                if match_sequences:
                    matches.append((i["OS"], match_sequences))

            return matches

        except re.error as e:
            return print(f"Regex error: {e}")

    # FINISHED
    def analyze_data(self):
        """
        Analyze UniProtKB data based on selected analysis option.

        Returns:
            None
        """

        option = self.option_var.get()

        # Per organism: top 10 Protein Sequences
        if option == "a":
            self.proteins_in_consensus_per_organism()

        # Per organism: found consensus sequences
        elif option == "b":
            self.consensus_per_organism()

        # Per accession code: corresponding sequence length
        elif option == "c":
            self.sequence_length_by_accession()

    def proteins_in_consensus_per_organism(self):
        """
        Analyze option 'a': Per organism, 
        count proteins with the corresponding consensus sequence.

        Returns:
            None
        """

        try:
            # a. Per organism, count how many proteins
            # the corresponding consensus sequence has.
            organism_protein_count = {}

            for entry in self.entries:
                organism = entry.get("OS", [""])[0]

                # a. Count proteins per organism
                if organism:
                    organism_protein_count[organism] = (
                        organism_protein_count.get(organism, 0) + 1
                    )

            top_10_organisms = dict(
                sorted(
                    organism_protein_count.items(),
                    key=lambda item: item[1],
                    reverse=True,
                )[:10]
            )

            self.show_plot(
                textual_data={'title': "Top 10 Protein Sequences "
                                       "per organism",
                              'x_label': "Organism",
                              'y_label': "Amount of Proteins"},
                data=top_10_organisms,
                row=1,
                col=1,
            )

        # Handling Exceptions
        except IndexError as e:
            print(f"Error while going through consensus sequences: {e}")

    def consensus_per_organism(self):
        """
        Analyze option 'b': Per organism,
        list the consensus sequences found.

        Returns:
            None
        """

        # b. Per organism, list the consensus sequences found.
        organism_consensus_sequences = {}

        #  b. List consensus sequences per organism
        for organism, consensus_sequences \
                in self.consensus_sequences:
            if organism[0] not in organism_consensus_sequences:
                organism_consensus_sequences[organism[0]] = []
            (organism_consensus_sequences[organism[0]]
             .extend(consensus_sequences))

        # Show consensus sequences per organism in de Text-widget
        self.text_output.delete(1.0,
                                tk.END)  # Remove previous contents

        # show all the organisms and found sequences
        # in the Textfield on Main App.
        for organism, sequences \
                in organism_consensus_sequences.items():
            # Makes sure for no double sequences in 'sequences' list
            sequences = set(sequences)
            output_text = (f"{organism}: {len(sequences)}"
                           f" consensus sequences\n")

            for i, sequence in enumerate(sequences, start=1):
                if i <= 20:
                    output_text += (f"  Consensus Sequence "
                                    f"{i}: {sequence}\n")
                else:
                    break

            self.text_output.insert(tk.END, output_text)

        # b. Print consensus sequences per organism
        for organism, sequence_count \
                in organism_consensus_sequences.items():
            print(f"{organism}: {len(sequence_count)}"
                  f" consensus sequences")

        # if you would like to see a plot, uncomment the next part:
        # self.show_plot(
        #     "total found consensus sequences per organism",
        #     "organism",
        #     "total found sequences",
        #     data_for_plot,
        #     1,
        #     2)

    def sequence_length_by_accession(self):
        """
        Analyze option 'c': Per accession code, 
        show the corresponding sequence length.

        Returns:
            None
        """

        accession_sequence_length = {}

        for entry in self.entries:
            accession = entry.get("AC", [""])[0]
            sequence = entry.get("SQ", [""])

            # Show sequence length per accession code
            if accession and sequence:
                sequence_length = len(sequence)
                accession_sequence_length[accession] = sequence_length

        # Print sequence length per accession code
        print("\n Sequence length per accession code:")
        for accession, length in accession_sequence_length.items():
            print(f"{accession}: {length} amino acids")

        self.update_text_widget(accession_sequence_length)

        # If you would like to see a plot, uncomment the next part:
        # self.show_plot(
        #     "Sequence length per accession code",
        #     "Accession code",
        #     "Sequence length",
        #     accession_sequence_length,
        #     1,
        #     3,
        # )

    # Update Text widget with sequence length per accession code
    def update_text_widget(self, data):
        """
        Update the text widget with sequence length per accession code.

        Args:
            data (dict): Dictionary containing accession codes
            and sequence lengths.

        Returns:
            None
        """
        self.text_output.delete(1.0, tk.END)  # Clear text

        # Print sequence length per accession code
        self.text_output.insert(tk.END,
                                "Sequence length per accession "
                                "code:\n")
        for accession, length in data.items():
            self.text_output.insert(tk.END, f"{accession}: "
                                            f"{length} amino acids\n")

    # Show Plot in Main App (Bio App)
    def show_plot(self, textual_data, data, row, col):
        """
        Display a plot in the main app window.

        Args:
            textual_data: contains the title, x_axis and y_axis of the
            plot.
            data (dict): Data for the plot.
            row (int): Row position in the Tkinter grid.
            col (int): Column position in the Tkinter grid.

        Returns:
            None
        """

        # Create the plot
        fig, ax = plt.subplots()
        fig.set_figheight(7)
        ax.bar(data.keys(), data.values())
        ax.set_xlabel(textual_data.get("x_label"))
        ax.set_ylabel(textual_data.get("y_label"))
        ax.set_title(textual_data.get("title"))
        for tick in ax.get_xticklabels():
            tick.set_fontsize(8)
        plt.xticks(rotation=25, ha="right")

        # Show plot on the GUI screen
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.get_tk_widget().grid(row=row, column=col)
        # fig.tight_layout()
        # canvas.get_tk_widget().pack()
        canvas.draw()


# Initialization
if __name__ == "__main__":
    root = tk.Tk()
    app = BioAPP(root)
    root.mainloop()
