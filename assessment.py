import tkinter as tk
import regex as re
import matplotlib.pyplot as plt
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import scrolledtext

class BioAPP:
    def __init__(self, master):
        self.master = master
        self.master.title("Bioinformatics App")

        self.label = tk.Label(master, text="Selecteer het invoerbestand:")
        self.label.pack(pady=10)

        self.entry = tk.Entry(master, width=40, text="Input Filepath goes here")
        self.entry.pack(pady=10)

        self.browse_button = tk.Button(
            master, text="Bladeren", command=self.browse_file
        )
        self.browse_button.pack(pady=10)

        self.load_button = tk.Button(
            master, text="Inladen", command=self.load_uniprotkb_file
        )
        self.load_button.pack(pady=20)

        self.option_var = tk.StringVar()

        self.radio_a = tk.Radiobutton(
            master,
            text="Per organisme: aantal eiwitten met consensus sequentie",
            variable=self.option_var,
            value="a",
        )
        self.radio_a.pack()
        self.radio_a["state"] = "disabled"  # Disabled

        self.radio_b = tk.Radiobutton(
            master,
            text="Per organisme: welke consensus sequenties er zijn gevonden",
            variable=self.option_var,
            value="b",
        )
        self.radio_b.pack()
        self.radio_b["state"] = "disabled"  # Disabled

        self.radio_c = tk.Radiobutton(
            master,
            text="Per accessiecode: bijbehorende sequentielengte",
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

        self.text_output = scrolledtext.ScrolledText(master, wrap="none", height=20, width=80)
        self.text_output.pack(pady=20)


        self.plot_frame = tk.Frame(master)
        self.plot_frame.pack()

    def browse_file(self):
        filepath = filedialog.askopenfilename(
            filetypes=[("Tekstbestanden", "*.txt"), ("Alle bestanden", "*.*")]
        )
        if filepath:
            self.entry.delete(0, tk.END)
            self.entry.insert(0, filepath)
            self.filepath = filepath

    def extract_entries(self):
        """
        Extracts UniProtKB entries from a file.

        Args:
            file_path (str): The path to the UniProtKB text file.

        Returns:
            list: List of dictionaries representing UniProtKB entries.
        """
        self.ids = []
        self.ac_numbers = []
        self.dates = []
        self.rec_names = []
        self.gn_names = []
        self.os_names = []
        self.sequence_headers = []
        self.sequences = []

        entries = []
        sequence_header = ""
        sequence = ""
        in_sequence_block = False

        with open(self.filepath, "r") as file:
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
                    # Start of the sequence block, set in_sequence_block flag
                    in_sequence_block = True
                    sequence_header = line[5:].strip()
                elif in_sequence_block:
                    # If in the sequence block, append line to the sequence
                    sequence += line.strip().replace(" ", "")
                else:
                    # If not in the sequence block, process the line as a regular entry field
                    key, value = line[:5].strip(), line[5:].strip()
                    if key not in entry:
                        entry[key] = []
                    entry[key].append(value)

        # return entries

        # # Read UniProtKB entries from the file
        # file_path = 'uniprotkb_xref_prosite_PS00800_2023_12_18.txt'
        # entries = extract_entries(file_path)

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

        # for i in self.sequences[:10]:
        #     print(i)

    # FINISHED: Load uniprotKB files and parse data
    def load_uniprotkb_file(self):
        try:
            self.extract_entries()
            self.entries = self.extract_data()
            self.consensus_sequences = self.find_consensus_sequences()
            if self.filepath:
                print(f"File Path: {self.filepath}")
                # print(accessions, sequences, IDS)
                # self.analyze_data(accessions, sequences)

                # The analyze_button + the select radio buttons will now be enabled
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

        except FileNotFoundError:
            print(f"Error: File '{self.filepath}' not Found!.")
        except Exception as e:
            print(f"Error while reading file: {e}")

    def find_consensus_sequences(self):
        # Zoeken naar consensus sequenties met regex
        prosite_pattern = r"[GSTNP][A-Z]{6}[FYVHR][IVN][KEP][A-Z]G[STIVKRQ]Y[DNQKRMV][EP][A-Z]{3}[LIMVA]"
        regex_pattern = r"[GSTNP]x{6}[FYVHR][IVN][KEP]xG[STIVKRQ]Y[DNQKRMV][EP]x(3)[LIMVA]"
        # Convert PROSITE pattern to regex pattern
        # regex_pattern =

        # print(regex_pattern)

        matches = []
        # for seq in self.sequences:
        #     matchs = re.findall(prosite_pattern, seq)
        #     for m in matchs:
        #         matches.append(m)
                
        for i in self.entries:
            match_sequences = [match.group() for match in re.finditer(prosite_pattern, i["SQ"])]
            if match_sequences:
                matches.append((i["OS"], match_sequences))

        return matches
        


    def extract_data(self):
        """
        Extracts UniProtKB entries from a file.

        Args:
            file_path (str): The path to the UniProtKB text file.

        Returns:
            list: List of dictionaries representing UniProtKB entries.
        """

        entries = []
        sequence_header = ""
        sequence = ""
        in_sequence_block = False

        with open(self.filepath, "r") as file:
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
                    # Start of the sequence block, set in_sequence_block flag
                    in_sequence_block = True
                    sequence_header = line[5:].strip()
                elif in_sequence_block:
                    # If in the sequence block, append line to the sequence
                    sequence += line.strip().replace(" ", "")
                else:
                    # If not in the sequence block, process the line as a regular entry field
                    key, value = line[:5].strip(), line[5:].strip()
                    if key not in entry:
                        entry[key] = []
                    entry[key].append(value)

        return entries

    # NEEDS WORK STILL
    def analyze_data(self):
        # Implementeer de gewenste analyse-opties (a, b, c)
        # Hier nemen we aan dat je de sequenties al hebt geladen en geparseerd.

        option = self.option_var.get()


        

        if option == "a":
            # Laat zien per organisme hoeveel eiwitten de bijbehorende consensus sequentie hebben.
            try:

                # a. Per organisme, tell how many proteins have the corresponding consensus sequence.
                organism_protein_count = {}

                for entry in self.entries:
                    organism = entry.get("OS", [""])[0]

                    # a. Count proteins per organism
                    if organism:
                        organism_protein_count[organism] = (
                            organism_protein_count.get(organism, 0) + 1
                        )

                # a. Print protein count per organism
                # print("a. Protein count per organism:")
                # for organism, count in organism_protein_count.items():
                    # print(f"{organism}: {count} proteins")
                    
                    
                top_10_organisms = dict(sorted(organism_protein_count.items(), key=lambda item: item[1], reverse=True)[:10])



                self.show_plot(
                    "Aantal eiwitten per organisme met consensus sequentie",
                    "Organisme",
                    "Aantal eiwitten",
                    top_10_organisms,
                    1,
                    1,
                )
            # Handlling Exceptions
            except Exception as e:
                print(f"Error while going through consensus sequences: {e}")

        elif option == "b":
            # b. Per organisme, list the consensus sequences found.
            organism_consensus_sequences = {}

            for entry in self.entries:
                organism = entry.get("OS", [""])[0]
                accession = entry.get("AC", [""])[0]
                sequence = entry.get("SQ", [""])

                # b. List consensus sequences per organism
                for organism, consensus_sequences in self.consensus_sequences:
                    if organism[0] not in organism_consensus_sequences:
                        organism_consensus_sequences[organism[0]] = []
                    organism_consensus_sequences[organism[0]].extend(consensus_sequences)

                # Toon consensus sequenties per organisme in de Text-widget
                self.text_output.delete(1.0, tk.END)  # Verwijder eerdere inhoud
                # for organism, sequences in organism_consensus_sequences.items():
                #                 self.text_output.insert(tk.END, f"{organism}: {len(sequences)} consensus sequences\n")
                #                 for i, sequence in enumerate(sequences, start=1):
                #                     self.text_output.insert(tk.END, f"  Consensus Sequence {i}: {sequence}\n")
                                    
                                    
                for organism, sequences in organism_consensus_sequences.items():
                    self.text_output.insert(tk.END, f"{organism}: {len(sequences)} consensus sequences\n")
                    
                    for i, sequence in enumerate(sequences, start=1):
                        if i <= 20:
                            self.text_output.insert(tk.END, f"  Consensus Sequence {i}: {sequence}\n")
                        else:
                            break

            # b. Print consensus sequences per organism
            # for organism, sequence_count in organism_consensus_sequences.items():
            #     print(f"{organism}: {sequence_count} consensus sequences")

            # Convert the dictionary to a list of tuples
            data_for_plot = dict(organism_consensus_sequences.items())

            # self.show_plot(
            #     "Aantal gevonden consensus sequenties per organisme",
            #     "Organisme",
            #     "Aantal gevonden sequenties",
            #     data_for_plot,
            #     1,
            #     2,
            # )

        elif option == "c":
            # Toon per accessiecode de bijbehorende sequentielengte
            
            # c. Per accessiecode, show the corresponding sequence length.
            accession_sequence_length = {}

            for entry in self.entries:
                organism = entry.get("OS", [""])[0]
                accession = entry.get("AC", [""])[0]
                sequence = entry.get("SQ", [""])

            # c. Show sequence length per accession code
            if accession and sequence:
                sequence_length = len(sequence)
                accession_sequence_length[accession] = sequence_length

            # c. Print sequence length per accession code
            # print("\nc. Sequence length per accession code:")
            # for accession, length in accession_sequence_length.items():
            #     print(f"{accession}: {length} amino acids")

            self.show_plot(
                "Sequentielengte per accessiecode",
                "Accessiecode",
                "Sequentielengte",
                accession_sequence_length,
                1,
                3,
            )

    # MAIN DEF FOR SHOWING PLOTS IN MAIN.APP
    def show_plot(self, title, x_label, y_label, data, row, col):
        # Maak de plot
        fig, ax = plt.subplots()
        ax.bar(data.keys(), data.values())
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        for tick in ax.get_xticklabels():
            tick.set_fontsize(8)
        plt.xticks(rotation=45, ha='right')


        # Toon de plot in het GUI-venster
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.get_tk_widget().grid(row=row, column=col)
        # fig.tight_layout()
        # Voeg de navigatietoolbar toe voor zoomen en verschuiven
        toolbar = NavigationToolbar2Tk(fig, self.plot_frame)
        toolbar.update()
        # canvas.get_tk_widget().pack()
        canvas.draw()


# Initialisatie
if __name__ == "__main__":
    root = tk.Tk()
    app = BioAPP(root)
    root.mainloop()
