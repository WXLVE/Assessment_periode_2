import tkinter as tk
import regex as re
import matplotlib.pyplot as plt
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


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
        
        with open(self.filepath, 'r') as file:
            entry = {}
            for line in file:
                line = line.strip()
                
                if line.startswith('//'):
                    # Add sequence header and sequence to the entry
                    entry['SQ Header'] = sequence_header.strip()
                    entry['SQ'] = sequence.strip()
                    entries.append(entry)
                    
                    # Reset variables for the next entry
                    entry = {}
                    sequence_header = ""
                    sequence = ""
                    in_sequence_block = False
                elif line.startswith('SQ   '):
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
            self.ids.append(entry.get('ID', [''])[0])
            self.ac_numbers.append(entry.get('AC', [''])[0])
            self.dates.append(entry.get('DT', []))
            self.rec_names.append(entry.get('DE', [''])[0])
            self.gn_names.append(entry.get('GN', [''])[0])
            self.os_names.append(entry.get('OS', [''])[0])
            self.sequence_headers.append(entry.get('SQ Header', ''))
            self.sequences.append(entry.get('SQ', ''))

# FINISHED: Load uniprotKB files and parse data
    def load_uniprotkb_file(self):
        try:
            self.extract_entries()
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

    def find_consensus_sequences(self, sequences):
        # Zoeken naar consensus sequenties met regex
        consensus_pattern = r"[GSTNP]-x(6)-[FYVHR]-[IVN]-[KEP]-x-G-[STIVKRQ]-Y-[DNQKRMV]-[EP]-x(3)-[LIMVA]"
        consensus_sequences = []

        for sequence in sequences:
            matches = re.finditer(consensus_pattern, sequence)
            for match in matches:
                consensus_sequences.append(match.group())

        print(consensus_sequences)
        return consensus_sequences

# NEEDS WORK STILL
    def analyze_data(self):
        # Implementeer de gewenste analyse-opties (a, b, c)
        # Hier nemen we aan dat je de sequenties al hebt geladen en geparseerd.

        option = self.option_var.get()

        sample_sequences = [
            "G-x(6)-F-IV-KE-G-S-Y-D-E-x(3)-L",
            "N-x(6)-H-IVN-KEP-G-S-T-Y-Q-M-x(3)-I",
            "T-x(6)-F-IVN-KE-G-S-T-Y-D-E-x(3)-LIMVA",
        ]

        consensus_sequences = self.find_consensus_sequences(sample_sequences)

        if option == "a":
            # Laat zien per organisme hoeveel eiwitten de bijbehorende consensus sequentie hebben.
            organism_counts = {}
            for sequence in sample_sequences:
                organism = sequence.split("-")[3]  # Bijvoorbeeld: "IVN"
                if organism in organism_counts:
                    organism_counts[organism] += 1
                else:
                    organism_counts[organism] = 1

            # print("Per organisme: aantal eiwitten met consensus sequentie")
            # for organism, count in organism_counts.items():
            #     print(f"{organism}: {count} eiwitten")

            # # Weergeef de resultaten in een staafdiagram
            # plt.bar(organism_counts.keys(), organism_counts.values())
            # plt.xlabel("Organisme")
            # plt.ylabel("Aantal eiwitten")
            # plt.title("Aantal eiwitten per organisme met consensus sequentie")
            # plt.show()

            self.show_plot(
                "Aantal eiwitten per organisme met consensus sequentie",
                "Organisme",
                "Aantal eiwitten",
                organism_counts,
            )

        elif option == "b":
            # Laat zien per organisme welke consensus sequenties er zijn gevonden.
            organism_sequences = {}
            for sequence in sample_sequences:
                organism = sequence.split("-")[3]  # Bijvoorbeeld: "IVN"
                if organism in organism_sequences:
                    organism_sequences[organism].append(sequence)
                else:
                    organism_sequences[organism] = [sequence]

            print("Per organisme: consensus sequenties")
            for organism, sequences in organism_sequences.items():
                print(f"{organism}: {sequences}")

            # Weergeef de resultaten in een staafdiagram
            # plt.bar(
            #     organism_sequences.keys(),
            #     [len(sequences) for sequences in organism_sequences.values()],
            # )
            # plt.xlabel("Organisme")
            # plt.ylabel("Aantal gevonden consensus sequenties")
            # plt.title("Aantal gevonden consensus sequenties per organisme")
            # plt.show()

            # self.show_plot("Aantal gevonden consensus sequenties per organisme", "Organisme", "Aantal gevonden sequenties", {org: len(sequences) for org, sequences in organism_sequences.items()})

        elif option == "c":
            # Toon per accessiecode de bijbehorende sequentielengte
            accessie_lengths = {}
            for sequence in sample_sequences:
                accessie = sequence.split("-")[1]  # Bijvoorbeeld: "x(6)"
                length = len(sequence)
                accessie_lengths[accessie] = length

            # print("Per accessiecode: sequentielengte")
            # for accessie, length in accessie_lengths.items():
            #     print(f"{accessie}: {length}")

            # # Weergeef de resultaten in een staafdiagram
            # plt.bar(accessie_lengths.keys(), accessie_lengths.values())
            # plt.xlabel("Accessiecode")
            # plt.ylabel("Sequentielengte")
            # plt.title("Sequentielengte per accessiecode")
            # plt.show()

            self.show_plot(
                "Sequentielengte per accessiecode",
                "Accessiecode",
                "Sequentielengte",
                accessie_lengths,
            )
# MAIN DEF FOR SHOWING PLOTS IN MAIN.APP
    def show_plot(self, title, x_label, y_label, data):
        # Maak de plot
        fig, ax = plt.subplots()
        ax.bar(data.keys(), data.values())
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)

        # Toon de plot in het GUI-venster
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.get_tk_widget().pack()
        canvas.draw()


# Initialisatie
if __name__ == "__main__":
    root = tk.Tk()
    app = BioAPP(root)
    root.mainloop()
