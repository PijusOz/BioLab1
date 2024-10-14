import numpy as np
from collections import defaultdict

# Definuojamos start ir stop aminorūgštys
start_amino_acid = 'M'  # Metioninas atitinka start kodoną ATG
stop_amino_acid = '*'   # Stop kodonai yra pažymėti '*'
amino_acids = 'ACDEFGHIKLMNPQRSTVWY*'  # Visos galimos aminorūgštys

# Kodonų lentelė baltymų sekos vertimui
codontab = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'
}

# Funkcija, kuri generuoja atvirkštinę komplementinę seką
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    try:
        return ''.join(complement[base] for base in reversed(seq))
    except KeyError:
        raise ValueError("Invalid nucleotide sequence")

# Funkcija, kuri verčia DNR seką į baltymų seką, naudojant kodonų lentelę
def translate_to_protein(dna_sequence):
    protein_sequence = []
    for i in range(0, len(dna_sequence) - 2, 3):  # Eina per seką trijų nukleotidų grupėmis
        codon = dna_sequence[i:i+3]
        if codon in codontab:
            amino_acid = codontab[codon]
            protein_sequence.append(amino_acid)
        else:
            raise ValueError(f"Invalid codon: {codon}")
    return ''.join(protein_sequence)

# Funkcija, kuri nuskaito FASTA failą ir grąžina seką ir pavadinimą
def read_fasta(filepath):
    try:
        with open(filepath, 'r') as fasta_file:
            lines = fasta_file.readlines()
            header = lines[0].strip()[1:]  # Read the header (without ">")
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>') and line.strip())
        return header, sequence
    except FileNotFoundError:
        raise FileNotFoundError(f"File {filepath} not found.")

# Funkcija, kuri ieško starto ir stopo pozicijų baltymų sekose
def find_protein_sequences(sequence):
    protein_sequence = translate_to_protein(sequence)
    return protein_sequence

# Funkcija, kuri skaičiuoja kodonų (aminorūgščių) dažnius be normalizavimo
def count_codon_frequencies(protein_sequence):
    codon_frequencies = defaultdict(int)

    # Skaičiuojame kodonų dažnius (absolute counts)
    for amino_acid in protein_sequence:
        codon_frequencies[amino_acid] += 1

    return dict(codon_frequencies)

# Funkcija, kuri skaičiuoja dikodonų dažnius be normalizavimo
def count_dicodon_frequencies(protein_sequence):
    dicodon_frequencies = defaultdict(int)

    # Skaičiuojame bendrą dikodonų skaičių (absolute counts)
    total_dicodons = len(protein_sequence) - 1
    for i in range(total_dicodons):
        dicodon = protein_sequence[i:i+2]
        dicodon_frequencies[dicodon] += 1

    return dict(dicodon_frequencies)

# Funkcija, kuri apjungia kodonų ir dikodonų dažnius iš tiesioginės ir atvirkštinės sekos
def combine_frequencies(direct_frequencies, reverse_frequencies):
    combined_frequencies = defaultdict(int)

    # Apjungti kodonų dažnius
    for key in set(direct_frequencies).union(reverse_frequencies):
        combined_frequencies[key] = direct_frequencies.get(key, 0) + reverse_frequencies.get(key, 0)

    return dict(combined_frequencies)

# Funkcija, kuri normalizuoja dažnius pagal bendrą skaičių
def normalize_frequencies(frequencies, total_count):
    return {key: value / total_count for key, value in frequencies.items()}

# Funkcija, kuri skaičiuoja Manhatano atstumą tarp dviejų kodonų/dikodonų dažnių vektorių
def manhattan_distance(freq1, freq2):
    all_keys = set(freq1.keys()).union(set(freq2.keys()))
    distance = 0
    for key in all_keys:
        distance += abs(freq1.get(key, 0) - freq2.get(key, 0))
    return distance

# Funkcija, kuri sukuria atstumų matricą iš kodonų ir dikodonų dažnių naudojant Manhatano atstumą
def create_distance_matrix(combined_frequencies_list):
    num_files = len(combined_frequencies_list)
    distance_matrix = np.zeros((num_files, num_files))
    
    for i in range(num_files):
        for j in range(i, num_files):
            dist = manhattan_distance(combined_frequencies_list[i], combined_frequencies_list[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist  # Matrica yra simetriška
    return distance_matrix

# Funkcija, kuri išsaugo atstumų matricą Phylip formatu
def save_phylip_matrix(distance_matrix, names, output_file):
    num_files = len(names)
    with open(output_file, 'w') as f:
        f.write(f"{num_files}\n")
        for i, name in enumerate(names):
            truncated_name = name[:11]
            distances = ' '.join(f"{dist:.6f}" for dist in distance_matrix[i])
            f.write(f"{truncated_name:<10} {distances}\n")

# Funkcija, kuri apdoroja kelis failus ir išsaugo dvi matricas (viena kodonams, kita dikodonams)
def process_multiple_files(filepaths):
    codon_frequencies_list = []
    dicodon_frequencies_list = []
    names = []

    for filepath in filepaths:
        print(f"\nApdorojamas failas: {filepath}")
        name, sequence = read_fasta(filepath)  # Now returning the name (header)
        names.append(name)  # Save the name from the header

        # Tiesioginė seka
        direct_protein = find_protein_sequences(sequence)
        direct_codon_frequencies = count_codon_frequencies(direct_protein)
        direct_dicodon_frequencies = count_dicodon_frequencies(direct_protein)

        # Normalizuoti dažnius
        total_codons = sum(direct_codon_frequencies.values())
        direct_codon_frequencies = normalize_frequencies(direct_codon_frequencies, total_codons)

        total_dicodons = sum(direct_dicodon_frequencies.values())
        direct_dicodon_frequencies = normalize_frequencies(direct_dicodon_frequencies, total_dicodons)

        # Atvirkštinė komplementinė seka
        reverse_sequence = reverse_complement(sequence)
        reverse_protein = find_protein_sequences(reverse_sequence)
        reverse_codon_frequencies = count_codon_frequencies(reverse_protein)
        reverse_dicodon_frequencies = count_dicodon_frequencies(reverse_protein)

        # Normalizuoti atvirkštinius dažnius
        total_reverse_codons = sum(reverse_codon_frequencies.values())
        reverse_codon_frequencies = normalize_frequencies(reverse_codon_frequencies, total_reverse_codons)

        total_reverse_dicodons = sum(reverse_dicodon_frequencies.values())
        reverse_dicodon_frequencies = normalize_frequencies(reverse_dicodon_frequencies, total_reverse_dicodons)

        # Apjungti tiesioginės ir atvirkštinės sekos dažnius
        combined_codon_frequencies = combine_frequencies(direct_codon_frequencies, reverse_codon_frequencies)
        combined_dicodon_frequencies = combine_frequencies(direct_dicodon_frequencies, reverse_dicodon_frequencies)

        # Pridedame kombinaciją į sąrašus
        codon_frequencies_list.append(combined_codon_frequencies)
        dicodon_frequencies_list.append(combined_dicodon_frequencies)

    # Sukuriame atstumų matricą kodonams
    codon_distance_matrix = create_distance_matrix(codon_frequencies_list)
    save_phylip_matrix(codon_distance_matrix, names, "codon_distance_matrix.phy")

    # Sukuriame atstumų matricą dikodonams
    dicodon_distance_matrix = create_distance_matrix(dicodon_frequencies_list)
    save_phylip_matrix(dicodon_distance_matrix, names, "dicodon_distance_matrix.phy")

# Sąrašas su 8 failais
filepaths = [
    'viruses/data/bacterial1.fasta',
    'viruses/data/bacterial2.fasta',
    'viruses/data/bacterial3.fasta',
    'viruses/data/bacterial4.fasta',
    'viruses/data/mamalian1.fasta',
    'viruses/data/mamalian2.fasta',
    'viruses/data/mamalian3.fasta',
    'viruses/data/mamalian4.fasta'
]

# Pagrindinė funkcija
def main():
    process_multiple_files(filepaths)

# Vykdoma programa
main()