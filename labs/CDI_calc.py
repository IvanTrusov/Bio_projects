from collections import defaultdict
import math

def validate_sequence(sequence):
    # Разбиваем последовательность на триплеты (кодоны), исключая те, что содержат 'N'
    sequence = "".join([sequence[i:i + 3] for i in range(0, len(sequence), 3) if 'N' not in sequence[i:i + 3]])

    # Если длина последовательности не кратна 3, обрезаем конец
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]

    return sequence
def read_fasta(file_path):
    """Читает FASTA-файл и возвращает словарь {header: sequence}."""
    sequences = {}
    with open(file_path, 'r') as file:
        header = ""
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header and sequence:
                    sequences[header] = sequence.upper()
                    sequence = ""
                header = line[1:]  # Убираем '>'
            else:
                sequence += line
        if header and sequence:
            sequences[header] = sequence.upper()
    return sequence.upper()


def calculate_cdi(sequences):
    """Вычисляет Codon Diversity Index (CDI) для всех последовательностей."""
    codon_counts = defaultdict(int)
    total_codons = 0


    # Разбиваем на кодоны и считаем их
    for i in range(0, len(sequences), 3):
        codon = sequences[i:i + 3]
        if len(codon) == 3:
            codon_counts[codon] += 1
            total_codons += 1


    # Рассчитываем CDI (аналог индекса Симпсона)
    cdi = 0.0
    for count in codon_counts.values():
        p = count / total_codons
        cdi += p * p

    cdi = 1 - cdi  # Чем ближе к 1, тем выше разнообразие
    return cdi


# Пример использования
if __name__ == "__main__":
    fasta_file = ""  # Замените на свой файл
    sequences = read_fasta(fasta_file)
    valid_seq = validate_sequence(sequences)

    if not valid_seq:
        print("Файл пуст или не содержит последовательностей.")
    else:
        cdi = calculate_cdi(valid_seq)
        print(f"Индекс разнообразия кодонов (CDI) = {cdi:.4f}")