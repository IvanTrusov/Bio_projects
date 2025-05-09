# Импортируем необходимые модули
from Bio import SeqIO  # Для работы с последовательностями в формате FASTA
from collections import defaultdict  # Для создания словарей со значениями по умолчанию


# Функция для валидации последовательности (удаление участков с N и выравнивание по триплетам)
def validate_sequence(sequence):
    # Разбиваем последовательность на триплеты (кодоны), исключая те, что содержат 'N'
    sequence = "".join([sequence[i:i + 3] for i in range(0, len(sequence), 3) if 'N' not in sequence[i:i + 3]])

    # Если длина последовательности не кратна 3, обрезаем конец
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]

    return sequence


# Функция для расчета частоты встречаемости кодонов
def calculate_codon_frequency(sequence):
    # Разбиваем последовательность на кодоны (триплеты)
    codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]

    # Создаем словарь для подсчета количества каждого кодона
    codon_counts = defaultdict(int)
    for codon in codons:
        codon_counts[codon] += 1

    # Общее количество кодонов в последовательности
    total = len(codons)

    # Рассчитываем частоту встречаемости каждого кодона
    frequencies = {cod: count / total for cod, count in codon_counts.items()}

    # Выводим топ-5 самых частых кодонов
    print(f'top_5_codons:')
    for codon, frequency in sorted(frequencies.items(), key=lambda x: -x[1])[:5]:
        print(f"{codon}: {frequency:.6f}")

    # Выводим количество уникальных кодонов
    print(f'unique_codons: {len(codon_counts)}')

    return frequencies


# Функция для расчета частоты встречаемости гексамеров (6-нуклеотидных последовательностей)
def calculate_hexamer_frequency(sequence):
    # Разбиваем последовательность на гексамеры
    hexamers = [sequence[i:i + 6] for i in range(len(sequence) - 5)]

    # Создаем словарь для подсчета количества каждого гексамера
    hex_counts = defaultdict(int)
    for hex in hexamers:
        hex_counts[hex] += 1

    # Общее количество гексамеров в последовательности
    total = len(hexamers)

    # Рассчитываем частоту встречаемости каждого гексамера
    frequencies = {hex: count / total for hex, count in hex_counts.items()}

    # Выводим топ-5 самых частых гексамеров
    print(f'top_5_hexamers:')
    for hexamer, frequency in sorted(frequencies.items(), key=lambda x: -x[1])[:5]:
        print(f"{hexamer}: {frequency:.6f}")

    # Выводим количество уникальных гексамеров
    print(f'unique_hexamers: {len(hex_counts)}')

    return frequencies


fasta_file = ""  # Имя файла с последовательностью

# Читаем FASTA-файл (берем первую запись, если их несколько)
record = next(SeqIO.parse(fasta_file, "fasta"))
# Преобразуем последовательность в строку в верхнем регистре
sequence = str(record.seq).upper()

# Валидируем последовательность (удаляем триплеты с N и выравниваем длину)
valid_seq = validate_sequence(sequence)

# Рассчитываем и выводим частоты кодонов
codon_frequency = calculate_codon_frequency(valid_seq)
# Рассчитываем и выводим частоты гексамеров
hexamer_frequencies = calculate_hexamer_frequency(valid_seq)
