# Импортируем необходимые библиотеки
from CAI import CAI  # Для расчета индекса адаптации кодонов (CAI)
from Bio import SeqIO  # Для работы с биологическими последовательностями
import pandas as pd  # Для работы с табличными данными и вывода результатов


# Функция для валидации последовательности (подготовка к анализу)
def validate_sequence(sequence):
    # Разбиваем последовательность на триплеты, исключая те, что содержат 'N' (неопределенные нуклеотиды)
    sequence = "".join([sequence[i:i + 3] for i in range(0, len(sequence), 3) if 'N' not in sequence[i:i + 3]])

    # Если длина последовательности не кратна 3, обрезаем лишние нуклеотиды в конце
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]

    return sequence


# Основной файл с последовательностью для анализа
fasta_file = "GCA_0007.fasta"

# Читаем FASTA-файл (берем первую запись)
record = next(SeqIO.parse(fasta_file, "fasta"))
# Преобразуем последовательность в строку в верхнем регистре
sequence = str(record.seq).upper()

# Валидируем последовательность
valid_seq = validate_sequence(sequence)

# Читаем референсные последовательности для расчета CAI
reference = [str(seq.seq) for seq in SeqIO.parse("ref.fasta", "fasta")]
results = []  # Список для хранения результатов анализа

# Анализируем генбанковский файл с аннотацией генома
for record in SeqIO.parse("genomic.gb", "genbank"):
    # Перебираем все аннотированные признаки в записи
    for feature in record.features:
        # Пропускаем все признаки, кроме CDS (кодирующих последовательностей)
        if feature.type != "CDS":
            continue

        # Работаем только с аннотированными как "hypothetical protein"
        if feature.qualifiers.get('product') != ['hypothetical protein']:
            # Извлекаем последовательность CDS
            cds_seq = str(feature.location.extract(record.seq).upper())

            # Валидируем последовательность CDS
            valid_seq = validate_sequence(cds_seq)
            # Рассчитываем индекс адаптации кодонов (CAI) относительно референса
            cai = CAI(valid_seq, reference=reference)

            # Собираем аннотации для вывода
            annotations = {
                'Product': feature.qualifiers.get('product', [''])[0],  # Название белка
                'Length': len(valid_seq),  # Длина валидированной последовательности
                'CAI': round(cai, 4)  # Индекс CAI, округленный до 4 знаков
            }
            results.append(annotations)

# Если найдены результаты, выводим их в виде таблицы
if results:
    df = pd.DataFrame(results)  # Создаем DataFrame из результатов
    print("\nCAI Results Summary:")
    # Выводим таблицу без индексов
    print(df.to_string(index=False))