# Импортируем необходимый модуль из библиотеки BioPython для работы с последовательностями
from Bio import SeqIO

# Функция для расчета GC-состава (процентное содержание гуанина и цитозина)
def gc_content(seq):
    # Считаем количество нуклеотидов G и C в последовательности
    gc = seq.count('G') + seq.count('C')
    # Возвращаем процентное содержание GC, округленное до 4 знаков после запятой
    return round(((gc / len(seq)) * 100), 4)

# Функция для поиска изохор - участков ДНК с относительно однородным составом
def find_isochores(gc_values, positions, threshold=1.0):
    isochores = []  # Список для хранения найденных изохор
    current_start = positions[0]  # Начальная позиция текущей изохоры
    current_gc = gc_values[0]     # GC-состав текущей изохоры

    # Проходим по всем значениям GC и их позициям
    for pos, gc in zip(positions[1:], gc_values[1:]):
        # Если разница в GC-составе превышает порог, значит началась новая изохора
        if abs(gc - current_gc) > threshold:
            # Добавляем найденную изохору (длина, средний GC-состав)
            isochores.append((pos-current_start, current_gc))
            # Обновляем начальную позицию и GC-состав для новой изохоры
            current_start = pos
            current_gc = gc
    # Добавляем последнюю изохору
    isochores.append((positions[-1] - current_start, current_gc))

    # Фильтруем изохоры, оставляя только те, что длиннее 1000 нуклеотидов
    return [isochor for isochor in filter(lambda x: x[0] > 1000, isochores)]

# Параметры скользящего окна
window_size = 1000  # Размер окна для анализа
step_size = 100     # Шаг перемещения окна

# Списки для хранения значений GC-состава и их позиций
gc_values = []
positions = []

# Читаем FASTA-файл с последовательностью
record = SeqIO.read("", "fasta") #your fasta
# Преобразуем последовательность в строку в верхнем регистре
sequence = str(record.seq).upper()

# Проходим по последовательности скользящим окном
for step in range(0, len(sequence) - window_size + 1, step_size):
    # Выделяем текущее окно последовательности
    window = sequence[step:step+window_size]
    # Рассчитываем GC-состав для окна и добавляем в список
    gc_values.append(gc_content(window))
    # Добавляем текущую позицию в список
    positions.append(step)

# Находим изохоры с порогом изменения GC-состава 0.5%
isochores = find_isochores(gc_values, positions, threshold=0.5)
# Выводим количество найденных изохор
print(len(isochores))
# Выводим список изохор (каждая изохора представлена кортежем (длина, GC-состав))
print(isochores)
