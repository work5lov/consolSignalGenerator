import numpy as np
import struct
import random
import datetime
from multiprocessing import Pool, cpu_count
from time import time
from tqdm import tqdm  # Импортируем tqdm для отображения прогресса

def signal_Gen(fs, period, Amp, step):
    out = {}
    end_time_sec = (period * step) / 2
    start_time_sec = 0
    f = 1 / (period * step)
    num_samples = round((end_time_sec - start_time_sec) * fs)
    time = [start_time_sec + i / fs for i in range(num_samples)]
    signal = [Amp * np.sin(2 * np.pi * f * t) for t in time]
    per = (period * step) / (step * 2 + 1)
    indices = [i for i in range(0, num_samples, round(per * fs))]
    signal_selected = [signal[i] for i in indices]
    c = step + 1
    y_vect = signal_selected[:c]
    x_vect = [i * ((period * step) / step) * 100 for i in range(c)]
    for i in range(len(x_vect)):
        out[x_vect[i]] = y_vect[i]
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(timestamp)
    print("Отсчеты сигнала:", out)
    return out

def generate_signal_chunk(args):
    """
    Генерирует часть сигнала для заданного диапазона индексов.
    """
    start_idx, end_idx, t_vect, a_vect, pulse_duration_in_samples, fs = args
    chunk_signal = np.zeros(end_idx - start_idx)

    indices1 = (t_vect * fs).astype(int)    
    
    for i in range(start_idx, end_idx):
        idx_local = i - start_idx  # Локальный индекс внутри части
        for j in range(len(indices1)):
            # indices1 = int(t_vect[j] * fs)            
            # Rising edge
            if indices1[j] <= i < int(indices1[j] + 0.2 * pulse_duration_in_samples):
                x1 = indices1[j]
                y1 = 0
                x2 = int(indices1[j] + 0.2 * pulse_duration_in_samples)
                y2 = a_vect[j]
                slope = (y2 - y1) / (x2 - x1)
                intercept = y1 - slope * x1
                chunk_signal[idx_local] = slope * i + intercept + random.uniform(0.01 * a_vect[j], 0.02 * a_vect[j])
            
            # Falling edge
            elif int(indices1[j] + 0.8 * pulse_duration_in_samples) <= i < int(indices1[j] + pulse_duration_in_samples):
                x1 = indices1[j] + 0.8 * pulse_duration_in_samples
                y1 = a_vect[j]
                x2 = int(indices1[j] + pulse_duration_in_samples)
                y2 = 0
                slope = (y2 - y1) / (x2 - x1)
                intercept = y1 - slope * x1
                chunk_signal[idx_local] = slope * i + intercept + random.uniform(0.01 * a_vect[j], 0.02 * a_vect[j])
            
            # Steady state
            elif int(indices1[j] + 0.2 * pulse_duration_in_samples) <= i < int(indices1[j] + 0.8 * pulse_duration_in_samples):
                chunk_signal[idx_local] = a_vect[j] + random.uniform(0.01 * a_vect[j], 0.02 * a_vect[j])
            
            # Noise outside pulses
            else:
                chunk_signal[idx_local] = random.uniform(-0.01, 0.01)
    
    return chunk_signal


def generate_signal(signal, pulse_duration, fs):
    """
    Генерация сигнала с использованием нескольких процессов.
    """
    t_vect = np.array([k for k in signal.keys()])
    a_vect = np.array([v for v in signal.values()])    
    
    # Calculate the length of the signal array
    num_samples = int((t_vect[-1] + pulse_duration) * fs)

    # Calculate the duration of the pulse in samples
    pulse_duration_in_samples = int(pulse_duration * fs)
    # print("Длительность импульса в отсчетах:", pulse_duration_in_samples)

    indices1 = (t_vect * fs).astype(int)
    print("Индексы:",indices1)

    print("Генерация сигнала...")

    # Determine the number of CPU cores
    num_cores = cpu_count()
    print(f"Используется {num_cores} ядер процессора.")

    # First level: split data into 10 parts
    first_level_parts = 10
    first_level_chunk_size = num_samples // first_level_parts
    first_level_chunks = []
    for part_idx in range(first_level_parts):
        start_idx = part_idx * first_level_chunk_size
        end_idx = (part_idx + 1) * first_level_chunk_size if part_idx != first_level_parts - 1 else num_samples
        first_level_chunks.append((start_idx, end_idx))

    # print("Индексы первого уровня:",first_level_chunks)

    # Second level: further split each part into chunks for each core
    chunks = []
    for start_idx, end_idx in first_level_chunks:
        second_level_chunk_size = (end_idx - start_idx) // num_cores
        for core_idx in range(num_cores):
            s_idx = start_idx + core_idx * second_level_chunk_size
            e_idx = start_idx + (core_idx + 1) * second_level_chunk_size if core_idx != num_cores - 1 else end_idx
            chunks.append((s_idx, e_idx, t_vect, a_vect, pulse_duration_in_samples, fs))

    # print("Индексы второго уровня:",chunks)

    # Use multiprocessing Pool to process chunks in parallel
    with Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap_unordered(generate_signal_chunk, chunks), total=len(chunks), desc="Прогресс генерации"))

    # Combine all chunks into the final signal
    signal = np.concatenate(results)

    # Print the first and last elements of the generated signal
    print(f"Первый элемент сигнала: {signal[0]}")
    print(f"Последний элемент сигнала: {signal[-1]}")

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(timestamp)
    return signal

def embed_data(final_Len, signal, start_time, fs, lowRand, highRand):
    num_samples = int(final_Len * fs)
    print("Итоговое число отсчетов:",num_samples)
    start_sample = int(start_time * fs)
    end_sample = start_sample + len(signal)
    print("Первый отсчет:",start_sample)
    print("Последний отсчет:",end_sample)
    final_signal = np.zeros(num_samples)
    prefix = "signal_P"    
    suffix = "S.bin"
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fileName = f"signal_{timestamp}.bin"
    print(f"\nЗапись в файл: {fileName}")
    
    with open(fileName, 'wb') as f:    
        # Оборачиваем цикл записи в tqdm для отображения прогресса
        for i in tqdm(range(num_samples), desc="Прогресс записи", unit="samples"):
            if start_sample <= i < end_sample:
                final_signal[i] = signal[i - start_sample]
                codes = adc_14bit(final_signal[i])
                f.write(struct.pack('H', codes))
                f.write(struct.pack('H', codes))
                f.flush()
            else:
                final_signal[i] = random.uniform(lowRand, highRand)
                codes = adc_14bit(final_signal[i])
                f.write(struct.pack('H', codes))
                f.write(struct.pack('H', codes))
                f.flush()
    
    print("Готово!")
    return final_signal   

def embed_data2(final_Len, signal, start_time, fs, lowRand, highRand):
    num_samples = int(final_Len * fs)     
    start_sample = int(start_time * fs)
    end_sample = start_sample + len(signal)
    final_signal = np.zeros(num_samples)
 
    for i in range(num_samples):
        if start_sample <= i < end_sample:
            final_signal[i] = signal[i - start_sample]
            codes = adc_14bit(final_signal[i])

        else:
            final_signal[i] = random.uniform(lowRand, highRand)
            codes = adc_14bit(final_signal[i])

    return final_signal   

if __name__ == "__main__":
    # Защита от рекурсивного запуска процессов
    from multiprocessing import freeze_support
    freeze_support()

    # Параметры сигнала
    r_Time_start = random.uniform(0, 30)
    # Записываем время начала выполнения программы
    start_timestamp = datetime.datetime.now()  # Для вывода времени начала
    timestamp = time()  # Для вычисления затраченного времени
    print(f"Начало выполнения: {start_timestamp.strftime('%Y-%m-%d %H:%M:%S')}")
    # Частота дискретизации
    fs = 50e6

    polez_signal = signal_Gen(fs,800e-6,1,10)
    # print(polez_signal)
    signalG = generate_signal(polez_signal, 10*10e-6, fs)
    def adc_14bit(voltage):
        # Преобразование напряжения в код АЦП
        code = int((voltage + 1) * 8191)
        return code
    
    dataG = embed_data(2,signalG,0,fs,-0.01,0.01)

    end_time = time()
    time_taken = end_time - timestamp
    time_delta = datetime.timedelta(seconds=time_taken)
    print(f"Затраченное время: {time_delta}")

