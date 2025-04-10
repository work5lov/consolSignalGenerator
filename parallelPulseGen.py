import numpy as np
import struct
import random
import datetime
import time
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import multiprocessing as mp

def signal_Gen(fs, period, Amp, step):
    out = {}
    end_time_sec = (period * step) / 2
    start_time_sec = 0
    f = 1 / (period * step)
    num_samples = round((end_time_sec - start_time_sec) * fs)
    
    # Генерация временных меток с прогрессом
    time = []
    for i in tqdm(range(num_samples), desc="Генерация временных меток"):
        time_val = start_time_sec + i / fs
        time.append(time_val)
        # time.sleep(0.0001)
    
    # Генерация сигнала с прогрессом
    signal = []
    for t in tqdm(time, desc="Генерация синусоиды"):
        val = Amp * np.sin(2 * np.pi * f * t)
        signal.append(val)
        # time.sleep(0.0001)
    
    per = (period * step) / (step * 2 + 1)
    indices = [i for i in range(0, num_samples, round(per * fs))]
    signal_selected = [signal[i] for i in indices]
    c = step + 1
    y_vect = signal_selected[:c]
    x_vect = [i * ((period * step) / step) * 100 for i in range(c)]
    
    for i in tqdm(range(len(x_vect)), desc="Формирование выходного словаря"):
        out[x_vect[i]] = y_vect[i]
    
    return out

# def process_signal_chunk(args):
#     (chunk_indices, indices1, a_vect, pulse_samples, fs) = args
#     signal_part = {}
    
#     for i in chunk_indices:
#         value = 0.0
#         for j in range(len(indices1)):
#             start = indices1[j]
#             rise_end = int(start + 0.2 * pulse_samples)
#             fall_start = int(start + 0.8 * pulse_samples)
#             end = start + pulse_samples

#             if start <= i < rise_end:
#                 x1, y1 = start, 0
#                 x2, y2 = rise_end, a_vect[j]
#                 slope = (y2 - y1)/(x2 - x1)
#                 val = slope*(i - x1) + y1 + random.uniform(0.01*a_vect[j], 0.02*a_vect[j])
#             elif fall_start <= i < end:
#                 x1, y1 = fall_start, a_vect[j]
#                 x2, y2 = end, 0
#                 slope = (y2 - y1)/(x2 - x1)
#                 val = slope*(i - x1) + y1 + random.uniform(0.01*a_vect[j], 0.02*a_vect[j])
#             elif rise_end <= i < fall_start:
#                 val = a_vect[j] + random.uniform(0.01*a_vect[j], 0.02*a_vect[j])
#             else:
#                 val = random.uniform(-0.01, 0.01)
            
#             if val != 0:
#                 value = val
        
#         signal_part[i] = value
    
#     return signal_part

# def generate_signal(signal_dict, pulse_duration, fs):
#     t_vect = np.array([k for k in signal_dict.keys()])
#     a_vect = np.array([v for v in signal_dict.values()])
#     indices1 = (t_vect * fs).astype(int)
#     num_samples = int((t_vect[-1] + pulse_duration) * fs)
#     pulse_samples = int(pulse_duration * fs)
    
#     # Разбиваем на блоки для параллельной обработки
#     chunk_size = max(1, num_samples // (4 * cpu_count()))
#     chunks = [list(range(i, min(i + chunk_size, num_samples))) 
#              for i in range(0, num_samples, chunk_size)]
    
#     # Подготовка параметров для процессов
#     params = [(chunk, indices1, a_vect, pulse_samples, fs) for chunk in chunks]
    
#     # Параллельная обработка с прогресс-баром
#     with Pool(processes=cpu_count()) as pool:
#         results = list(tqdm(
#             pool.imap_unordered(process_signal_chunk, params),
#             total=len(chunks),
#             desc="Параллельная генерация импульсов"
#         ))
    
#     # Сборка результата
#     signal = np.zeros(num_samples)
#     for part in results:
#         for idx, val in part.items():
#             signal[idx] = val
    
#     return signal

def adc_14bit(voltage):
    code = int((voltage + 1) * 8191)
    return max(0, min(16383, code))

# def generate_block(params):
#     (block_start, block_duration, fs, signal_info, 
#      low, high) = params
    
#     signal_array, start_sample, end_sample = signal_info
#     num_samples = int(block_duration * fs)
#     buffer = bytearray()
    
#     for i in range(num_samples):
#         global_sample = int(block_start * fs) + i
#         if start_sample <= global_sample < end_sample:
#             value = signal_array[global_sample - start_sample]
#         else:
#             value = random.uniform(low, high)
        
#         code = adc_14bit(value)
#         buffer.extend(struct.pack('H', code)*2)
        
#     return buffer

# def generate_and_write_signal(filename, total_duration, signal_array, 
#                              start_time, fs, low=-0.01, high=0.01,
#                              block_size=0.5):
#     start_sample = int(start_time * fs)
#     end_sample = start_sample + len(signal_array)
#     signal_info = (signal_array, start_sample, end_sample)
    
#     num_blocks = int(np.ceil(total_duration / block_size))
#     params = []
    
#     for block_idx in range(num_blocks):
#         block_start = block_idx * block_size
#         current_duration = min(block_size, total_duration - block_start)
#         params.append((
#             block_start,
#             current_duration,
#             fs,
#             signal_info,
#             low,
#             high
#         ))
    
#     with Pool(processes=cpu_count()) as pool, \
#          tqdm(total=num_blocks, desc="Запись блоков", unit="block") as pbar:
        
#         with open(filename, 'wb') as f:
#             for buffer in pool.imap_unordered(generate_block, params):
#                 f.write(buffer)
#                 f.flush()
#                 pbar.update(1)

# if __name__ == "__main__":
#     # Параметры генерации
#     fs = 50e6  # 50 МГц
#     total_duration = 5  # Общая длительность записи (сек)
#     start_time = 1  # Начало вставки сигнала (сек)
    
#     # Генерация полезного сигнала с прогрессом
#     print("Начало генерации исходного сигнала")
#     polez_signal = signal_Gen(fs, 800e-6, 1, 10)
#     signalG = generate_signal(polez_signal, 10*10e-6, fs)
    
#     # Генерация имени файла
#     timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
#     filename = f"signal_{timestamp}.bin"
    
#     # Параллельная генерация и запись
#     start = time.time()
#     generate_and_write_signal(
#         filename=filename,
#         total_duration=total_duration,
#         signal_array=signalG,
#         start_time=start_time,
#         fs=fs,
#         low=-0.01,
#         high=0.01,
#         block_size=0.5
#     )
#     print(f"\nВремя выполнения: {time.time() - start:.2f} сек")

def generate_signal_chunk(params):
    """
    Генерирует часть сигнала для указанного временного диапазона
    """
    (start_idx, end_idx, pulse_duration, fs, 
     t_vect, a_vect, indices1, pulse_samples) = params

    chunk = np.zeros(end_idx - start_idx)
    
    for i in range(start_idx, end_idx):
        value = 0.0
        for j in range(len(indices1)):
            start = indices1[j]
            if i < start or i >= start + pulse_samples:
                continue
                
            rise_end = start + int(0.2 * pulse_samples)
            fall_start = start + int(0.8 * pulse_samples)
            end = start + pulse_samples

            if i < rise_end:
                x1, y1 = start, 0
                x2, y2 = rise_end, a_vect[j]
                slope = (y2 - y1)/(x2 - x1)
                val = slope*(i - x1) + y1 + random.uniform(0.01*a_vect[j], 0.02*a_vect[j])
            elif i >= fall_start:
                x1, y1 = fall_start, a_vect[j]
                x2, y2 = end, 0
                slope = (y2 - y1)/(x2 - x1)
                val = slope*(i - x1) + y1 + random.uniform(0.01*a_vect[j], 0.02*a_vect[j])
            else:
                val = a_vect[j] + random.uniform(0.01*a_vect[j], 0.02*a_vect[j])
            
            if val != 0:
                value = val
                
        chunk[i - start_idx] = value if value != 0 else random.uniform(-0.01, 0.01)
    
    return chunk

def generate_and_write_signal_pipeline(
    filename, total_duration, signal_params,
    start_time, fs, low=-0.01, high=0.01,
    block_size=0.5, chunk_duration=0.1
):
    # Распаковываем параметры сигнала
    (pulse_duration, t_vect, a_vect, indices1) = signal_params
    
    # Глобальные параметры
    pulse_samples = int(pulse_duration * fs)
    start_sample = int(start_time * fs)  # Начало сигнала в сэмплах
    end_sample = int((start_time + total_duration) * fs)  # Конец сигнала в сэмплах
    total_samples = end_sample - start_sample
    
    # Размер чанка в сэмплах
    chunk_samples = int(chunk_duration * fs)
    
    # Прогресс-бар
    with tqdm(total=total_samples // 2, desc="Генерация и запись", 
              unit='sample', unit_scale=True) as pbar:
        with Pool(cpu_count()) as pool, open(filename, 'wb') as f:
            for chunk_start in range(0, total_samples, chunk_samples):
                # Определяем границы чанка
                chunk_end = min(chunk_start + chunk_samples, total_samples)
                current_chunk_size = chunk_end - chunk_start
                
                # Генерируем параметры для чанка
                chunk_params = (
                    chunk_start, chunk_end, pulse_duration, fs,
                    t_vect, a_vect, indices1, pulse_samples
                )
                
                # Асинхронная генерация чанка
                async_result = pool.apply_async(generate_signal_chunk, (chunk_params,))
                chunk = async_result.get()
                
                # Преобразование в байты
                buffer = bytearray()
                for value in chunk:
                    code = adc_14bit(value)
                    buffer.extend(struct.pack('H', code) * 2)  # Запись дважды
                
                # Запись в файл
                f.write(buffer)
                f.flush()
                
                # Обновление прогресса
                pbar.update(current_chunk_size // 2)  # Делим на 2, так как данные дублируются

if __name__ == "__main__":
    fs = 50e6  # 50 МГц
    total_duration = 30  # сек
    start_time = 10  # сек
    
    # Генерация параметров сигнала
    polez_signal = signal_Gen(fs, 800e-6, 1, 10)
    t_vect = np.array([k for k in polez_signal.keys()])
    a_vect = np.array([v for v in polez_signal.values()])
    indices1 = (t_vect * fs).astype(int)
    pulse_duration = 10*10e-6
    
    # Собираем параметры для передачи
    signal_params = (pulse_duration, t_vect, a_vect, indices1)
    
    # Генерация имени файла
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"signal_{timestamp}.bin"
    
    # Запуск конвейера
    start = time.time()
    generate_and_write_signal_pipeline(
        filename=filename,
        total_duration=total_duration,
        signal_params=signal_params,
        start_time=start_time,
        fs=fs,
        low=-0.01,
        high=0.01,
        block_size=0.5,
        chunk_duration=0.1  # Размер чанка для обработки (сек)
    )
    print(f"\nВремя выполнения: {time.time() - start:.2f} сек")