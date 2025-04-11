import numpy as np
import matplotlib.pyplot as plt
import sys
import struct
from tqdm import tqdm

def calculate_sample_indices(file_size, fs, num_channels=2, sample_fraction=0.1):
    """Вычисление индексов необходимых данных"""
    bytes_per_sample = 2 * num_channels
    total_samples = file_size // bytes_per_sample
    num_samples_needed = int(total_samples * sample_fraction)
    indices = np.linspace(0, total_samples-1, num_samples_needed, dtype=int)
    return np.sort(indices), total_samples

def read_selected_data(file_path, indices, num_channels=2):
    """
    Чтение только выбранных данных из файла для двух каналов
    :return: (channel_1_data, channel_2_data)
    """
    bytes_per_sample = 2 * num_channels
    channel_1_data = []
    channel_2_data = []
    
    with open(file_path, 'rb') as f:
        for idx in tqdm(indices, desc="Чтение данных", unit="sample"):
            f.seek(idx * bytes_per_sample)
            raw_data = f.read(bytes_per_sample)
            sample = struct.unpack(f'{num_channels}H', raw_data)
            channel_1_data.append(sample[0])  # Первый канал
            channel_2_data.append(sample[1])  # Второй канал
            
    return np.array(channel_1_data, dtype=np.uint16), np.array(channel_2_data, dtype=np.uint16)

def adc_to_voltage(codes):
    """Преобразование кодов АЦП в напряжение"""
    return np.clip(codes / 8191 - 1, -1, 1)

def plot_signals(voltage_1, voltage_2, times, fs):
    """Визуализация сигналов для двух каналов"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    
    # График для первого канала
    ax1.plot(times, voltage_1, 'b-', lw=0.7)
    ax1.set_title('Канал 1')
    ax1.set_ylabel('Напряжение')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # График для второго канала
    ax2.plot(times, voltage_2, 'r-', lw=0.7)
    ax2.set_title('Канал 2')
    ax2.set_xlabel('Время, сек')
    ax2.set_ylabel('Напряжение')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.show()

def main(file_path, sample_fraction=0.1):
    FS = 50e6  # Частота дискретизации
    NUM_CHANNELS = 2  # Кол-во каналов
    
    try:
        # Получаем размер файла
        file_size = None
        with open(file_path, 'rb') as f:
            f.seek(0, 2)  # Переходим в конец файла
            file_size = f.tell()  # Размер файла в байтах
        
        if not file_size:
            raise ValueError("Файл пустой или недоступен")
        
        # Вычисляем необходимые индексы
        print("Вычисление индексов...")
        needed_indices, total_samples = calculate_sample_indices(
            file_size, FS, NUM_CHANNELS, sample_fraction)
        
        print(f"Общее количество сэмплов: {total_samples}")
        print(f"Выбрано для визуализации: {len(needed_indices)} сэмплов")
        
        # Чтение только нужных данных
        print("Чтение выбранных данных...")
        channel_1_data, channel_2_data = read_selected_data(file_path, needed_indices, NUM_CHANNELS)
        
        # Преобразование в напряжение
        print("Обработка данных...")
        voltage_1 = adc_to_voltage(channel_1_data)
        voltage_2 = adc_to_voltage(channel_2_data)
        
        # Временные метки
        times = needed_indices / (FS * NUM_CHANNELS / 2)
        
        # Визуализация
        print(f"Отображение графиков...")
        plot_signals(voltage_1, voltage_2, times, FS)
        
    except Exception as e:
        print(f"Ошибка: {str(e)}")

if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     print("Использование: python visualize.py <файл> [канал=0]")
    #     sys.exit(1)
        
    # file_path = sys.argv[1]
    # sample_fraction = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    file_path = "build/Desktop_Qt_5_15_2_MinGW_32_bit-Debug/11-04-2025_14-32-17.bin"
    sample_fraction=0.01
    
    try:
        main(file_path, sample_fraction)
    except Exception as e:
        print(f"Ошибка: {str(e)}")
        sys.exit(1)