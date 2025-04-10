import numpy as np
import matplotlib.pyplot as plt
import sys
import struct
from tqdm import tqdm

def read_binary_file(file_path, num_channels=2):
    """Чтение бинарного файла с АЦП данными"""
    with open(file_path, 'rb') as f:
        data = f.read()
    total_samples = len(data) // (2 * num_channels)
    return np.array(struct.unpack(f'{total_samples*num_channels}H', data), 
                    dtype=np.uint16).reshape(-1, num_channels)

def adc_to_voltage(codes):
    """Преобразование кодов АЦП в напряжение"""
    return np.clip(codes / 8191 - 1, -1, 1)

def sample_data(data, fraction=0.1):
    """Равномерная выборка данных"""
    num_samples = int(len(data) * fraction)
    indices = np.linspace(0, len(data)-1, num_samples, dtype=int)
    return data[indices], indices

def plot_signal(voltage, times, fs):
    """Визуализация сигнала"""
    plt.figure(figsize=(14, 7))
    plt.plot(times, voltage, 'b-', lw=0.7)
    plt.title(f'Визуализация сигнала ({len(voltage)} точек, {fs/1e6} МГц)')
    plt.xlabel('Время, сек')
    plt.ylabel('Напряжение')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

def main(file_path, channel=0, sample_fraction=0.1):
    # Конфигурация
    FS = 50e6  # Частота дискретизации из основного скрипта
    NUM_CHANNELS = 2  # Кол-во каналов в файле
    
    # Чтение данных
    print("Чтение файла...")
    adc_data = read_binary_file(file_path, NUM_CHANNELS)
    
    # Выбор канала
    channel_data = adc_data[:, channel]
    
    # Преобразование в напряжение
    print("Обработка данных...")
    voltage = adc_to_voltage(channel_data)
    
    # Выборка
    print("Выборка данных...")
    sampled_voltage, indices = sample_data(voltage, sample_fraction)
    
    # Временные метки
    times = indices / FS
    
    # Визуализация
    print(f"Отображение {len(sampled_voltage)} точек")
    plot_signal(sampled_voltage, times, FS)

if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     print("Использование: python visualize.py <файл> [канал=0]")
    #     sys.exit(1)
        
    # file_path = sys.argv[1]
    # channel = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    file_path = "signal_20250410_085106.bin"
    channel = 0
    
    try:
        main(file_path, channel)
    except Exception as e:
        print(f"Ошибка: {str(e)}")
        sys.exit(1)