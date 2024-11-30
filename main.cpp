#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <fstream>
#include <iomanip> // Для std::setprecision
#include <random> // Для генерации случайных чисел

// Общий мьютекс для синхронизации вывода
std::mutex mtx;

// Функция для заполнения части вектора
void fill_vector_part(std::vector<double>& iteration_vec, int global_start, int local_start, int count, double fs) {
    for (int i = 0; i < count; ++i) {
        iteration_vec[local_start + i] = static_cast<double>(global_start + i) / fs;
    }
}

// Функция для умножения части вектора на 4 с учетом диапазона и генерации шума
void multiply_vector_part(std::vector<double>& iteration_vec, int local_start, int count, double signal_start, double signal_end) {
    // Генератор случайных чисел
    std::random_device rd; // Получаем случайное начальное значение
    std::mt19937 gen(rd()); // Инициализация генератора
    std::uniform_real_distribution<double> dis(-0.01, 0.01); // Распределение для генерации шума

    for (int i = 0; i < count; ++i) {
        // Получаем текущее значение
        double value = iteration_vec[local_start + i];
//        iteration_vec[local_start + i] *= 1;
        // Умножаем на 1000, если значение находится в заданном диапазоне
        if (value >= signal_start && value <= signal_end) {
//            std::cout << "Generating sinusoid at time: " << value << " " << iteration_vec[0] << iteration_vec [count] << std::endl;
//            iteration_vec[local_start + i] = 1 * sin(2 * M_PI * 1000 * value); // Умножаем на 1000
            iteration_vec[local_start + i] = 1;
        } else {
            iteration_vec[local_start + i] = dis(gen); // Добавляем шум
        }
    }
}

void writeToFileLE(std::ofstream& file, size_t num_samples, const std::vector<double>& signal_vect) {
    const int blockSize = 1024;

    for (size_t i = 0; i < num_samples; i += blockSize) {
        size_t block_end = std::min(i + blockSize, num_samples);
        for (size_t j = i; j < block_end; ++j) {
            // Преобразуем значение в int16_t и кодируем в Little Endian
            int16_t value = static_cast<int16_t>(std::round((signal_vect[j] + 1) * 8191));

            // Переставляем байты для Little Endian
            uint8_t bytes[2];
            bytes[0] = static_cast<uint8_t>(value & 0xFF);        // Младший байт
            bytes[1] = static_cast<uint8_t>((value >> 8) & 0xFF);  // Старший байт

            // Записываем байты в файл дважды (16-битный формат)
            file.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
            file.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));  // Дублируем для стерео/двухканальной записи
        }
    }
}

int main() {

    const double total_time = 10.0;
    const double fs_p = 50e6;
    const int total_elements = static_cast<int>(total_time * fs_p);  // Общее количество элементов
    const int num_threads = 5;  // Количество потоков
    const double elements_per_iteration = total_elements * 0.01;  // Элементов за итерацию
    const double elements_per_thread = elements_per_iteration / num_threads;  // Количество элементов на поток
    const int bar_width = 20;  // Ширина прогресс-бара

    // Определяем границы полезного сигнала
    double signal_start = 1; // Начало полезного сигнала (пример)
    double signal_end = 5;   // Конец полезного сигнала (пример)

    // Открываем файл для записи
    std::string filename  = "generated_data_" + std::to_string(total_time) + "s.bin";
    std::ofstream output_file(filename, std::ios::binary | std::ios::app);
    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open file for writing!" << std::endl;
        return 1;
    }

//    std::ofstream output_file1("generated_data.txt", std::ios::out | std::ios::trunc);
//    if (!output_file.is_open()) {
//        std::cerr << "Error: Unable to open file for writing!" << std::endl;
//        return 1;
//    }

//    std::ofstream output_file2("generated_dataNSTR.txt", std::ios::out | std::ios::trunc);
//    if (!output_file.is_open()) {
//        std::cerr << "Error: Unable to open file for writing!" << std::endl;
//        return 1;
//    }

    // Засекаем время перед началом обработки
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int iteration = 0; iteration < total_elements / elements_per_iteration; ++iteration) {
        // Промежуточный вектор для текущей итерации
        std::vector<double> iteration_vec(static_cast<int>(elements_per_iteration));

        // Глобальный стартовый индекс для текущей итерации
        int global_start_index = iteration * elements_per_iteration;

        // Запускаем потоки для заполнения вектора значениями
        std::vector<std::thread> threads;

        for (int i = 0; i < num_threads; ++i) {
            int local_start_index = i * elements_per_thread;  // Локальный индекс в iteration_vec
            int thread_global_start = global_start_index + local_start_index;  // Глобальный старт
            threads.emplace_back(fill_vector_part, std::ref(iteration_vec), thread_global_start, local_start_index, elements_per_thread, fs_p);
        }

        // Ждем завершения всех потоков
        for (auto& t : threads) {
            t.join();
        }

//        std::cout << iteration_vec.at(5) << " " << signal_start << " " << signal_end << std::endl;

        // Теперь умножаем элементы на 4 или 1000 с помощью потоков
        threads.clear(); // Очищаем вектор потоков для повторного использования

        for (int i = 0; i < num_threads; ++i) {
            int local_start_index = i * elements_per_thread;
            threads.emplace_back(multiply_vector_part, std::ref(iteration_vec), local_start_index, elements_per_thread, signal_start, signal_end);
        }

        // Ждем завершения всех потоков
        for (auto& t : threads) {
            t.join();
        }

        // Синхронизируем запись результатов в файл
        {
            std::lock_guard<std::mutex> lock(mtx);

            // Вызов функции записи в файл
            writeToFileLE(output_file, iteration_vec.size(), iteration_vec);

//            for (const auto& val : iteration_vec) {
//                output_file1 << val << " ";
//            }

//            for (const auto& val : iteration_vec) {
//                output_file2 << val << "\n";
//            }
        }

        // Вычисляем время, прошедшее с начала выполнения
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = current_time - start_time;
        int minutes = static_cast<int>(elapsed_time.count()) / 60;
        int seconds = static_cast<int>(elapsed_time.count()) % 60;

        // Выводим процент завершения в виде прогресс-бара
        double progress = static_cast<double>(iteration + 1) / (total_elements / elements_per_iteration);
        int pos = static_cast<int>(bar_width * progress); // Позиция заполненной части

        // Заполняем прогресс-бар символами
        std::string bar(bar_width, ' '); // Изначально бар пуст
        for (int i = 0; i < pos; ++i) {
            bar[i] = '#'; // Заполненная часть
        }

        // Меняем последний символ в зависимости от остатка от деления на 5
        if (pos > 0 && pos < bar_width) {
            switch ((iteration + 1) % 5) {
                case 1: bar[pos] = '#'; break;
                case 2: bar[pos] = ' '; break;
                case 3: bar[pos] = '#'; break;
                case 4: bar[pos] = ' '; break;
            }
        }

        // Форматируем вывод прогресс-бара
        std::cout << "[" << bar << "] " << std::fixed << std::setprecision(1) << (progress * 100) << "% "
                  << "(" << minutes << "m " << seconds << "s)" << "\r"; // Отображение процента и времени
        std::cout.flush(); // Сбрасываем буфер вывода
    }

    // Засекаем время после обработки
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;

    // Выводим время обработки вектора
    std::cout << "\nProcessed vector size: " << total_elements << ", processing time: " << duration.count() << " seconds" << std::endl;

    // Закрываем файл
    output_file.close();
    // Закрываем файл
//    output_file1.close();
//    output_file2.close();

    return 0;
}
