#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <fstream>
#include <iomanip> // Для setprecision
#include <random> // Для генерации случайных чисел
#include <sstream>
#include <cmath>

using namespace std;
using namespace chrono;

// Общий мьютекс для синхронизации вывода
mutex mtx;

/*void signal_Gen(double fs,
                double period,
                double Amp,
                int step,
                vector<double>& x_vect,
                vector<double>& y_vect) {
    double end_time_sec = (period * step) / 2;
    double start_time_sec = 0;
    double f = 1 / (period * step);
    int num_samples = static_cast<int>(round((end_time_sec - start_time_sec) * fs));

    vector<double> time(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        time[i] = start_time_sec + i / fs;
    }

    vector<double> signal(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        signal[i] = Amp * sin(2 * M_PI * f * time[i]);
    }

    double per = (period * step) / (step * 2 + 1);
    vector<int> indices;
    for (int i = 0; i < num_samples; i += static_cast<int>(round(per * fs))) {
        indices.push_back(i);
    }

    vector<double> signal_selected;
    for (int i : indices) {
        if (i < num_samples) {
            signal_selected.push_back(signal[i]);
        }
    }

    int c = step + 1;
    y_vect = vector<double>(signal_selected.begin(), signal_selected.begin() + c);
    x_vect.resize(c); // Убедимся, что размер вектора x_vect соответствует y_vect
    for (int i = 0; i < c; ++i) {
        x_vect[i] = (i * ((period * step) / step) * 100) * fs;
    }
}

// Функция для заполнения части вектора
void fill_vector_partO(vector<double>& iteration_vec,
                      int global_start,
                      int local_start,
                      int count,
                      double fs) {
    for (int i = 0; i < count; ++i) {
        iteration_vec[local_start + i] = static_cast<double>(global_start + i) / fs;
    }
}

void fill_vector_part(vector<double>& iteration_vec,
                      vector<double>& x_vect,
                      vector<double>& y_vect,
                      int global_start,
                      int local_start,
                      int count,
                      double pulse_duartion_in_samples)
{
    // Генератор случайных чисел
    random_device rd; // Получаем случайное начальное значение
    mt19937 gen(rd()); // Инициализация генератора
    uniform_real_distribution<double> dis(-0.01, 0.01); // Распределение для генерации шума
    for (int i = 0; i < count; ++i) {
        for(int j = 0; j < x_vect.size(); ++j){
            if (x_vect[j] <= global_start + local_start + i && global_start + local_start + i < x_vect[j] + 0.2 * pulse_duartion_in_samples) {
                uint64_t x1 = x_vect[j];
                double y1 = 0;
                uint64_t x2 = x_vect[j] + 0.2 * pulse_duartion_in_samples;
                double y2 = y_vect[j];
                double slope = (y2 - y1) / (x2 - x1);
                double intercept = y1 - slope * x1;
                iteration_vec[local_start + i] = slope * i + intercept + dis(gen);
                // cout << "rising" << endl;
            }
            else if(x_vect[j]+0.8*pulse_duartion_in_samples <= global_start + local_start + i < x_vect[j]+pulse_duartion_in_samples){
                uint64_t x1 = x_vect[j]+0.8*pulse_duartion_in_samples;
                double y1 = y_vect[j];
                uint64_t x2 = x_vect[j]+pulse_duartion_in_samples;
                double y2 = 0;
                double slope = (y2 - y1) / (x2 - x1);
                double intercept = y1 - slope * x1;
                iteration_vec[local_start + i] = slope * i + intercept + dis(gen);
                // cout << "falling" << endl;
            }
            else if(x_vect[j]+0.2*pulse_duartion_in_samples <= global_start + local_start + i < x_vect[j]+0.8*pulse_duartion_in_samples){
                iteration_vec[local_start + i] = y_vect[j] + dis(gen);
                // cout << "pulse" << endl;
            }
            // else if()
            else {
                iteration_vec[local_start + i] = dis(gen);
            }
        }
    }
}

void generate_signal(vector<double>& signal, const vector<double>& t_vect, const vector<double>& a_vect, int pulse_duration_in_samples, int global_start_index, int local_start_index, int elements_per_thread) {
    // Инициализация генератора случайных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> noise(-0.01, 0.01);

    // Вычисление индексов начала импульсов
    vector<int> indices1(t_vect.size());
    for (size_t i = 0; i < t_vect.size(); ++i) {
        indices1[i] = static_cast<int>(t_vect[i]);
    }

    // Цикл по элементам текущего потока
    for (int i = local_start_index; i < local_start_index + elements_per_thread; ++i) {
        int global_index = global_start_index + i;

        // Цикл по каждому импульсу
        for (size_t j = 0; j < indices1.size(); ++j) {
            if (indices1[j] <= global_index && global_index < indices1[j] + pulse_duration_in_samples) {
                if (global_index < indices1[j] + 0.2 * pulse_duration_in_samples) {
                    // Расчет наклона и пересечения для фронта импульса
                    double x1 = indices1[j];
                    double y1 = 0;
                    double x2 = indices1[j] + 0.2 * pulse_duration_in_samples;
                    double y2 = a_vect[j];
                    double slope = (y2 - y1) / (x2 - x1);
                    double intercept = y1 - slope * x1;
                    signal[i] = slope * global_index + intercept + noise(gen);
                } else if (global_index >= indices1[j] + 0.8 * pulse_duration_in_samples) {
                    // Расчет наклона и пересечения для спада импульса
                    double x1 = indices1[j] + 0.8 * pulse_duration_in_samples;
                    double y1 = a_vect[j];
                    double x2 = indices1[j] + pulse_duration_in_samples;
                    double y2 = 0;
                    double slope = (y2 - y1) / (x2 - x1);
                    double intercept = y1 - slope * x1;
                    signal[i] = slope * global_index + intercept + noise(gen);
                } else {
                    // Установить значение сигнала на амплитуду импульса
                    signal[i] = a_vect[j] + noise(gen);
                }
            } else if (signal[i] != 0) {
                continue;
            } else {
                signal[i] = noise(gen);
            }
        }
    }
}

// Функция для умножения части вектора на 4 с учетом диапазона и генерации шума
void multiply_vector_part(vector<double>& iteration_vec,
                          int local_start,
                          int count,
                          double signal_start,
                          double signal_end,
                          double fs) {
    // Генератор случайных чисел
    random_device rd; // Получаем случайное начальное значение
    mt19937 gen(rd()); // Инициализация генератора
    uniform_real_distribution<double> dis(-0.01, 0.01); // Распределение для генерации шума

    for (int i = 0; i < count; ++i) {
        // Получаем текущее значение
        double value = iteration_vec[local_start + i];
//        iteration_vec[local_start + i] *= 1;
        // Умножаем на 1000, если значение находится в заданном диапазоне
        if (value >= signal_start && value <= signal_end) {
//            cout << "Generating sinusoid at time: " << value << " " << iteration_vec[0] << iteration_vec [count] << endl;
//            iteration_vec[local_start + i] = 1 * sin(2 * M_PI * 1000 * value); // Умножаем на 1000
            iteration_vec[local_start + i] = 1;
        } else {
            iteration_vec[local_start + i] = dis(gen); // Добавляем шум
        }
    }
}

// Функция для преобразования напряжения в код АЦП
int adc_14bit(double voltage) {
    return static_cast<int>((voltage + 1) * 8191);
}

// Функция для встраивания данных в сигнал
// void embed_data2(vector<double>& final_signal, const vector<double>& signal, double start_sample, double end_sample, double lowRand, double highRand, int global_start_index, int local_start_index, int elements_per_thread) {
//     // Инициализация генератора случайных чисел
//     random_device rd;
//     mt19937 gen(rd());
//     uniform_real_distribution<> rand_dist(lowRand, highRand);

//     // Цикл по элементам текущего потока
//     for (int i = local_start_index; i < local_start_index + elements_per_thread; ++i) {
//         int global_index = global_start_index + i;

//         if (start_sample <= global_index && global_index < end_sample) {
//             // Встраиваем сигнал
//             final_signal[i] = signal[global_index - start_sample];
//             // adc_14bit(final_signal[i]); // Вызов функции adc_14bit
//         } else {
//             // Генерируем случайное значение
//             final_signal[i] = rand_dist(gen);
//             // adc_14bit(final_signal[i]); // Вызов функции adc_14bit
//         }
//     }
// }

void writeToFileLE(ofstream& file,
                   size_t num_samples,
                   const vector<double>& signal_vect) {
    const int blockSize = 1024;

    for (size_t i = 0; i < num_samples; i += blockSize) {
        size_t block_end = min(i + blockSize, num_samples);
        for (size_t j = i; j < block_end; ++j) {
            // Преобразуем значение в int16_t и кодируем в Little Endian
            int16_t value = static_cast<int16_t>(round((signal_vect[j] + 1) * 8191));

            // Переставляем байты для Little Endian
            uint8_t bytes[2];
            bytes[0] = static_cast<uint8_t>(value & 0xFF);        // Младший байт
            bytes[1] = static_cast<uint8_t>((value >> 8) & 0xFF);  // Старший байт

            // Записываем байты в файл дважды (16-битный формат)
            file.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
            file.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));  // Дублируем для стерео/двухканальной записи
            // file.flush();
        }        
    }
}

void set_new_indeces(vector<double>& x_vect, double pulse_duartion_in_samples) {
    for (int i = 0; i < x_vect.size(); ++i) {
        x_vect[i] = x_vect[i] + pulse_duartion_in_samples;
    }
}*/

string generateFileName() {
    // Получаем текущее время
    auto now = chrono::system_clock::now();
    time_t time = chrono::system_clock::to_time_t(now);

    // Преобразуем в локальное время
    tm localTime;
#ifdef _WIN32
    localtime_s(&localTime, &time); // Windows
#else
    localtime_r(&time, &localTime); // Linux/Unix
#endif

    // Форматируем имя файла
    ostringstream fileName;
    fileName << put_time(&localTime, "%d-%m-%Y_%H-%M-%S");
    return fileName.str();
}

void show_progress_percentage(const time_point<steady_clock>& start_time, int current, int total) {
    auto current_time = steady_clock::now();
    chrono::duration<double> elapsed_time = current_time - start_time;
    int minutes = static_cast<int>(elapsed_time.count()) / 60;
    int seconds = static_cast<int>(elapsed_time.count()) % 60;
    const int bar_width = 20;  // Ширина прогресс-бара
    double progress = static_cast<double>(current + 1) / total * 100.0;
    int pos = static_cast<int>(bar_width * progress); // Позиция заполненной части

    // Заполняем прогресс-бар символами
    string bar(bar_width, ' '); // Изначально бар пуст
    for (int i = 0; i < pos; ++i) {
        bar[i] = '#'; // Заполненная часть
    }

    // Меняем последний символ в зависимости от остатка от деления на 5
    if (pos > 0 && pos < bar_width) {
        switch ((current + 1) % 5) {
        case 1: bar[pos] = '#'; break;
        case 2: bar[pos] = ' '; break;
        case 3: bar[pos] = '#'; break;
        case 4: bar[pos] = ' '; break;
        }
    }

    // Форматируем вывод прогресс-бара
    cout << " [" << bar << "] " << fixed << setprecision(1) << (progress * 100) << "% "
         << "(" << minutes << "m " << seconds << "s)" << "\r"; // Отображение процента и времени
}

// Функция для генерации случайного числа в диапазоне
double random_uniform(double low, double high) {
    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<> dis(low, high);
    return dis(gen);
}

// Преобразование напряжения в 14-битный код АЦП
uint16_t adc_14bit(double voltage) {
    return static_cast<uint16_t>((voltage + 1) * 8191);
}

// Генерация сигнала (аналог signal_Gen)
vector<pair<double, double>> signal_Gen(double fs, double period, double Amp, int step) {
    vector<pair<double, double>> out;
    double end_time_sec = (period * step) / 2;
    double start_time_sec = 0;
    double f = 1 / (period * step);
    int num_samples = static_cast<int>((end_time_sec - start_time_sec) * fs);

    vector<double> time(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        time[i] = start_time_sec + i / fs;
    }

    vector<double> signal(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        signal[i] = Amp * sin(2 * M_PI * f * time[i]);
    }

    double per = (period * step) / (step * 2 + 1);
    vector<int> indices;
    for (int i = 0; i < num_samples; i += static_cast<int>(round(per * fs))) {
        indices.push_back(i);
    }

    vector<double> signal_selected;
    for (int i : indices) {
        if (i < num_samples) {
            signal_selected.push_back(signal[i]);
        }
    }

    int c = step + 1;
    vector<double> y_vect = vector<double>(signal_selected.begin(), signal_selected.begin() + c);
    vector<double> x_vect; x_vect.resize(c); // Убедимся, что размер вектора x_vect соответствует y_vect
    for (int i = 0; i < c; ++i) {
        x_vect[i] = (i * ((period * step) / step) * 100);
        out.emplace_back(x_vect[i], y_vect[i]);
    }

    return out;
}

// Обработка части сигнала
vector<double> generate_signal_chunk(int start_idx, int end_idx, const vector<double>& t_vect, const vector<double>& a_vect,
                                     int pulse_duration_in_samples, double fs) {
    vector<double> chunk_signal(end_idx - start_idx);
    vector<int> indices1(t_vect.size());

    for (size_t i = 0; i < t_vect.size(); ++i) {
        indices1[i] = static_cast<int>(t_vect[i] * fs);
    }

    for (int i = start_idx; i < end_idx; ++i) {
        int idx_local = i - start_idx;
        for (size_t j = 0; j < indices1.size(); ++j) {
            // Rising edge
            if (indices1[j] <= i && i < indices1[j] + static_cast<int>(0.2 * pulse_duration_in_samples)) {
                double x1 = indices1[j], y1 = 0;
                double x2 = indices1[j] + 0.2 * pulse_duration_in_samples, y2 = a_vect[j];
                double slope = (y2 - y1) / (x2 - x1);
                double intercept = y1 - slope * x1;
                chunk_signal[idx_local] = slope * i + intercept + random_uniform(0.01 * a_vect[j], 0.02 * a_vect[j]);
            }
            // Falling edge
            else if (indices1[j] + static_cast<int>(0.8 * pulse_duration_in_samples) <= i &&
                     i < indices1[j] + pulse_duration_in_samples) {
                double x1 = indices1[j] + 0.8 * pulse_duration_in_samples, y1 = a_vect[j];
                double x2 = indices1[j] + pulse_duration_in_samples, y2 = 0;
                double slope = (y2 - y1) / (x2 - x1);
                double intercept = y1 - slope * x1;
                chunk_signal[idx_local] = slope * i + intercept + random_uniform(0.01 * a_vect[j], 0.02 * a_vect[j]);
            }
            // Steady state
            else if (indices1[j] + static_cast<int>(0.2 * pulse_duration_in_samples) <= i &&
                     i < indices1[j] + static_cast<int>(0.8 * pulse_duration_in_samples)) {
                chunk_signal[idx_local] = a_vect[j] + random_uniform(0.01 * a_vect[j], 0.02 * a_vect[j]);
            }

            else if (chunk_signal[idx_local] != 0) continue;

            // Noise outside pulses
            else if (chunk_signal[idx_local] == 0) {
                chunk_signal[idx_local] = random_uniform(-0.01, 0.01);
            }
        }
    }

    return chunk_signal;
}

// Основная функция генерации сигнала
vector<double> generate_signal(const vector<pair<double, double>>& signal, double pulse_duration, double fs) {
    vector<double> t_vect, a_vect;
    for (const auto& p : signal) {
        t_vect.push_back(p.first);
        a_vect.push_back(p.second);
    }

    int num_samples = static_cast<int>((t_vect.back() + pulse_duration) * fs);
    int pulse_duration_in_samples = static_cast<int>(pulse_duration * fs);

    vector<double> final_signal(num_samples, 0.0);
    int num_threads = thread::hardware_concurrency();
    vector<thread> threads;
    vector<vector<double>> results(num_threads);

    int chunk_size = num_samples / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        int start_idx = i * chunk_size;
        int end_idx = (i == num_threads - 1) ? num_samples : start_idx + chunk_size;
        threads.emplace_back([&, i, start_idx, end_idx]() {
            results[i] = generate_signal_chunk(start_idx, end_idx, t_vect, a_vect, pulse_duration_in_samples, fs);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    auto start_time = steady_clock::now(); // Время начала сборки
    for (int i = 0; i < num_threads; ++i) {
        int start_idx = i * chunk_size;
        int end_idx = (i == num_threads - 1) ? num_samples : start_idx + chunk_size;
        for (int j = start_idx; j < end_idx; ++j) {
            final_signal[j] = results[i][j - start_idx];
        }
        // show_progress_percentage(start_time, i, num_threads);
    }

    return final_signal;
}

// Запись данных в бинарный файл
void embed_data(double final_Len, const vector<double>& signal, double start_time, double fs, double lowRand, double highRand) {
    int num_samples = static_cast<int>(final_Len * fs);
    int start_sample = static_cast<int>(start_time * fs);
    int end_sample = start_sample + signal.size();

    string fileDate = generateFileName();
    string filename  = fileDate + ".bin";

    cout << fileDate << endl;

    ofstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Error file opening!" << endl;
        return;
    }
    auto start_time_progress = steady_clock::now(); // Время начала записи
    for (int i = 0; i < num_samples; ++i) {
        double value = (start_sample <= i && i < end_sample) ? signal[i - start_sample] : random_uniform(lowRand, highRand);
        uint16_t code = adc_14bit(value);
        file.write(reinterpret_cast<const char*>(&code), sizeof(uint16_t));
        file.write(reinterpret_cast<const char*>(&code), sizeof(uint16_t));
        // show_progress_percentage(start_time_progress, i, num_samples);
    }

    file.close();
    cout << "Writing ended!" << endl;
}

int main(int argc, char* argv[]) {
    // Значения по умолчанию
    double total_time = 2.0;
    double signal_start = 0.0;

    // Проверяем наличие аргументов командной строки
    if (argc >= 2) {
        try {
            total_time = stod(argv[1]); // Первый аргумент - total_time
        } catch (const exception& e) {
            cerr << "Invalid total_time argument. Using default value: " << total_time << endl;
        }
    }
    if (argc >= 3) {
        try {
            signal_start = stod(argv[2]); // Второй аргумент - signal_start
        } catch (const exception& e) {
            cerr << "Invalid signal_start argument. Using default value: " << signal_start << endl;
        }
    }

    double amplitude = 1.0;
    double frequency = 1000.0;

    if (argc >= 4) amplitude = stod(argv[4]);
    if (argc >= 5) frequency = stod(argv[5]);

    cout << "Amplitude: " << amplitude << ", Frequency: " << frequency << " Hz" << endl;

    // Вывод параметров для проверки
    cout << "Parameters:" << endl;
    cout << "Total Time: " << total_time << " seconds" << endl;
    cout << "Signal Start: " << signal_start << " seconds" << endl;

    // const double total_time = 5.0;
    const double fs_p = 50e6;
    const int total_elements = static_cast<int>(total_time * fs_p);  // Общее количество элементов
    const int num_threads = 5;  // Количество потоков
    const double elements_per_iteration = total_elements * 0.01;  // Элементов за итерацию
    const double elements_per_thread = elements_per_iteration / num_threads;  // Количество элементов на поток
    const int bar_width = 20;  // Ширина прогресс-бара

    const double pulse_duartion = 10*10e-6;
    int pulse_duartion_in_samples = pulse_duartion * fs_p;
    const double period = 800e-6; // Пример значения
    const int step = 10; // Пример значения

    double lowRand = -0.01;
    double highRand = 0.01;

    // Генерация сигнала
    auto start = chrono::high_resolution_clock::now();
    double start_time = chrono::duration<double>(chrono::high_resolution_clock::now() - start).count();

    auto polez_signal = signal_Gen(fs_p, 800e-6, amplitude, 10);
    vector<double> t_vect, a_vect;
    for (const auto& p : polez_signal) {
        t_vect.push_back(p.first);
        a_vect.push_back(p.second);
    }

    vector<double> signalG = generate_signal(polez_signal, pulse_duartion, fs_p);
    embed_data(total_time, signalG, signal_start, fs_p, lowRand, highRand);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    cout << "Elapsed time: " << elapsed.count() << " secounds" << endl;

    return 0;
}
