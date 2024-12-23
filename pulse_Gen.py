import numpy as np
import struct
import random
import datetime
import time

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
    return out

def generate_signal(signal, pulse_duration, fs):
    t_vect = np.array([k for k in signal.keys()])
    a_vect = np.array([v for v in signal.values()])
    # Calculate the indices of the signal corresponding to the start of each pulse
    indices1 = (t_vect * fs).astype(int)
    # Calculate the length of the signal array
    num_samples = int((t_vect[-1] + pulse_duration) * fs)

    # Create an array of zeros for the signal
    signal = np.zeros(num_samples)

    # Calculate the duration of the pulse in samples
    pulse_duration_in_samples = int(pulse_duration * fs)

    # Loop over each sample in the signal
    for i in range(num_samples):
        # Loop over each pulse
        for j in range(len(indices1)):
            # Check if the sample is in the rising edge of the pulse
            if indices1[j] <= i < int(indices1[j] + 0.2 * pulse_duration_in_samples):
                # Calculate the slope and intercept of the rising edge
                x1 = indices1[j]
                y1 = 0
                x2 = int(indices1[j] + 0.2 * pulse_duration_in_samples)
                y2 = a_vect[j]
                slope = (y2 - y1) / (x2 - x1)
                intercept = y1 - slope * x1                
                # Calculate the signal value for the current sample
                signal[i] = slope * i + intercept+random.uniform(0.01*a_vect[j], 0.02*a_vect[j])                
                
            # Check if the sample is in the falling edge of the pulse
            elif int(indices1[j] + 0.8 * pulse_duration_in_samples) <= i < int(indices1[j] + pulse_duration_in_samples):
                # Calculate the slope and intercept of the falling edge
                x1 = indices1[j] + 0.8 * pulse_duration_in_samples
                y1 = a_vect[j]
                x2 = int(indices1[j] + pulse_duration_in_samples)
                y2 = 0
                slope = (y2 - y1) / (x2 - x1)
                intercept = y1 - slope * x1

                # Calculate the signal value for the current sample
                signal[i] = slope * i + intercept+random.uniform(0.01*a_vect[j], 0.02*a_vect[j])                
            # Check if the sample is in the steady state of the pulse
            elif int(indices1[j] + 0.2 * pulse_duration_in_samples) <= i < int(indices1[j] + 0.8 * pulse_duration_in_samples):
                # Set the signal value to the amplitude of the pulse
                signal[i] = a_vect[j]+random.uniform(0.01*a_vect[j], 0.02*a_vect[j])                
            # Set the signal value to zero for all other samples
            elif signal[i] != 0:
                continue
            else:
                signal[i] = random.uniform(-0.01, 0.01)                
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(timestamp)
    return signal

def embed_data(final_Len, signal, start_time, fs, lowRand, highRand):
    num_samples = int(final_Len * fs)     
    start_sample = int(start_time * fs)
    end_sample = start_sample + len(signal)
    final_signal = np.zeros(num_samples)
    prefix = "signal_P"    
    suffix = "S.bin"
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fileName = "{}{}{}{}".format(prefix, final_Len, suffix, timestamp)
    progress_counter = 0
    print(fileName)
    with open(fileName, 'wb') as f:    
        for i in range(num_samples):
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
    f.close()
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
    
r_Time_start = random.uniform(0, 30)
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(timestamp)
# Частота дискретизации
fs = 50e6

polez_signal = signal_Gen(fs,800e-6,1,10)
print(polez_signal)
signalG = generate_signal(polez_signal, 10*10e-6, fs)
def adc_14bit(voltage):
    # Преобразование напряжения в код АЦП
    code = int((voltage + 1) * 8191)
    return code

dataG = embed_data(5,signalG,1,fs,-0.01,0.01)

# end_time = time.time()
# time_taken = end_time - timestamp
# time_delta = datetime.timedelta(seconds=time_taken)
# print(f"Time taken: {time_delta}")
