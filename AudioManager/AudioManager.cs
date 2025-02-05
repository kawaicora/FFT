using NAudio.Wave;
using System;
using System.Numerics;
using System.Collections.Generic;
using static AudioManager.FFT;
using NAudio.CoreAudioApi;
using NAudio.Wave.SampleProviders;
using System.Linq;
using NAudio.CoreAudioApi.Interfaces;
using System.Runtime.InteropServices;

namespace AudioManager
{
    public class SlidingWindow
    {
        private int windowSize;
        private float[] window;
        private int currentIndex;

        public int WindowSize
        {
            get { return windowSize; }
        }
        // Constructor
        public SlidingWindow(int size)
        {
            windowSize = size;
            window = new float[windowSize];
            currentIndex = 0;

            // Initialize window data to 0
            for (int i = 0; i < windowSize; i++)
            {
                window[i] = 0;
            }
        }
        // Add new data to the window
        public void AddData(float[] newData)
        {
            int newDataLength = newData.Length;

            if (currentIndex + newDataLength <= windowSize)
            {
                Array.Copy(newData, 0, window, currentIndex, newDataLength);
                currentIndex += newDataLength;
            }
            else
            {
                int shiftAmount = Math.Min(newDataLength, windowSize);
                Array.Copy(window, shiftAmount, window, 0, windowSize - shiftAmount);
                int remainingSpace = windowSize - shiftAmount;
                Array.Copy(newData, 0, window, remainingSpace, Math.Min(shiftAmount, newDataLength));
                currentIndex = windowSize;
            }
        }

        // Get current window data
        public float[] GetWindow()
        {
            return (float[])window.Clone();
        }
    }

    // Object Pool for Complex Arrays
    public class ComplexArrayPool
    {
        private Queue<Complex[]> _pool;
        private int _arraySize;
        private int _maxSize;  // Track the max size manually

        public ComplexArrayPool(int arraySize, int poolSize)
        {
            _arraySize = arraySize;
            _maxSize = poolSize;  // Set the maximum pool size
            _pool = new Queue<Complex[]>(poolSize);
            for (int i = 0; i < poolSize; i++)
            {
                _pool.Enqueue(new Complex[arraySize]);
            }
        }

        public Complex[] Rent()
        {
            if (_pool.Count > 0)
            {
                return _pool.Dequeue();
            }
            return new Complex[_arraySize]; // Allocate a new array if pool is empty
        }
        public void Return(Complex[] array)
        {
            if (_pool.Count < _maxSize)  // Check if the pool is not full
            {
                _pool.Enqueue(array);
            }
        }
    }

    public class DirectSound
    {
        private WaveInEvent _waveIn;
        private DirectSoundOut _directSoundOut;
        private BufferedWaveProvider _waveProvider;
        private int _sampleRate;
        private WaveInCapabilities? _device;
        public delegate void OnDataRecevice(float[] buffer);
        public OnDataRecevice? onDataRecevice;

        public DirectSound(string deviceName, int sampleRate = 48000)
        {
            _sampleRate = sampleRate;
            if (string.IsNullOrEmpty(deviceName))
            {
                _device = WaveInEvent.GetCapabilities(0); // 获取默认设备
            }
            else
            {
                _device = GetWaveInDeviceByName(deviceName);
            }

            if (_device == null)
            {
                throw new Exception("this._device is null");
            }

            // 初始化 WaveIn 捕获
            _waveIn = new WaveInEvent
            {
                DeviceNumber = GetDeviceNumberByName(deviceName),
                WaveFormat = new WaveFormat(_sampleRate, 16, 2)
            };
            _waveIn.DataAvailable += OnDataAvailable;

            // 初始化缓冲区
            _waveProvider = new BufferedWaveProvider(new WaveFormat(_sampleRate, 16, 2))
            {
                DiscardOnBufferOverflow = true
            };

            // 初始化 DirectSound 输出
            _directSoundOut = new DirectSoundOut();
            _directSoundOut.Init(_waveProvider);
        }

        private void OnDataAvailable(object sender, WaveInEventArgs e)
        {
            byte[] buffer = e.Buffer;
            int bytesRecorded = e.BytesRecorded;
            float[] floatBuffer = new float[bytesRecorded / 4]; // 假设 32 位浮点数
            for (int i = 0; i < floatBuffer.Length; i++)
            {
                floatBuffer[i] = BitConverter.ToSingle(buffer, i * 4);
            }
            onDataRecevice?.Invoke(floatBuffer);
        }

        public DirectSoundOut directSoundOut
        {
            get
            {
                return _directSoundOut;
            }
        }

        public void Stop()
        {
            _waveIn.StopRecording();
            _directSoundOut.Stop();
        }

        public void Play()
        {
            _waveIn.StartRecording();
            _directSoundOut.Play();
        }

        public static string[] GetAllDirectSoundDevices()
        {
            try
            {
                var deviceCount = WaveInEvent.DeviceCount;
                var devices = new List<string>();
                for (int i = 0; i < deviceCount; i++)
                {
                    var deviceInfo = WaveInEvent.GetCapabilities(i);
                    devices.Add(deviceInfo.ProductName);
                }
                return devices.ToArray();
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

        public WaveInCapabilities? GetWaveInDeviceByName(string deviceName)
        {
            try
            {
                var deviceCount = WaveInEvent.DeviceCount;
                for (int i = 0; i < deviceCount; i++)
                {
                    var deviceInfo = WaveInEvent.GetCapabilities(i);
                    if (deviceInfo.ProductName == deviceName)
                    {
                        return deviceInfo;
                    }
                }
                return null;
            }
            catch (Exception ex)
            {
                _ = ex;
                return null;
            }
        }

        private int GetDeviceNumberByName(string deviceName)
        {
            var deviceCount = WaveInEvent.DeviceCount;
            for (int i = 0; i < deviceCount; i++)
            {
                var deviceInfo = WaveInEvent.GetCapabilities(i);
                if (deviceInfo.ProductName == deviceName)
                {
                    return i;
                }
            }
            return 0; // 默认设备
        }
    }

    public class WASAPI
    {
        private WasapiCapture _wasapiCapture;
        private WasapiOut _wasapiOut;
        private BufferedWaveProvider _waveProvider;
        private int _sampleRate;
        private MMDevice? _device;
        public delegate void OnDataRecevice(float[] buffer);
        public OnDataRecevice? onDataRecevice;

        public WASAPI(string deviceName, int sampleRate = 48000)
        {
            _sampleRate = sampleRate;
            if (deviceName == "")
            {
                var deviceEnumerator = new MMDeviceEnumerator();
                _device = deviceEnumerator.GetDefaultAudioEndpoint(DataFlow.All, Role.Multimedia);
            }
            else
            {
                _device = GetMMDeviceByName(deviceName);
            }

            if (_device == null)
            {
                throw new Exception("this._device is null");

            }
            // 初始化 WASAPI 捕获
            _wasapiCapture = new WasapiCapture(_device);
            _wasapiCapture.DataAvailable += OnDataAvailable;

            // 初始化缓冲区
            _waveProvider = new BufferedWaveProvider(new WaveFormat(_sampleRate, 16, 2))
            {
                DiscardOnBufferOverflow = true
            };
            // 初始化 WASAPI 输出
            _wasapiOut = new WasapiOut(_device, AudioClientShareMode.Shared, false, 100);

        }

        private void OnDataAvailable(object sender, WaveInEventArgs e)
        {
            byte[] buffer = e.Buffer;
            int bytesRecorded = e.BytesRecorded;
            float[] floatBuffer = new float[bytesRecorded / 4]; // 假设 32 位浮点数
            for (int i = 0; i < floatBuffer.Length; i++)
            {
                floatBuffer[i] = BitConverter.ToSingle(buffer, i * 4);
            }
            if (onDataRecevice != null)
            {
                onDataRecevice(floatBuffer);
            }
            else
            {
                throw new Exception("onDataRecevice is null");
            }
        }

        public WasapiOut wasapiOut
        {
            get
            {
                return _wasapiOut;
            }
        }

        public void Stop()
        {
            _wasapiCapture.StopRecording();
            _wasapiOut.Stop();
        }

        public void Play()
        {
            _wasapiCapture.StartRecording();
            _wasapiOut.Play();
        }

        public static string[] GetAllWASAPIDevices()
        {
            try
            {
                var deviceEnumerator = new MMDeviceEnumerator();
                var deviceEnumeratorArray = deviceEnumerator.EnumerateAudioEndPoints(DataFlow.All, DeviceState.Active).ToArray();

                string[] devices = new string[deviceEnumeratorArray.Length];
                for (int i = 0; i < deviceEnumeratorArray.Length; i++)
                {
                    devices[i] = deviceEnumeratorArray[i].FriendlyName;
                }
                return devices;
            }
            catch (Exception ex)
            {
                throw (ex);
                //return new string[0];
            }
        }

        public MMDevice? GetMMDeviceByName(string deviceName)
        {
            try
            {
                var deviceEnumerator = new MMDeviceEnumerator();
                var deviceEnumeratorArray = deviceEnumerator.EnumerateAudioEndPoints(DataFlow.All, DeviceState.Active);
                foreach (var device in deviceEnumeratorArray)
                {
                    if (device.FriendlyName == deviceName)
                    {
                        return device;
                    }
                }
                return null;
            }
            catch (Exception ex)
            {
                _ = ex;
                return null;
            }
        }
    }



    public class ASIO
    {

        private AsioOut _asioOut;               // ASIO Output
        private BufferedWaveProvider _waveProvider;
        private int _sampleRate;
        //private string _asioDriverName;
        public delegate void OnDataRecevice(float[] buffer);
        public OnDataRecevice? onDataRecevice;
        public ASIO(string asioDriverName, int sampleRate = 48000)
        {
            _sampleRate = sampleRate;
            // Initialize ASIO driver
            _asioOut = new AsioOut(asioDriverName);

            _asioOut.AudioAvailable += OnAudioAvailable;

            // Initialize buffer
            _waveProvider = new BufferedWaveProvider(new WaveFormat(_sampleRate, 24, 2))
            {
                DiscardOnBufferOverflow = true
            };

            _asioOut.InitRecordAndPlayback(null, 2, _sampleRate);
        }

        private void OnAudioAvailable(object sender, AsioAudioAvailableEventArgs e)
        {
            float[] buffer = new float[_sampleRate];
            int size = e.GetAsInterleavedSamples(buffer);
            Array.Resize(ref buffer, size);
            if (onDataRecevice != null)
            {
                onDataRecevice(buffer);
            }



        }
        public AsioOut asioOut
        {
            get
            {
                return _asioOut;
            }
        }

        public void Stop() => _asioOut.Stop();
        public void Play() => _asioOut.Play();

        public void Play(int inputChannelOffset)
        {
            _asioOut.InputChannelOffset = inputChannelOffset;
            Play();
        }

        public string[] GetAllInputChannel()
        {
            string[] inputChannelList = new string[_asioOut.DriverInputChannelCount];
            for (int i = 0; i < inputChannelList.Length; i++)
            {
                inputChannelList[i] = _asioOut.AsioInputChannelName(i);

            }
            return inputChannelList;
        }
        public string[] GetAllOutputChannel()
        {
            string[] outputChannelList = new string[_asioOut.DriverOutputChannelCount];
            for (int i = 0; i < outputChannelList.Length; i++)
            {
                outputChannelList[i] = _asioOut.AsioOutputChannelName(i);

            }
            return outputChannelList;
        }
        public static string[] GetAllAsioDevice()
        {
            return AsioOut.GetDriverNames(); ;
        }
    }
    public class FFT
    {
        public enum FillDataType
        {
            SILDINGWINDOW,
            NONE
  
        }



        private SlidingWindow _slidingWindowLeft = null;
        private SlidingWindow _slidingWindowRight = null;
        private FillDataType _fillDataType = FillDataType.SILDINGWINDOW;
        private ComplexArrayPool _complexArrayPool;
        private int _zeroPadSize = 16384;
        private int _slidingWindowSize = 4096;
        private int _sampleRate = 48000;

        Complex[]? leftComplex = null;
        Complex[]? rightComplex = null;

        public struct SpectrumData
        {
            public float[] frequencies;
            public float[] amplitudesLeft;
            public float[] amplitudesRight;
            public float volumeLeft;
            public float volumeRight;
        }
        private SpectrumData _spectrumData;

        // Constructor
        public FFT(int sampleRate = 48000)
        {

            _sampleRate = sampleRate;
          
            _spectrumData = new SpectrumData();


        }
        public FillDataType fillDataType
        {
            set
            {
                _fillDataType = value;
            }
            get
            {
                return _fillDataType;
            }
        }
        public int zeroPadSize
        {
            set
            {
                _zeroPadSize = value;
            }
        }

        public int slidingWindowSize
        {
            set
            {
                _slidingWindowSize = value;
            }
        }
        ~FFT()
        {

        }

        public SpectrumData spectrumData { get { return _spectrumData; } }



        // Compute FFT using Cooley-Tukey algorithm
        public Complex[] ComputeFFT(Complex[] data)
        {
            int N = data.Length;

            if (N <= 1)
                return data;

            Complex[] even = new Complex[N / 2];
            Complex[] odd = new Complex[N / 2];
            for (int i = 0; i < N / 2; i++)
            {
                even[i] = data[2 * i];
                odd[i] = data[2 * i + 1];
            }

            even = ComputeFFT(even);
            odd = ComputeFFT(odd);

            Complex[] result = new Complex[N];
            for (int k = 0; k < N / 2; k++)
            {
                Complex t = Complex.Exp(-2 * Math.PI * Complex.ImaginaryOne * k / N) * odd[k];
                result[k] = even[k] + t;
                result[k + N / 2] = even[k] - t;
            }

            return result;
        }

        // 计算音量
        private float CalculateVolume(float[] samples, int channelIndex)
        {
            float volume = 0.0f;
            for (int i = channelIndex; i < samples.Length; i += 2)
            {
                volume += samples[i] * samples[i];
            }
            volume = (float)Math.Sqrt(volume / (samples.Length / 2));
            return volume;
        }

        // Convert FFT results to magnitude
        public float[] ConvertToMagnitude(Complex[] fftResult)
        {
            int N = fftResult.Length;
            float[] magnitude = new float[N];

            for (int i = 0; i < N; i++)
            {
                magnitude[i] = (float)fftResult[i].Magnitude;
            }

            return magnitude;
        }
        // Apply a window function (Hamming window)

        public void ApplyGaussianWindow(Complex[] data)
        {
            int N = data.Length;
            float sigma = N / 8.0f;  // Adjust this value based on your needs
            float mean = N / 2.0f;

            for (int i = 0; i < N; i++)
            {
                float window = (float)Math.Exp(-0.5 * Math.Pow((i - mean) / sigma, 2));
                data[i] *= window;
            }
        }
        public void ApplyKaiserWindow(Complex[] data, float beta = 14.0f)
        {
            int N = data.Length;

            for (int i = 0; i < N; i++)
            {
                float window = I0(beta * (float)Math.Sqrt(1 - Math.Pow(2.0 * i / (N - 1) - 1, 2))) / I0(beta);
                data[i] *= window;
            }
        }



        // 凯泽窗函数
        public void ApplyKaiserWindow(Complex[] data)
        {
            int N = data.Length;
            float beta = 14.0f; // 调整 beta 值以控制旁瓣衰减

            // 应用凯泽窗函数
            for (int i = 0; i < N; i++)
            {
                float window = I0(beta * (float)Math.Sqrt(1 - Math.Pow(2.0 * i / (N - 1) - 1, 2))) / I0(beta);
                data[i] *= window;
            }
        }

        // Modified Bessel function of the first kind (I0) for Kaiser window
        public float I0(float x)
        {
            float sum = 1.0f;
            float term = 1.0f;
            int n = 1;

            while (term > 1e-15)
            {
                term *= x * x / (4.0f * n * n);
                sum += term;
                n++;
            }

            return sum;
        }
        public void ApplyFlattopWindow(Complex[] data)
        {
            int N = data.Length;

            // Flattop窗的系数
            float a0 = 1.0f;
            float a1 = 1.93f;
            float a2 = 1.29f;
            float a3 = 0.388f;
            float a4 = 0.032f;

            // 遍历每个样本，应用Flattop窗
            for (int i = 0; i < N; i++)
            {
                // Flattop窗的公式
                float window = a0 - a1 * (float)Math.Cos(2 * Math.PI * i / (N - 1))
                                + a2 * (float)Math.Cos(4 * Math.PI * i / (N - 1))
                                - a3 * (float)Math.Cos(6 * Math.PI * i / (N - 1))
                                + a4 * (float)Math.Cos(8 * Math.PI * i / (N - 1));

                // 将窗函数乘到数据上
                data[i] *= window;
            }
        }
        public void ApplyHammingWindow(Complex[] data)
        {
            int N = data.Length;

            // 遍历每个样本，应用汉明窗
            for (int i = 0; i < N; i++)
            {
                // 汉明窗的公式
                float window = 0.54f - 0.46f * (float)Math.Cos(2 * (float)Math.PI * i / (N - 1));

                // 将窗函数乘到数据上
                data[i] *= window;
            }
        }



        // 执行FFT
        public Complex[] PerformFFT(Complex[] data, int targetLength)
        {

            ApplyFlattopWindow(data);//平顶窗
            return ComputeFFT(data);//计算FFT
        }


        // FFT结果复数转float
        public float[] ConvertToDecibels(Complex[] fftResult)
        {
            int N = fftResult.Length;
            float[] decibels = new float[N];

            for (int i = 0; i < N; i++)
            {
                float magnitude = (float)fftResult[i].Magnitude;
                if (magnitude > 0)
                {
                    decibels[i] = 20 * (float)Math.Log10(magnitude);
                }
                else
                {
                    decibels[i] = -100;
                }
            }

            return decibels;
        }

        // 计算点对应频率
        public float[] CalculateFrequencies(int N, int fs)
        {
            float[] frequencies = new float[N];
            for (int i = 0; i < N; i++)
            {
                frequencies[i] = i * fs / N;
            }
            return frequencies;
        }
        // 平滑数据
        private void SmoothData(float[] data, int smoothingFactor)
        {
            if (smoothingFactor <= 1) return;  // No smoothing needed if factor is 1 or less

            float[] smoothedData = new float[data.Length];
            float sum = 0f;

            for (int i = 0; i < data.Length; i++)
            {
                sum += data[i];

                // 将总和保持在“smoothingFactor”元素的窗口内
                if (i >= smoothingFactor)
                {
                    sum -= data[i - smoothingFactor];
                }

                // 计算移动平均线
                smoothedData[i] = sum / Math.Min(i + 1, smoothingFactor);
            }

            // 将平滑后的数据复制回原始数组
            Array.Copy(smoothedData, data, data.Length);
        }

        // 处理 FFT
        private void ProcessFFT(float[] leftChannel, float[] rightChannel)
        {
            if (_spectrumData.frequencies == null)
            {
                _spectrumData.frequencies = new float[(leftChannel.Length + rightChannel.Length) / 2];
                _spectrumData.frequencies = CalculateFrequencies((leftChannel.Length + rightChannel.Length) / 2, _sampleRate);
            }
            if (_spectrumData.amplitudesLeft == null)
            {
                _spectrumData.amplitudesLeft = new float[leftChannel.Length];
            }
            if (_spectrumData.amplitudesRight == null)
            {
                _spectrumData.amplitudesRight = new float[rightChannel.Length];
            }
            if (_complexArrayPool == null)
            {
                _complexArrayPool = new ComplexArrayPool((leftChannel.Length + rightChannel.Length) / 2, 10);  // Pool size can be adjusted

            }

            // Initialize the ComplexArrayPool to avoid allocating new arrays frequently


            leftComplex = _complexArrayPool.Rent();
            rightComplex = _complexArrayPool.Rent();

            for (int i = 0; i < (leftChannel.Length + rightChannel.Length) / 2; i++)
            {
                leftComplex[i] = new Complex(leftChannel[i], 0);
                rightComplex[i] = new Complex(rightChannel[i], 0);
            }

            Complex[] leftFFT = PerformFFT(leftComplex, 0);
            Complex[] rightFFT = PerformFFT(rightComplex, 0);

            float[] leftDecibels = ConvertToMagnitude(leftFFT);
            float[] rightDecibels = ConvertToMagnitude(rightFFT);

            for (int i = 0; i < _spectrumData.frequencies.Length; i++)
            {
                _spectrumData.amplitudesLeft[i] = leftDecibels[i];
                _spectrumData.amplitudesRight[i] = rightDecibels[i];
            }

            _complexArrayPool.Return(leftComplex);
            _complexArrayPool.Return(rightComplex);
        }

        public void OnDataRecevice(float[] buffer)
        {

            int size = buffer.Length;

            float[] leftChannel = new float[size / 2];
            float[] rightChannel = new float[size / 2];

            for (int i = 0; i < size / 2; i++)
            {
                leftChannel[i] = buffer[i * 2];
                rightChannel[i] = buffer[i * 2 + 1];
            }

            _spectrumData.volumeLeft = CalculateVolume(leftChannel, 0);
            _spectrumData.volumeRight = CalculateVolume(rightChannel, 1);
            if (_slidingWindowLeft == null) {
                _slidingWindowLeft = new SlidingWindow(_slidingWindowSize);
            }
            if (_slidingWindowRight == null) {
                _slidingWindowRight = new SlidingWindow(_slidingWindowSize);
            }
            
            float[] extCL = new float[_zeroPadSize];
            float[] extCR = new float[_zeroPadSize];
            switch (_fillDataType)
            {
                case FillDataType.SILDINGWINDOW:
                    _zeroPadSize = _slidingWindowSize;
                    _slidingWindowLeft.AddData(leftChannel);
                    _slidingWindowRight.AddData(rightChannel);
                 
                    ProcessFFT(_slidingWindowLeft.GetWindow(), _slidingWindowRight.GetWindow());
              
               
                    break;
                case FillDataType.NONE:
                    ProcessFFT(leftChannel, rightChannel);
                    break;
             
            }
        }
    }
}
