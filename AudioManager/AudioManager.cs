using NAudio.Wave;
using System;
using System.Numerics;
using System.Collections.Generic;

using NAudio.CoreAudioApi;

using System.Linq;

using WaveFormat = NAudio.Wave.WaveFormat;

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
        WaveFormat waveFormat;
        public DirectSound(string deviceName, int sampleRate = 48000)
        {
            _sampleRate = sampleRate;
            waveFormat = new WaveFormat(_sampleRate, 16, 2);
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
                WaveFormat = waveFormat
            };
            _waveIn.DataAvailable += OnDataAvailable;

            // 初始化缓冲区
            _waveProvider = new BufferedWaveProvider(waveFormat)
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

            // 创建 WaveFormat 对象
            
            int bitsPerSample = waveFormat.BitsPerSample;

            // 根据 bitsPerSample 进行转换
            if (bitsPerSample == 16)
            {
                float[] floatBuffer = new float[bytesRecorded / 2];
                for (int i = 0; i < floatBuffer.Length; i++)
                {
                    short value = BitConverter.ToInt16(buffer, i * 2);
                    floatBuffer[i] = value / 32768f; // 将 16 位整数转换为浮点数
                }
                onDataRecevice?.Invoke(floatBuffer);
            }
            else if (bitsPerSample == 32)
            {
                float[] floatBuffer = new float[bytesRecorded / 4];
                for (int i = 0; i < floatBuffer.Length; i++)
                {
                    float value = BitConverter.ToSingle(buffer, i * 4);
                    if (float.IsNaN(value))
                    {
                        value = 0.0f; // 处理 NaN 值
                    }
                    floatBuffer[i] = value;
                }
                onDataRecevice?.Invoke(floatBuffer);
            }
            else
            {
                throw new NotSupportedException($"不支持的位深度: {bitsPerSample}");
            }
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
        WaveFormat waveFormat;
        public WASAPI(string deviceName, int sampleRate = 48000)
        {
            _sampleRate = sampleRate;
            waveFormat = new WaveFormat(_sampleRate, 16, 2);
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
                var deviceEnumerator = new MMDeviceEnumerator();
                _device = deviceEnumerator.GetDefaultAudioEndpoint(DataFlow.Render, Role.Multimedia);


            }
            // 初始化 WASAPI 捕获
            _wasapiCapture = new WasapiCapture(_device);
            _wasapiCapture.DataAvailable += OnDataAvailable;

            // 初始化缓冲区
            _waveProvider = new BufferedWaveProvider(waveFormat)
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

            // 创建 WaveFormat 对象

            int bitsPerSample = waveFormat.BitsPerSample;
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
        WaveFormat waveFormat;
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
    public class CubicSpline
    {
        private readonly double[] a;
        private readonly double[] b;
        private readonly double[] c;
        private readonly double[] d;
        private readonly double[] x;

        public CubicSpline(double[] x, double[] y)
        {
            int n = x.Length;
            this.x = x;
            a = new double[n];
            b = new double[n];
            c = new double[n];
            d = new double[n];

            double[] h = new double[n - 1];
            double[] alpha = new double[n - 1];
            for (int i = 0; i < n - 1; i++)
            {
                h[i] = x[i + 1] - x[i];
                if (i > 0)
                {
                    alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]);
                }
            }

            double[] l = new double[n];
            double[] mu = new double[n];
            double[] z = new double[n];
            l[0] = 1;
            mu[0] = 0;
            z[0] = 0;
            for (int i = 1; i < n - 1; i++)
            {
                l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }
            l[n - 1] = 1;
            z[n - 1] = 0;
            c[n - 1] = 0;

            for (int j = n - 2; j >= 0; j--)
            {
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
                d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
                a[j] = y[j];
            }
        }

        public double Interpolate(double xVal)
        {
            int i = Array.BinarySearch(x, xVal);
            if (i < 0)
            {
                i = ~i - 1;
            }
            double dx = xVal - x[i];
            return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }
    }



    public class FFT
    {
        public enum FillDataType
        {
            SILDINGWINDOW,
            ZEROPAD,
            SILDINGWINDOW_AND_ZEROPAD,
            ZEROPAD_AND_SILDINGWINDOW,
            NONE
  
        }



        private SlidingWindow _slidingWindowLeft = null;
        private SlidingWindow _slidingWindowRight = null;
        private FillDataType _fillDataType = FillDataType.ZEROPAD_AND_SILDINGWINDOW;
        private ComplexArrayPool _complexArrayPool;
        //默认配置 是平衡的 
        private int _zeroPadSize = 8196;
        private int _slidingWindowSize = 16384;
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


        // 计算 FFT 使用 Cooley-Tukey 算法
        public Complex[] ComputeFFT(Complex[] data)
        {
            // 获取输入数据的长度
            int N = data.Length;

            // 如果数据长度小于等于 1，直接返回原数据
            // 因为长度为 1 或 0 的序列不需要进行 FFT 计算
            if (N <= 1)
                return data;

            // 创建两个数组，分别存储输入数据中的偶数项和奇数项
            Complex[] even = new Complex[N / 2];
            Complex[] odd = new Complex[N / 2];

            // 将输入数据的偶数项和奇数项分别存储到 even 和 odd 数组中
            for (int i = 0; i < N / 2; i++)
            {
                even[i] = data[2 * i];
                odd[i] = data[2 * i + 1];
            }

            // 递归调用 ComputeFFT 方法，对偶数项和奇数项分别进行 FFT 计算
            even = ComputeFFT(even);
            odd = ComputeFFT(odd);

            // 创建一个长度为 N 的数组，用于存储最终的 FFT 结果
            Complex[] result = new Complex[N];

            // 合并偶数项和奇数项的 FFT 结果
            for (int k = 0; k < N / 2; k++)
            {
                // 计算旋转因子 t
                // Complex.Exp 是计算复数的指数函数
                // -2 * Math.PI * Complex.ImaginaryOne * k / N 是旋转因子的指数部分
                Complex t = Complex.Exp(-2 * Math.PI * Complex.ImaginaryOne * k / N) * odd[k];

                // 计算结果的前半部分
                result[k] = even[k] + t;

                // 计算结果的后半部分
                result[k + N / 2] = even[k] - t;
            }

            // 返回最终的 FFT 结果
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

            //ApplyFlattopWindow(data);//平顶窗
                                     //ApplyHammingWindow(data);
                                     //ApplyKaiserWindow(data);

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
            
            float[] zeroPadLeftChannel = new float[_zeroPadSize];
            float[] zeroPadRightChannel = new float[_zeroPadSize];
            switch (_fillDataType)
            {
                case FillDataType.ZEROPAD:
                    Array.Copy(leftChannel, zeroPadLeftChannel, leftChannel.Length);
                    Array.Copy(rightChannel, zeroPadRightChannel, rightChannel.Length);
                    ProcessFFT(zeroPadLeftChannel, zeroPadRightChannel);
                    break;
                case FillDataType.SILDINGWINDOW:
                    _zeroPadSize = _slidingWindowSize;
                    _slidingWindowLeft.AddData(leftChannel);
                    _slidingWindowRight.AddData(rightChannel);
                    ProcessFFT(_slidingWindowLeft.GetWindow(), _slidingWindowRight.GetWindow());
                    break;

                case FillDataType.SILDINGWINDOW_AND_ZEROPAD:
                    if(_slidingWindowSize > _zeroPadSize){
                        throw new Exception("先进行滑动窗口尺寸不能大于零填充尺寸");
                        //break;
                    }
                    _slidingWindowLeft.AddData(leftChannel);
                    _slidingWindowRight.AddData(rightChannel);
                    Array.Copy(_slidingWindowLeft.GetWindow(), zeroPadLeftChannel, _slidingWindowLeft.GetWindow().Length);
                    Array.Copy(_slidingWindowRight.GetWindow(), zeroPadRightChannel, _slidingWindowRight.GetWindow().Length);
                    ProcessFFT(zeroPadLeftChannel, zeroPadRightChannel);
                    break;

                case FillDataType.ZEROPAD_AND_SILDINGWINDOW:
                    if (  _zeroPadSize > _slidingWindowSize)
                    {
                        throw new Exception("先进行零填充,零填充数组尺寸不能大于滑动窗口尺寸");
                        //break;
                    }
                    Array.Copy(leftChannel, zeroPadLeftChannel, leftChannel.Length);
                    Array.Copy(rightChannel, zeroPadRightChannel, rightChannel.Length);
                    _slidingWindowLeft.AddData(zeroPadLeftChannel);
                    _slidingWindowRight.AddData(zeroPadRightChannel);
                    ProcessFFT(_slidingWindowLeft.GetWindow(), _slidingWindowRight.GetWindow());
                    break;

                case FillDataType.NONE:
                    ProcessFFT(leftChannel, rightChannel);
                    break;
             
            }
        }


        public static float[] InterpolateData(float[] x, float[] y, int targetLength)
        {
            float[] interpolatedData = new float[targetLength];
            double step = (x[x.Length - 1] - x[0]) / (targetLength - 1);

            for (int i = 0; i < targetLength; i++)
            {
                double xVal = x[0] + i * step;
                int j = Array.BinarySearch(x, (float)xVal);
                if (j < 0)
                {
                    j = ~j - 1;
                }

                if (j >= x.Length - 1)
                {
                    interpolatedData[i] = y[x.Length - 1];
                }
                else
                {
                    double t = (xVal - x[j]) / (x[j + 1] - x[j]);
                    interpolatedData[i] = (float)((1 - t) * y[j] + t * y[j + 1]);
                }
            }

            // 确保频率数组的最后一个元素不为 0
            if (interpolatedData[targetLength - 1] == 0)
            {
                interpolatedData[targetLength - 1] = interpolatedData[targetLength - 2];
            }

            return interpolatedData;
        }
    }






}
