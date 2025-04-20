using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Shapes;
using Avalonia.Media;
using Avalonia.Markup.Xaml;
using AudioManager;
using static AudioManager.FFT;
using System;
using System.Threading.Tasks;
using System.Threading;
using Avalonia.Data;
using Avalonia.Collections;
namespace AvaloniaFFTApp.Views
{
    public partial class MainView : UserControl
    {
        enum AudioType { 
            ASIO,
            DirectSound,
            WASAPI
        }
        private AudioType audioType; 
        private FFT fft;
        private ASIO asio;
        private WASAPI wasapi;
        DirectSound directSound;
        private float volumeLeft;
        private float volumeRight;
        private float[] frequencies;
        private float[] amplitudesLeft;
        private float[] amplitudesRight;
        private int sampleRate = 48000;
        double yScale = 4;
        //int zeroPadSize = 8196;
        //int windowSize = 16384;
        // 设定要绘制的频率范围
        const double minFrequency = 0;
        const double maxFrequency = 10000;
        bool drawRectangles = false;  //绘制矩形还是波形
        double rectangleWidth = 6;
        private IBinding iBinding;
        private const int consoleWidth = 800;
        private const int consoleHeight = 600;
        public MainView()
        {
            InitializeComponent();
            audioType = AudioType.DirectSound;

            var audioTypeComboBox = this.Find<ComboBox>("AudioTypeSelectBox");

            if (audioTypeComboBox == null) {
                return;
            }
            audioTypeComboBox.ItemsSource = null;
            if (audioTypeComboBox.ItemsSource == null)
            {
                audioTypeComboBox.ItemsSource = new AvaloniaList<string>();
            }
            var items = audioTypeComboBox.ItemsSource as AvaloniaList<string>;

            if (items == null) {
                return;
            }
            items.Clear();
            items.Add("DirectSound");
            items.Add("WASAPI");
            items.Add("ASIO");

            // 处理窗口大小变化
            this.GetObservable(BoundsProperty).Subscribe(OnBoundsChanged);
        }
        private void AudioTypeSelectBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {

            if (sender is ComboBox comboBox)
            {
                var selectedItem = comboBox.SelectedItem as string;

                switch (selectedItem) {
                    case "ASIO":
                        audioType = AudioType.ASIO;
                        break;
                    case "WASAPI":
                        audioType = AudioType.WASAPI;
                        break;
                    case "DirectSound":
                        audioType = AudioType.DirectSound;
                        break;
                }   

                var deviceComboBox = this.Find<ComboBox>("DeviceSelectBox");
                if (deviceComboBox == null)
                {
                    return;
                }
                deviceComboBox.ItemsSource = null;
                // 绑定 ComboBox 的 ItemsSource 属性到 ViewModel 的 ComboBoxItems 属性
                if (deviceComboBox.ItemsSource == null)
                {
                    deviceComboBox.ItemsSource = new AvaloniaList<string>();
                }
                
                var items = deviceComboBox.ItemsSource as AvaloniaList<string>;
                
                if (items != null)
                {
                    items.Clear();
                    switch (audioType)
                    {
                        case AudioType.ASIO:
                            foreach (var item in ASIO.GetAllAsioDevice())
                            {
                                //添加音频设备选项
                                items.Add(item);
                            }
                            break;
                        case AudioType.WASAPI:
                            foreach (var item in WASAPI.GetAllWASAPIDevices())
                            {
                                //添加音频设备选项
                                items.Add(item);

                            }
                            break;
                        case AudioType.DirectSound:
                            foreach (var item in DirectSound.GetAllDirectSoundDevices())
                            {
                                items.Add(item);
                            }
                            break;
                    }

                }

            }
        }

        private void DeviceSelectBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (sender is ComboBox comboBox)
            {
                string? selectedItem = comboBox.SelectedItem as string;
                fft = new FFT(sampleRate); //初始化FFT
                //不定义使用默认配置
                fft.addWindow = false;
                fft.scaleSpectrumData = true;
                switch (audioType)
                {
                    //选择不同的设备类型 执行不同的初始化与获取下一步所需数据
                    case AudioType.ASIO:
                       
                        if (selectedItem == null) {
                            return;
                        }
                        asio = new ASIO(selectedItem, sampleRate);//初始化ASIO
                        asio.onDataRecevice = fft.OnDataRecevice;  //绑定FFT方法到ASIO委托变量


                        var channelComboBox = this.Find<ComboBox>("ChannelSelectBox");
                        // 绑定 ComboBox 的 ItemsSource 属性到 ViewModel 的 ComboBoxItems 属性
                        if (channelComboBox != null && channelComboBox.ItemsSource == null)
                        {
                            channelComboBox.ItemsSource = new AvaloniaList<string>();
                            var items = channelComboBox.ItemsSource as AvaloniaList<string>;
                            if (items == null) {

                                return;
                            }
                            items.Clear();
                            var outputChannels = asio.GetAllOutputChannel();
                            foreach (var item in outputChannels)
                            {
                                if (!string.IsNullOrEmpty(item) && items != null)
                                {
                                    items.Add(item);
                                }

                            }
                            channelComboBox.IsVisible = true;
                            var channelSelectText = this.Find<TextBlock>("ChannelSelectText");
                            if (channelSelectText != null)
                            {
                                channelSelectText.IsVisible = true;
                            }
                        }


                        break;
                    case AudioType.WASAPI:
                   
                        
                        wasapi = new WASAPI(selectedItem, sampleRate);
                        wasapi.onDataRecevice += fft.OnDataRecevice;
                        wasapi.Stop();
                        wasapi.Play();
                        StartDrawTask();
                        break;
                    case AudioType.DirectSound:
                        directSound = new DirectSound(selectedItem, sampleRate);
                        directSound.onDataRecevice += fft.OnDataRecevice;
                        directSound.Stop();
                        directSound.Play();
                        StartDrawTask();
                        break;
                }
            }
        }
        private void ChannelSelectBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (sender is ComboBox comboBox)
            {
                var selectedItem = comboBox.SelectedItem as string;
                // 这里添加选定选项后的回调逻辑
       
                var currentChannelIndexOffset = Array.IndexOf(asio.GetAllOutputChannel(), selectedItem);
                asio.Stop();
                asio.Play(currentChannelIndexOffset);
                StartDrawTask();



            }
        }

        private void StartDrawTask() {
            // 启动一个后台任务来更新 FFT 数据
            Task.Run(async () =>
            {
                while (true)
                {
                    try
                    {
                        UpdateData(fft.spectrumData);
                        Thread.Sleep(1000 / 15);
                    }
                    catch (Exception ex)
                    {
                    }
                }
            });
        }
        private void InitializeComponent()
        {
            AvaloniaXamlLoader.Load(this);
        }
        private void UpdateData(SpectrumData spectrumData)
        {
            volumeLeft = spectrumData.volumeLeft;
            volumeRight = spectrumData.volumeRight;
            frequencies = spectrumData.frequencies;
    
            amplitudesLeft = spectrumData.amplitudesLeft;
            amplitudesRight = spectrumData.amplitudesRight;
            // 在 UI 线程更新图形
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                DrawGraph();
            });
        }
        private void OnBoundsChanged(Rect bounds)
        {
            DrawGraph();
        }


   


        private void DrawGraph()
        {
            var canvas = this.FindControl<Canvas>("canvas");
            if (canvas != null)
            {
                if (frequencies != null && amplitudesLeft != null)
                {
                    //DrawFFTBarChart(canvas, frequenciesLeft, amplitudesLeft, minFrequency, maxFrequency, rectangleWidth, yScale);
                    DrawFFT(canvas, frequencies, amplitudesLeft, minFrequency, maxFrequency, rectangleWidth, yScale);
                }
            }
        }




        private void DrawFFT(Canvas canvas, float[] frequencies, float[] amplitudes, double minFrequency, double maxFrequency, double rectangleWidth, double yScale)
        {
            if (canvas == null || frequencies == null || amplitudes == null)
                return;
            canvas.Children.Clear();
            if (frequencies == null || amplitudesLeft == null)
                return; // 如果数据无效，直接返回



            int dataLength = frequencies.Length;

            // 计算 X 轴缩放比例：画布宽度 / (最大频率 - 最小频率)
            double xScale = canvas.Bounds.Width / (maxFrequency - minFrequency);

            // Y 轴缩放比例：画布高度 / 2

            // 对振幅数据进行插值处理
            //int targetLength = dataLength * 4; // 目标点数，可以根据需要调整
            //float[] interpolatedAmplitudesLeft = InterpolateData(frequenciesLeft, amplitudesLeft, targetLength);
            //float[] interpolatedFrequenciesLeft = InterpolateData(frequenciesLeft, frequenciesLeft, targetLength);


            if (drawRectangles)
            {
                // 绘制方块
                for (int i = 0; i < frequencies.Length; i++)
                {
                    // 计算当前数据点的频率值
                    double frequency = frequencies[i];

                    if (frequency >= minFrequency && frequency <= maxFrequency)
                    {
                        // 映射频率值到 X 坐标
                        double x = (frequency - minFrequency) * xScale;

                        // 映射振幅值到 Y 坐标
                        double height = amplitudesLeft[i] * yScale;
                        double y = canvas.Bounds.Height - height;

                        var rectangle = new Rectangle
                        {
                            Fill = Brushes.Red,
                            Width = rectangleWidth,
                            Height = height,
                            Margin = new Thickness(x, y, 0, 0)
                        };
                        canvas.Children.Add(rectangle);
                    }
                }
            }
            else
            {
                // 创建用于绘制频率图和振幅图的 Polyline
                var amplitudeLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = new Avalonia.Collections.AvaloniaList<Point>()
                };

                // 绘制振幅线条：筛选并映射到画布
                for (int i = 0; i < frequencies.Length; i++)
                {
                    // 计算当前数据点的频率值
                    double frequency = frequencies[i];

                    if (frequency >= minFrequency && frequency <= maxFrequency)
                    {
                        // 映射频率值到 X 坐标
                        double x = (frequency - minFrequency) * xScale;

                        // 映射振幅值到 Y 坐标
                        double y = canvas.Bounds.Height - 5 - amplitudesLeft[i] * yScale;

                        amplitudeLine.Points.Add(new Point(x, y));
                    }
                }
                // 添加到画布
                canvas.Children.Add(amplitudeLine);
            }
        }

    }
    }

