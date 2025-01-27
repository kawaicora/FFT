using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Shapes;
using Avalonia.Input;
using Avalonia.Media;
using Avalonia.Markup.Xaml;
using Avalonia.Rendering;
using Avalonia.VisualTree;
using Avalonia.Platform;
using AudioManager;
using static AudioManager.FFT;
using System;
using System.Threading.Tasks;
using System.Threading;
using AvaloniaFFTApp.ViewModels;
using Avalonia.Data;
using Avalonia.Collections;
using Avalonia.Controls.Templates;
using DynamicData;
namespace AvaloniaFFTApp.Views
{
    public partial class MainView : UserControl
    {
        enum AudioType { 
            ASIO,
            WASAPI
        }
        private AudioType audioType; 
        private FFT fft;
        private ASIO asio;
        private WASAPI wasapi;
        private float[] frequenciesLeft;
        private float[] frequenciesRight;
        private float[] amplitudesLeft;
        private float[] amplitudesRight;
        private int sampleRate = 48000;
        double yScale = 4;
        int zeroPadSize = 8196;
        int windowSize = 4096;
        // 设定要绘制的频率范围
        const double minFrequency = 0;
        const double maxFrequency = 24000;
        bool drawRectangles = false;
        double rectangleWidth = 6;
        private IBinding iBinding;
        private const int consoleWidth = 800;
        private const int consoleHeight = 600;
        public MainView()
        {
            InitializeComponent();
            audioType = AudioType.WASAPI;

            var audioTypeComboBox = this.Find<ComboBox>("AudioTypeSelectBox");

            if (audioTypeComboBox == null) {
                return;
            }

            if (audioTypeComboBox.ItemsSource == null)
            {
                audioTypeComboBox.ItemsSource = new AvaloniaList<string>();
            }
            var items = audioTypeComboBox.ItemsSource as AvaloniaList<string>;

            if (items == null) {
                return;
            }
            items.Clear();
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
                }

                var deviceComboBox = this.Find<ComboBox>("DeviceSelectBox");
                if (deviceComboBox == null)
                {
                    return;
                }
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
                fft.slidingWindowSize = windowSize;
                fft.zeroPadSize = zeroPadSize;
                fft.fillDataType = FillDataType.SILDINGWINDOW;
                switch (audioType)
                {

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
            frequenciesLeft = spectrumData.frequencies;
            frequenciesRight = spectrumData.frequencies;
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
                // 清空画布
                canvas.Children.Clear();

                if (frequenciesLeft == null || amplitudesLeft == null)
                    return; // 如果数据无效，直接返回

                int dataLength = frequenciesLeft.Length;
                
                // 计算 X 轴缩放比例：画布宽度 / (最大频率 - 最小频率)
                double xScale = canvas.Bounds.Width / (maxFrequency - minFrequency);

                // Y 轴缩放比例：画布高度 / 2
                

                if (drawRectangles)
                {
                    // 绘制方块
                    for (int i = 0; i < dataLength; i++)
                    {
                        // 计算当前数据点的频率值
                        double frequency = frequenciesLeft[i];

                        if (frequency >= minFrequency && frequency <= maxFrequency)
                        {
                            // 映射频率值到 X 坐标
                            //double x = (frequency - minFrequency) * xScale;

                            double x = (frequency - minFrequency ) * xScale ;
                            // 映射振幅值到 Y 坐标
                            double height = amplitudesLeft[i] * yScale;
                            double y = canvas.Bounds.Height  - height;

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
                    for (int i = 0; i < dataLength; i++)
                    {
                        // 计算当前数据点的频率值
                        double frequency = frequenciesLeft[i];

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

                //// 可选：添加坐标轴
                //var axis = new Line
                //{
                //    Stroke = Brushes.Black,
                //    StrokeThickness = 1,
                //    StartPoint = new Point(0, canvas.Bounds.Height - 5),
                //    EndPoint = new Point(canvas.Bounds.Width, canvas.Bounds.Height - 5)
                //};
                //canvas.Children.Add(axis);
            }
        }
    }
}

