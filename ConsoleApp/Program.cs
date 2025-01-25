using AudioManager;
// See https://aka.ms/new-console-template for more information

foreach (var item in WASAPI.GetAllWASAPIDevices()) {
    Console.WriteLine(item);
}