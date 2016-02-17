using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace ControlPanel
{
    /// <summary>
    /// MainWindow.xaml에 대한 상호 작용 논리
    /// </summary>
    public partial class MainWindow : Window
    {
        NamedPipeServer PipeServer;
        bool start = false;

        public MainWindow()
        {
            InitializeComponent();

            // 1 이 버퍼에 쓰기
            PipeServer = new NamedPipeServer(@"\\.\pipe\NamePipeGoogolTech", 1);
            PipeServer.Start();
            start = true;
        }

        protected override void OnClosed(EventArgs e)
        {
            PipeServer.StopServer();
            base.OnClosed(e);
        }

        private void Solve_Click(object sender, RoutedEventArgs e)
        {
            textBox.Text = "Solve button click";
            if(start)
            {
                PipeServer.SendMessage("apply=", PipeServer.clientse);
            }
        }
    }
}
