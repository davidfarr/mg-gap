using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mg_gap
{
    class Window
    {
        private int chromosome = 0;
        private List<SNP> snps_inwindow;
        private int window_id = 0;
        private int windowSize = 0;
        
        public int Chromosome { get { return chromosome; } set { chromosome = value; } }
        public List<SNP> SNP_Window { get { return snps_inwindow; } set { snps_inwindow = value; } }
        public int WindowID { get { return window_id; } set { window_id = value; } }
        public int S_Size { get { return windowSize; } set { windowSize = value; } }
    }
}
