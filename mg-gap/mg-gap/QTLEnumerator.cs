using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mg_gap
{
    class QTLEnumerator
    {
        private int chromosome = 0;
        private int startrange = 0;
        private int endrange = 0;
        private string gene = string.Empty;
        private string padj = string.Empty;
        private string qtlp = string.Empty;
        private string desc = string.Empty;

        public int Chromosome { get { return chromosome; } set { chromosome = value; } }
        public int StartRange { get { return startrange; } set { startrange = value; } }
        public int EndRange { get { return endrange; } set { endrange = value; } }
        public string Gene { get { return gene; } set { gene = value; } }
        public string P_Adj { get { return padj; } set { padj = value; } }
        public string QTL_P { get { return qtlp; } set { qtlp = value; } }
        public string Description { get { return desc; } set { desc = value; } }


    }
}
