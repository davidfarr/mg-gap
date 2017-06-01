using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mg_gap
{
    class FDR
    {
        public static List<string> assessment (string filepath)
        {
            List<string> treated_list = new List<string>();
            SortedList<string, double> sorted = new SortedList<string, double>();

            using (var fileStream = File.OpenRead(filepath))
            using (var streamReader = new StreamReader(fileStream))
            {
                String line;
                while ((line = streamReader.ReadLine()) != null)
                {
                    //pull the bs list into another list
                    string[] linearray = line.Split('\t');
                    sorted.Add(linearray[0], Convert.ToDouble(linearray[3]));
                }
            }

            //we have now a list sorted from the file based on auto-sorted p values
            foreach(var pair in sorted)
            {

            }



            return treated_list;
        }
    }
}
