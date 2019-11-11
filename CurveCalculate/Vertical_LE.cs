using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 线元法竖曲线线路
    /// </summary>
    class Vertical_LE
    {
        public List<VCurve_LE > curves;
        private double startL;
        private double endL;
        private double L;

        public Vertical_LE ()
        {
            curves = new List<VCurve_LE >();
        }

        public void Add(VCurve_LE curve)
        {
            curves.Add(curve);
            Init();
        }

        /// <summary>
        /// 里程参数值更新
        /// </summary>
        private void Flush()
        {
            startL = curves[0].K;
            endL = curves[curves.Count - 1].K;
            L = endL - startL;
        }
        /// <summary>
        /// 竖曲线要素初始化
        /// </summary>
        private void Init()
        {
            if (curves.Count <= 1)
                return;
            double i1, i2;
            i1 = (curves[1].H - curves[0].H) / (curves[1].K - curves[0].K);
            i2 = i1;
            curves[0].Add_i(i1, i2);
            for (int i = 1; i < curves.Count - 1; i++)
            {
                i1 = i2;
                i2 = (curves[i + 1].H - curves[i].H) / (curves[i + 1].K - curves[i].K);
                curves[i].Add_i(i1, i2);
            }
            i1 = i2;
            curves[curves.Count - 1].Add_i(i1, i2);
            Flush();
        }
        /// <summary>
        /// 根据里程搜索所在曲线
        /// </summary>
        /// <param name="L">里程</param>
        /// <returns>曲线顺序号</returns>
        private int SearchCurveByMilleage(double L)
        {
            if (L < startL || L > endL)
                return -1;
            int k = -2;
            double minL = this.L;
            for (int i = 1; i < curves.Count - 1; i++)
            {
                if (curves[i].i1 == curves[i].i2)
                    continue;
                if (Math.Abs(L - curves[i].K) < minL)
                {
                    minL = Math.Abs(L - curves[i].K);
                    k = i;
                }
            }
            return k;
        }
        /// <summary>
        /// 根据里程计算高程
        /// </summary>
        /// <param name="milleage">里程</param>
        /// <param name="H">高程</param>
        public void HightByMilleage(double milleage, out double H)
        {
            int k = SearchCurveByMilleage(milleage);
            if (k < 0)
                throw new Exception("Vertical:SearchCurveByMilleage:里程值不在竖曲线范围内");
            milleage += 0.00000001;
            H = curves[k].HightComputeByK(milleage);
        }
    }
}
