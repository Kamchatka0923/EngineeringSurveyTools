using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 平曲线转换（交点法转线元法，线元法转交点法）
    /// </summary>
    public class HCurveTransformation
    {
        private Horizontal_LE LEs;
        private Horizontal_JD JDs;
        /// <summary>
        /// 构造函数
        /// </summary>
        public HCurveTransformation()
        {
            LEs = new Horizontal_LE();
            JDs = new Horizontal_JD();
        }
        /// <summary>
        /// 添加平曲线交点信息
        /// </summary>
        /// <param name="K">交点里程</param>
        /// <param name="X">交点X坐标</param>
        /// <param name="Y">交点Y坐标</param>
        /// <param name="A">方位角（弧度制）</param>
        /// <param name="alpha">转向角（弧度制）,左偏为负，右偏为正</param>
        /// <param name="R">曲线半径</param>
        /// <param name="Ls">缓和曲线长度</param>
        public void AddJD(double K, double X, double Y, double A, double alpha, double Ls, double R)
        {
            HCurve_JD hc = new HCurve_JD(K, X, Y, A, alpha, Ls, R);
            JDs.Add(hc);
            JDs.ToLE(out Horizontal_LE LE);
            LEs.curves.Clear();
            foreach (var hcle in LE.curves)
            {
                HCurve_LE le1 = new HCurve_LE(hcle.type, hcle.K0, hcle.azimuth, hcle.x0, hcle.y0, hcle.R, hcle.L);
                LEs.Add(le1);
            }
        }
        /// <summary>
        /// 输出线元法信息
        /// </summary>
        /// <returns></returns>
        public double[,] OutputLE()
        {
            double[,] LE = new double[LEs.curves.Count, 7];
            for(int i=0;i<LEs.curves.Count;i++)
            {
                LE[i, 0] = LEs.curves[i].type;
                LE[i, 1] = LEs.curves[i].K0;
                LE[i, 2] = LEs.curves[i].azimuth;
                LE[i, 3] = LEs.curves[i].x0;
                LE[i, 4] = LEs.curves[i].y0;
                LE[i, 5] = LEs.curves[i].R;
                LE[i, 6] = LEs.curves[i].L;
            }
            return LE;
        }

    }
}
