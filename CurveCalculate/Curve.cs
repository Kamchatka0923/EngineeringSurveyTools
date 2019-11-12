using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 道路曲线
    /// </summary>
    public class Curve
    {
        private Horizontal_LE Horizontal;
        private Horizontal_JD HorizontalJD;
        private Vertical_LE Vertical;
        private Milleage HK;
        private Milleage VK;

        public Curve()
        {
            Horizontal = new Horizontal_LE();
            HorizontalJD = new Horizontal_JD();
            Vertical = new Vertical_LE();
            HK = new Milleage();
            VK = new Milleage();
        }
        /// <summary>
        /// 添加平曲线线元
        /// </summary>
        /// <param name="type">线型（直线：0，直缓圆：1，圆缓直：2，圆弧：3）</param>
        /// <param name="K0">起始里程</param>
        /// <param name="A">起始方位角（弧度）</param>
        /// <param name="x0">起始北坐标</param>
        /// <param name="y0">起始东坐标</param>
        /// <param name="R">转弯半径</param>
        /// <param name="L">平曲线线元长</param>
        public void AddH_LE(int type, double K0, double A, double x0, double y0, double R, double L)
        {
            HK.AddHa(K0, L);
            K0 = HK.GetContinueK(K0);
            HCurve_LE hc = new HCurve_LE(type, K0, A, x0, y0, R, L);
            Horizontal.Add(hc);
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
        public void AddH_JD(double K, double X, double Y, double A, double alpha, double Ls, double R)
        {
            HCurve_JD hc = new HCurve_JD(K, X, Y, A, alpha, Ls, R);
            HorizontalJD.Add(hc);
            HorizontalJD.ToLE(out Horizontal_LE LE);
            Horizontal.curves.Clear();
            foreach(var hcle in LE.curves)
            {
                HK.AddHa(hcle.K0, hcle.K1, hcle.L);
                double K0 = HK.GetContinueK(hcle.K0);
                HCurve_LE le1 = new HCurve_LE(hcle.type, K0, hcle.azimuth, hcle.x0, hcle.y0, hcle.R, hcle.L);
                Horizontal.Add(le1);
            }
        }
        /// <summary>
        /// 添加断链
        /// </summary>
        /// <param name="type">断链点类型</param>
        /// <param name="bK">断链前里程，非负值</param>
        /// <param name="aK">断链后里程，非负值</param>
        /// <param name="pre">断链前冠号</param>
        /// <param name="suf">断链后冠号</param>
        /// <param name="name">断链点名</param>
        public void AddB(int type, double bK, double aK,
            string pre, string suf, string name = "")
        {
            BrokenChainage B = new BrokenChainage(type, bK, aK, pre, suf, name);
            HK.AddHB(B);
            VK.AddVB(B);
        }

        /// <summary>
        /// 添加竖曲线线元
        /// </summary>
        /// <param name="K">里程</param>
        /// <param name="H">高程</param>
        /// <param name="R">半径</param>
        public void AddV(double K, double H, double R)
        {
            VK.AddVa(K);
            K = VK.GetContinueK(K);
            VCurve_LE vc = new VCurve_LE(K, H, R);
            Vertical.Add(vc);
        }

        /// <summary>
        /// （里程、偏距、高差）计算（x,y,z）
        /// </summary>
        /// <param name="K">里程</param>
        /// <param name="O">偏距</param>
        /// <param name="H">高差</param>
        /// <param name="x">x</param>
        /// <param name="y">y</param>
        /// <param name="z">z</param>
        public void KOH2xyz(double K,double O,double H,out double x,out double y,out double z)
        {
            K = HK.GetContinueK(K);
            Horizontal.MilleageOffset2XY(K, O, out x, out y);
            Vertical.HightByMilleage(K, out double z0);
            z = H + z0;
        }
        /// <summary>
        /// （里程、偏距）计算（x,y）
        /// </summary>
        /// <param name="K">里程</param>
        /// <param name="O">偏距</param>
        /// <param name="x">x</param>
        /// <param name="y">y</param>
        public void KO2xy(double K, double O, out double x, out double y)
        {
            K = HK.GetContinueK(K);
            Horizontal.MilleageOffset2XY(K, O, out x, out y);
        }
        /// <summary>
        /// （x,y,z）计算（里程、偏距、高差）
        /// </summary>
        /// <param name="x">x</param>
        /// <param name="y">y</param>
        /// <param name="z">z</param>
        /// <param name="K">里程</param>
        /// <param name="O">偏距</param>
        /// <param name="H">高差</param>
        public void xyz2KOH(double x,double y,double z,out double K,out double O,out double H)
        {
            Horizontal.XY2MilleageOffset(x, y, out K, out O);
            Vertical.HightByMilleage(K, out double z0);
            H = z - z0;
            K = HK.GetActualK(K);
        }
        /// <summary>
        /// （x,y）计算（里程、偏距）
        /// </summary>
        /// <param name="x">x</param>
        /// <param name="y">y</param>
        /// <param name="K">里程</param>
        /// <param name="O">偏距</param>
        public void xy2KO(double x, double y, out double K, out double O)
        {
            Horizontal.XY2MilleageOffset(x, y, out K, out O);
            K = HK.GetActualK(K);
        }

    }
}
