using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringSurveyTools.CurveCalculate
{
    class Horizontal_JD
    {
        private List<HCurve_JD> curves;
        private struct Kinfo
        {
            internal double K0;//起始里程
            internal double K1;//结束里程
            internal int index0;//起始序号
            internal int index1;//结束序号
        }
        private List<Kinfo> Krange;
        public Horizontal_JD()
        {
            curves = new List<HCurve_JD>();
            Krange = new List<Kinfo>();
        }

        /// <summary>
        /// 添加平曲线交点信息
        /// </summary>
        /// <param name="curve">平曲线交点信息</param>
        public void Add(HCurve_JD curve)
        {
            if (curves.Count == 0)
            {
                curves.Add(curve);
                Kinfo k1 = new Kinfo();
                k1.K0 = curve.K - curve.Ls;
                k1.index0 = 0;
                k1.K1 = curve.K + curve.Ls;
                k1.index1 = 0;
                Krange.Add(k1);
            }
            else if 
                (curve.azimuth == curves[curves.Count - 1].azimuth+ curves[curves.Count - 1].alpha)
            {
                curves.Add(curve);
                Kinfo k1 = new Kinfo();
                k1.K0 = Krange[Krange.Count - 1].K0;
                k1.index0 = Krange[Krange.Count - 1].index0;
                k1.K1 = curve.K + curve.Ls;
                k1.index1 = curves.Count - 1;
                Krange.Add(k1);
            }
            else             //曲线不连续时
            {
                curves.Add(curve);
                Kinfo k1 = new Kinfo();
                k1.K0 = curve.K - curve.Ls;
                k1.index0 = curves.Count - 1;
                k1.K1 = curve.K + curve.Ls;
                k1.index1 = curves.Count - 1;
                Krange.Add(k1);
            }
        }
        /// <summary>
        /// 交点法平曲线转线元法平曲线
        /// </summary>
        /// <param name="LE">线元法平曲线</param>
        public void ToLE(out Horizontal_LE LE)
        {
            LE = new Horizontal_LE();
            HCurve_LE c0;
            HCurve_LE c1;
            HCurve_LE c2;
            HCurve_LE c3;
            for(int i=0;i<curves.Count;i++)
            {
                if (i == 0)
                    continue;
                if(curves[i].R==0)
                {
                    c0 = new HCurve_LE(0, curves[i - 1].Khz, curves[i].azimuth,
                        curves[i - 1].Xhz, curves[i - 1].Yhz, 0, curves[i].Kzh - curves[i - 1].Khz);
                    LE.Add(c0);
                    continue;
                }
                c0 = new HCurve_LE(0, curves[i - 1].Khz, curves[i].azimuth,
                        curves[i - 1].Xhz, curves[i - 1].Yhz, 0, curves[i].Kzh - curves[i - 1].Khz);
                c1 = new HCurve_LE(1, curves[i].Kzh, curves[i].azimuth, curves[i].Xzh,
                    curves[i].Yzh, curves[i].R, curves[i].Khy - curves[i].Kzh);
                c3 = new HCurve_LE(3, curves[i].Khy, curves[i].Ahy, curves[i].Xhy,
                    curves[i].Yhy, curves[i].R, curves[i].Kyh - curves[i].Khy);
                c2 = new HCurve_LE(2, curves[i].Kyh, curves[i].Ayh, curves[i].Xyh,
                    curves[i].Yyh, curves[i].R, curves[i].Khz - curves[i].Kyh);
                LE.Add(c0);
                LE.Add(c1);
                LE.Add(c3);
                LE.Add(c2);
            }
        }

    }
}
