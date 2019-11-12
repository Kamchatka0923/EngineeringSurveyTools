using System;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 交点法平曲线
    /// </summary>
    class HCurve_JD
    {
        public double K;    //交点里程
        public double X;    //交点X坐标
        public double Y;    //交点Y坐标
        public double azimuth;//方位角（弧度制）
        public double alpha;//转向角（弧度制）,左偏为负，右偏为正
        public double Ls;  //缓和曲线长度
        public double R;    //曲线半径

        public double Kzh;         //直缓点里程
        public double Khy;         //缓圆点里程
        public double Kyh;         //圆缓点里程
        public double Khz;         //缓直点里程

        public double Ahy;           //缓圆点切线角
        public double Ayh;           //圆缓点切线角

        public double Xzh;         //直缓点X坐标
        public double Yzh;         //直缓点Y坐标
        public double Xhy;         //缓圆点X坐标
        public double Yhy;         //缓圆点Y坐标
        public double Xyh;         //圆缓点X坐标
        public double Yyh;         //圆缓点Y坐标
        public double Xhz;         //缓直点X坐标
        public double Yhz;         //缓直点Y坐标
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="K">交点里程</param>
        /// <param name="X">交点X坐标</param>
        /// <param name="Y">交点Y坐标</param>
        /// <param name="A">方位角（弧度制）</param>
        /// <param name="alpha">转向角（弧度制）,左偏为负，右偏为正</param>
        /// <param name="R">曲线半径</param>
        /// <param name="Ls">缓和曲线长度</param>
        public HCurve_JD(double K,double X,double Y,double A,double alpha,double Ls,double R)
        {
            this.K = K;
            this.X = X;
            this.Y = Y;
            this.azimuth = A;
            this.alpha = alpha;
            this.R = R;
            this.Ls = Ls;
            if (alpha != 0)
            {
                CurvePara(Ls, R, alpha, out double m, out double P, out double betta0,
                out double TH, out double LH, out double EH, out double q);
                ZDK(K, TH, Ls, LH,
                    out Kzh, out Khy, out double Kqz, out Kyh, out Khz);
                double[,] ZD = CurveCordinate(Ls, R, alpha, azimuth, TH, X, Y);
                Xzh = ZD[0, 0];
                Yzh = ZD[0, 1];
                Xhy = ZD[1, 0];
                Yhy = ZD[1, 1];
                Xyh = ZD[2, 0];
                Yyh = ZD[2, 1];
                Xhz = ZD[3, 0];
                Yhz = ZD[3, 1];
                double ksi = (alpha > 0) ? 1 : -1;
                Ahy = azimuth + ksi * Ls / (2 * R);
                Ayh = azimuth + alpha - ksi * Ls / (2 * R);
            }
            else
            {
                Kzh = K;
                Khy = K;
                Kyh = K;
                Khz = K;
                Ahy = A;
                Ayh = A;
                Xzh = X;
                Xhy = X;
                Xyh = X;
                Xhz = X;
                Yzh = Y;
                Yhy = Y;
                Yyh = Y;
                Yhz = Y;

            }
        }
        /// <summary>
        /// 计算曲线参数
        /// </summary>
        /// <param name="Ls">缓和曲线长</param>
        /// <param name="R">圆曲线半径</param>
        /// <param name="alpha">转角（弧度制）</param>
        /// <param name="m">切垂距</param>
        /// <param name="P">圆曲线内移量</param>
        /// <param name="betta0">缓和曲线切线角</param>
        /// <param name="TH">切线长</param>
        /// <param name="LH">曲线长</param>
        /// <param name="EH">外矢距</param>
        /// <param name="q">切曲差</param>
        private void CurvePara(double Ls, double R, double alpha,
            out double m, out double P, out double betta0, out double TH, out double LH,
            out double EH, out double q)
        {
            alpha = Math.Abs(alpha);
            m = Ls / 2 - Math.Pow(Ls, 3) / (240 * R * R);
            P = Ls * Ls / (240 * R);
            betta0 = Ls / (2 * R);
            TH = m + (R + P) * Math.Tan(alpha / 2);
            LH = (alpha - 2 * betta0) * R + 2 * Ls;
            EH = (R + P) / Math.Cos(alpha / 2) - R;
            q = 2 * TH - LH;
        }
        /// <summary>
        /// 计算曲线主点里程
        /// </summary>
        /// <param name="Kjd">交点里程</param>
        /// <param name="TH">切线长</param>
        /// <param name="Ls">缓和曲线长</param>
        /// <param name="LH">曲线长</param>
        /// <param name="Kzh">直缓点里程</param>
        /// <param name="Khy">缓圆点里程</param>
        /// <param name="Kqz">曲中点里程</param>
        /// <param name="Kyh">圆缓点里程</param>
        /// <param name="Khz">缓直点里程</param>
        private void ZDK(double Kjd, double TH, double Ls, double LH,
            out double Kzh, out double Khy, out double Kqz, out double Kyh, out double Khz)
        {
            double LT = LH - Ls - Ls;
            Kzh = Kjd - TH;
            Khy = Kzh + Ls;
            Kqz = Kzh + LH / 2;
            Kyh = Khy + LT;
            Khz = Kyh + Ls;
        }
        /// <summary>
        /// 曲线主点坐标计算
        /// </summary>
        /// <param name="Ls">缓和曲线长</param>
        /// <param name="R">圆曲线半径</param>
        /// <param name="alpha">转向角（弧度制）,左偏为负，右偏为正</param>
        /// <param name="azi">方位角（弧度制）</param>
        /// <param name="TH">切线长</param>
        /// <param name="Xjd">交点x坐标</param>
        /// <param name="Yjd">交点y坐标</param>
        /// <returns>曲线主点坐标（ZH，HY，YH，HZ）</returns>
        private double[,] CurveCordinate(double Ls, double R, double alpha, double azi, double TH, double Xjd, double Yjd)
        {
            double[,] ZD = new double[4, 2];
            double Xzh, Xhy, Xyh, Xhz, Yzh, Yhy, Yyh, Yhz;
            Xzh = Xjd - TH * Math.Cos(azi);
            Yzh = Yjd - TH * Math.Sin(azi);
            Xhz = Xjd + TH * Math.Cos(azi + alpha);
            Yhz = Yjd + TH * Math.Sin(azi + alpha);

            Xhy = Ls - Math.Pow(Ls, 3) / (40 * R * R);
            Yhy = Ls * Ls / (6 * R);
            Xyh = -Ls - Math.Pow(-Ls, 3) / (40 * R * R);
            Yyh = Ls * Ls / (6 * R);

            ZD[0, 0] = Xzh;
            ZD[0, 1] = Yzh;
            ZD[3, 0] = Xhz;
            ZD[3, 1] = Yhz;
            if (alpha < 0)
            {
                alpha = Math.Abs(alpha);
                Yhy = -Yhy;
                Yyh = -Yyh;
            }
            ZD[1, 0] = Xzh + Xhy * Math.Cos(azi) - Yhy * Math.Sin(azi);
            ZD[1, 1] = Yzh + Xhy * Math.Sin(azi) + Yhy * Math.Cos(azi);
            ZD[2, 0] = Xhz + Xyh * Math.Cos(azi + alpha) + Yyh * Math.Sin(azi + alpha);
            ZD[2, 1] = Yhz + Xyh * Math.Sin(azi + alpha) - Yyh * Math.Cos(azi + alpha);

            return ZD;
        }
    }
}
